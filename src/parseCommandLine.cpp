#include "include/TandemTwister.hpp"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <cstdlib>
#include <cstdio>

// for ANSI styling
#ifdef _WIN32
#include <io.h>
#define TT_ISATTY _isatty
#define TT_FILENO _fileno
#else
#include <unistd.h>
#define TT_ISATTY isatty
#define TT_FILENO fileno
#endif

namespace {

constexpr int kOptionColumnWidth = 48;

struct CommandInfo {
    const char* name;
    const char* description;
};

struct OptionInfo {
    std::vector<std::string> flags;  // e.g., {"-b", "--bam"}
    const char* description;
    const char* section;  // e.g., "Required Arguments", "Alignment Parameters"
    bool required;
    bool command_specific;  // true if only for specific commands
    std::vector<std::string> commands;  // empty = all commands, or {"germline", "somatic"}
};

// Command information
const std::array<CommandInfo, 3> kCommandInfos = {{
    {"germline", "Genotyping tandem repeats from long-read alignments."},
    {"somatic",  "Somatic expansion profiling using long-read alignments."},
    {"assembly", "Genotyping tandem repeats from aligned assembly input."}
}};

// Data-driven option definitions - single source of truth
const std::vector<OptionInfo> kAllOptions = {
    // Commands
    {{"germline", "--germline"}, "Selects the germline analysis workflow", "Commands", false, false, {}},
    {{"somatic", "--somatic"}, "Selects the somatic analysis workflow", "Commands", false, false, {}},
    {{"assembly", "--assembly"}, "Selects the assembly analysis workflow", "Commands", false, false, {}},
    
    // Required Arguments (common)
    {{"-b", "--bam"}, "Path to the BAM file containing aligned reads", "Required Arguments", true, false, {}},
    {{"-r", "--ref"}, "Reference FASTA file (.fa / .fna)", "Required Arguments", true, false, {}},
    {{"-m", "--motif_file"}, "Motif definition file (BED / TSV / CSV)", "Required Arguments", true, false, {}},
    {{"-o", "--output_file"}, "Output file for region, motif, and haplotype copy numbers", "Required Arguments", true, false, {}},
    {{"-s", "--sex"}, "Sample sex (0 | female, 1 | male)", "Required Arguments", true, false, {}},
    {{"-sn", "--sample"}, "Optional sample identifier", "Required Arguments", false, false, {}},
    {{"-rt", "--reads_type"}, "Sequencing platform / read type (default: CCS)", "Required Arguments", false, false, {}},
    
    // Required Arguments (germline/somatic specific)
    {{"-t", "--threads"}, "Number of threads (default: 1)", "Required Arguments", false, true, {"germline", "somatic"}},
    
    // Alignment Parameters
    {{"-mml", "--min_match_ratio_l"}, "Minimum match ratio for long motifs (default: 0.5)", "Alignment Parameters", false, false, {}},
    
    // Read Extraction (germline/somatic only) - Note: -s is used for --sex, so statistics uses long form only
    {{"--output_file_statistics"}, "Optional phasing summary output file", "Read Extraction", false, true, {"germline", "somatic"}},
    {{"-pad", "--padding"}, "Padding around the STR region when extracting reads (default: 0)", "Read Extraction", false, true, {"germline", "somatic"}},
    {{"-kcr", "--keepCutReads"}, "Retain reads trimmed during preprocessing (default: false)", "Read Extraction", false, true, {"germline", "somatic"}},
    {{"-minR", "--minReadsInRegion"}, "Minimum spanning reads required per region (default: 2)", "Read Extraction", false, true, {"germline", "somatic"}},
    {{"-btg", "--bamIsTagged"}, "Treat BAM as pre-tagged/phased (default: false)", "Read Extraction", false, true, {"germline", "somatic"}},
    {{"-qs", "--quality_score"}, "Minimum read quality score (default: 10, max: 60)", "Read Extraction", false, true, {"germline", "somatic"}},
    
    // Reference-based Correction (germline/somatic only)
    {{"-rtr", "--refineTrRegions"}, "Refine tandem repeat regions (default: false)", "Reference-based Correction Parameters", false, true, {"germline", "somatic"}},
    {{"-tanCon", "--tandem_run_threshold"}, "Maximum bases for merging tandem-repeat runs (default: 2 × motif size)", "Reference-based Correction Parameters", false, true, {"germline", "somatic"}},
    
    // Read-based Correction (germline/somatic only)
    {{"-cor", "--correct"}, "Enable genotype correction (CCS default: false; CLR/ONT default: true)", "Read-based Correction Parameters", false, true, {"germline", "somatic"}},
    {{"-crs", "--consensus_ratio_str"}, "Minimum consensus ratio for STR correction (default: 0.3)", "Read-based Correction Parameters", false, true, {"germline", "somatic"}},
    {{"-crv", "--consensus_ratio_vntr"}, "Minimum consensus ratio for VNTR correction (default: 0.3)", "Read-based Correction Parameters", false, true, {"germline", "somatic"}},
    {{"-roz", "--removeOutliersZscore"}, "Remove outliers before phasing (default: false)", "Read-based Correction Parameters", false, true, {"germline", "somatic"}},
    
    // Clustering (germline/somatic only)
    {{"-seps", "--start_eps_str"}, "Initial epsilon for STR clustering (default: 0.2)", "Clustering", false, true, {"germline", "somatic"}},
    {{"-sepv", "--start_eps_vntr"}, "Initial epsilon for VNTR clustering (default: 0.2)", "Clustering", false, true, {"germline", "somatic"}},
    {{"-minPF", "--minPts_frac"}, "Minimum read fraction per cluster (default: 0.12)", "Clustering", false, true, {"germline", "somatic"}},
    {{"-nls", "--noise_limit_str"}, "Noise limit for STR clustering (default: 0.2)", "Clustering", false, true, {"germline", "somatic"}},
    {{"-nlv", "--noise_limit_vntr"}, "Noise limit for VNTR clustering (default: 0.35)", "Clustering", false, true, {"germline", "somatic"}},
    {{"-ci", "--cluster_iter"}, "Number of clustering iterations (default: 40)", "Clustering", false, true, {"germline", "somatic"}},
    
    // General
    {{"-v", "--verbose"}, "Verbosity level (0: error, 1: critical, 2: info, 3: debug)", "General", false, false, {}},
    {{"-h", "--help"}, "Display this help or version information", "General", false, false, {}},
    {{"--version"}, "Display version information", "General", false, false, {}},
};

// Convert string to lowercase
std::string toLower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), 
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

// Print divider
void printDivider() {
    std::cerr << std::string(80, '-') << std::endl;
}

// Print title
void printTitle(const std::string& title) {
    printDivider();
    std::cerr << title << std::endl;
    printDivider();
}

// Print section header
void printSectionHeader(const std::string& headline) {
    std::cerr << std::endl;
    std::cerr << headline << std::endl;
    std::cerr << std::string(headline.size(), '-') << std::endl;
}

// Print option
void printOption(const std::string& option, const std::string& description) {
    std::cerr << "  " << std::left << std::setw(kOptionColumnWidth) << option << description << std::endl;
}

// Check if ANSI styling is supported
bool supportsAnsiStyling() {
    static const bool kSupported = []() {
        auto envTruthy = [](const char* name) -> bool {
            const char* value = std::getenv(name);
            if (!value || !*value) return false;
            return std::string(value) != "0";
        };

        if (envTruthy("TT_FORCE_COLOR") || envTruthy("TT_FORCE_BOLD") ||
            envTruthy("FORCE_COLOR") || envTruthy("CLICOLOR_FORCE")) {
            return true;
        }

        if (std::getenv("NO_COLOR") != nullptr) return false;
        if (!TT_ISATTY(TT_FILENO(stderr))) return false;
        const char* term = std::getenv("TERM");
        return term != nullptr && std::string(term) != "dumb";
    }();
    return kSupported;
}

// Maybe bold text
std::string maybeBold(const std::string& text, bool emphasize) {
    if (!supportsAnsiStyling()) return text;
    if (emphasize) {
        return std::string("\033[1;97m") + text + "\033[0m";
    }
    return std::string("\033[1;37m") + text + "\033[0m";
}

// Style label
std::string styleLabel(const std::string& label) {
    if (!supportsAnsiStyling()) return label;
    return std::string("\033[1;97m") + label + "\033[0m";
}

// Print version
void printVersion() {
    printTitle("TandemTwister");
    const std::string indent = "  ";
    std::cerr << indent << styleLabel("Purpose") << ": A tool for genotyping tandem repeats from long reads and aligned genome input" << std::endl;
    std::cerr << indent << styleLabel("Version") << ": 0.1.0" << std::endl;
    std::cerr << indent << styleLabel("Author") << ": Lion Ward Al Raei" << std::endl;
    std::cerr << indent << styleLabel("Contact") << ": Lionward.alraei@gmail.com" << std::endl;
    std::cerr << indent << styleLabel("Institute") << ": Max Planck Institute for Molecular Genetics" << std::endl;
    std::cerr << std::endl;
}

// Print license notice
void printLicenseNotice() {
    printSectionHeader("License & Warranty");
    std::cerr << "  " << styleLabel("License") << ": BSD 3-Clause" << std::endl; 
    std::cerr << "  " << styleLabel("Warranty") << ": Provided \"AS IS\" without warranty of any kind." << std::endl;
    std::cerr << "  " << styleLabel("Usage") << ": For research purposes only; not for diagnostic or clinical use." << std::endl;
}

// Print commands overview
void printCommandsOverview(const std::string& activeCommand = "") {
    printSectionHeader("Commands");
    for (const auto& info : kCommandInfos) {
        const bool highlight = !activeCommand.empty() && activeCommand == info.name;
        printOption(maybeBold(info.name, highlight), info.description);
    }
}

// Get options for a specific command and section
std::vector<const OptionInfo*> getOptionsForCommand(const std::string& command, const std::string& section) {
    std::vector<const OptionInfo*> result;
    for (const auto& opt : kAllOptions) {
        // Skip if section doesn't match
        if (opt.section != section) continue;
        
        // Check if option applies to this command
        if (opt.command_specific) {
            if (opt.commands.empty() || 
                std::find(opt.commands.begin(), opt.commands.end(), command) == opt.commands.end()) {
                continue;
            }
        }
        result.push_back(&opt);
    }
    return result;
}

// Print options for a section
void printOptionsForSection(const std::string& command, const std::string& section, 
                           const std::string& activeCommand = "") {
    auto options = getOptionsForCommand(command, section);
    if (options.empty()) return;
    
    printSectionHeader(section);
    
    for (const auto* opt : options) {
        // Build option string (e.g., "-b, --bam")
        std::string optStr;
        for (size_t i = 0; i < opt->flags.size(); ++i) {
            if (i > 0) optStr += ", ";
            optStr += opt->flags[i];
        }
        
        // Special handling for command options
        if (section == "Commands") {
            const bool highlight = !activeCommand.empty() && 
                std::find(opt->flags.begin(), opt->flags.end(), "--" + activeCommand) != opt->flags.end();
            optStr = maybeBold(opt->flags[0], highlight);
            if (opt->flags.size() > 1) {
                optStr += " [" + opt->flags[1] + "]";
            }
        }
        
        printOption(optStr, opt->description);
    }
}

// Print command-specific documentation
void printCommandDoc(const std::string& command, const std::string& activeCommand = "") {
    const bool highlight = activeCommand == command || activeCommand.empty();
    
    // Usage
    printSectionHeader("Usage");
    std::string usage = "tandemtwister " + maybeBold(command, highlight) + 
                       " [options] -b <input_bam> -m <motif_file> -r <reference_fasta> -o <output_file> -s <sample_sex>";
    printOption(usage, "");
    
    // Define sections in order
    std::vector<std::string> sections = {
        "Required Arguments",
        "Alignment Parameters",
        "Read Extraction",
        "Reference-based Correction Parameters",
        "Read-based Correction Parameters",
        "Clustering",
        "General"
    };
    
    // Print each section
    for (const auto& section : sections) {
        printOptionsForSection(command, section, activeCommand);
    }
}

// Print full documentation
void printFullDoc(const std::string& activeCommand = "") {
    
    printSectionHeader("Usage");
    printOption("tandemtwister [global options] <command> [command options]", "");
    
    printCommandsOverview(activeCommand);
    
    printSectionHeader("Getting Started");
    for (const auto& cmd : kCommandInfos) {
        const bool highlight = activeCommand == cmd.name;
        printOption("tandemtwister " + maybeBold(cmd.name, highlight) + " --help", 
                   "Display " + std::string(cmd.name) + "-specific options.");
    }
    
    printOptionsForSection("", "General", activeCommand);
    printLicenseNotice();
}

}  // namespace

// Set command
void setCommand(std::unordered_map<std::string, std::string>& args,
                const std::string& command,
                const std::unordered_map<std::string, std::string>& commandFlags) {
    args["command"] = command;
    for (const auto& entry : commandFlags) {
        if (entry.second == command) {
            args.emplace(entry.first, "");
        }
    }
}

// Through error
void TandemTwister::through_error(std::unordered_map<std::string, std::string>& args, const std::string& message) {
    std::cerr << message << std::endl;
    printUsage(args);
    exit(1);
}

// Open file
bool TandemTwister::openFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return false;
    }
    file.close();
    return true;
}

// Parse command line - simplified version
std::unordered_map<std::string, std::string> TandemTwister::parseCommandLine(int argc, char** argv) {
    std::unordered_map<std::string, std::string> args;
    if (argc <= 1) {
        printUsage(args);
        exit(1);
    }

    const std::unordered_set<std::string> commandNames = {"germline", "somatic", "assembly"};
    const std::unordered_map<std::string, std::string> commandFlags = {
        {"--germline", "germline"},
        {"--somatic", "somatic"},
        {"--assembly", "assembly"}
    };

    // Build accepted arguments set from kAllOptions
    std::unordered_set<std::string> accepted_arguments = {"-h", "--help", "--version", "-v", "--verbose"};
    for (const auto& opt : kAllOptions) {
        for (const auto& flag : opt.flags) {
            accepted_arguments.insert(flag);
        }
    }

    int index = 1;
    std::string command;

    // Parse command
    if (index < argc && argv[index][0] != '-') {
        const std::string candidate = toLower(argv[index]);
        if (commandNames.count(candidate) != 0) {
            command = candidate;
            setCommand(args, command, commandFlags);
            ++index;
        } else {
            std::cerr << "Error: Unrecognized command " << argv[index] << std::endl;
            printFullDoc(command);
            exit(1);
        }
    }

    // Parse arguments
    for (int i = index; i < argc; ++i) {
        if (argv[i][0] == '-') {
            if (accepted_arguments.count(argv[i]) == 0) {
                std::cerr << "Error: Unrecognized argument " << argv[i] << std::endl;
                printUsage(args);
                exit(1);
            }
            
            const std::string key = argv[i];
            
            // Check if it's a command flag
            const auto flagIt = commandFlags.find(key);
            if (flagIt != commandFlags.end()) {
                if (!command.empty() && command != flagIt->second) {
                    std::cerr << "Error: Multiple commands specified (" << command << " and " << flagIt->second << ")." << std::endl;
                    printUsage(args);
                    exit(1);
                }
                command = flagIt->second;
                setCommand(args, command, commandFlags);
            }

            // Get value if available
            if (i + 1 < argc && argv[i + 1][0] != '-' &&
                accepted_arguments.count(argv[i + 1]) == 0) {
                args[key] = argv[i + 1];
                ++i;
            } else {
                args[key] = "";
            }
        } else {
            const std::string token = toLower(argv[i]);
            if (commandNames.count(token) != 0) {
                if (!command.empty() && command != token) {
                    std::cerr << "Error: Multiple commands specified (" << command << " and " << token << ")." << std::endl;
                    printUsage(args);
                    exit(1);
                }
                command = token;
                setCommand(args, command, commandFlags);
            } else {
                std::cerr << "Error: Unrecognized argument " << argv[i] << std::endl;
                printUsage(args);
                exit(1);
            }
        }
    }

    // Handle missing command
    if (command.empty()) {
        for (const auto& entry : commandFlags) {
            if (args.find(entry.first) != args.end()) {
                command = entry.second;
                setCommand(args, command, commandFlags);
                break;
            }
        }
    }

    const bool hasHelp = args.find("-h") != args.end() || args.find("--help") != args.end();
    const bool hasVersion = args.find("--version") != args.end();

    if (hasVersion && command.empty()) {
        printVersion();
        exit(0);
    }

    if (hasHelp && command.empty()) {
        printUsage(args);
        exit(0);
    }

    if (command.empty()) {
        std::cerr << "Error: Missing analysis command. Use one of: germline, somatic, assembly." << std::endl;
        printFullDoc(command);
        exit(1);
    }

    args["command"] = command;

    if (hasVersion) {
        printVersion();
        exit(0);
    }

    if (hasHelp) {
        printUsage(args);
        exit(0);
    }

    // Build command line string
    std::string command_line;
    for (int i = 0; i < argc; i++) {
        command_line += argv[i];
        command_line += " ";
    }
    this->commandline = command_line;
    return args;
}

// Print usage - simplified and data-driven
void TandemTwister::printUsage(std::unordered_map<std::string, std::string> &args) {
    const auto commandIt = args.find("command");
    std::string command = (commandIt != args.end()) ? commandIt->second : "";

    // Handle legacy flags
    if (command.empty()) {
        if (args.find("--germline") != args.end()) command = "germline";
        else if (args.find("--somatic") != args.end()) command = "somatic";
        else if (args.find("--assembly") != args.end()) command = "assembly";
    }

    printVersion();

    if (command == "germline") {
        printSectionHeader("Germline Analysis");
        std::cerr << "Genotyping tandem repeats from long reads." << std::endl;
        printCommandsOverview("germline");
        printCommandDoc("germline", "germline");
        printLicenseNotice();
    } else if (command == "somatic") {
        printSectionHeader("Somatic Analysis");
        std::cerr << "Genotyping tandem repeats from long reads." << std::endl;
        printCommandsOverview("somatic");
        printCommandDoc("somatic", "somatic");
        printLicenseNotice();
    } else if (command == "assembly") {
        printSectionHeader("Assembly Analysis");
        std::cerr << "Genotyping tandem repeats from aligned genome input." << std::endl;
        printCommandDoc("assembly", "assembly");
        printLicenseNotice();
    } else {
        printFullDoc(command);
    }
}
