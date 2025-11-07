#include "include/TandemTwister.hpp"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <cstdlib>
#include <cstdio>

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

const std::array<CommandInfo, 3> kCommandInfos = {{
    {"germline", "Genotyping tandem repeats from long-read alignments."},
    {"somatic",  "Somatic expansion profiling using long-read alignments."},
    {"assembly", "Genotyping tandem repeats from aligned assembly input."}
}};

std::string toLower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

void printDivider() {
    std::cerr << std::string(80, '-') << std::endl;
}

void printTitle(const std::string& title) {
    printDivider();
    std::cerr << title << std::endl;
    printDivider();
}

void printSectionHeader(const std::string& headline) {
    std::cerr << std::endl;
    std::cerr << headline << std::endl;
    std::cerr << std::string(headline.size(), '-') << std::endl;
}

void printOption(const std::string& option, const std::string& description) {
    std::cerr << "  " << std::left << std::setw(kOptionColumnWidth) << option << description << std::endl;
}

bool supportsAnsiStyling() {
    static const bool kSupported = []() {
        auto envTruthy = [](const char* name) -> bool {
            const char* value = std::getenv(name);
            if (!value || !*value) {
                return false;
            }
            return std::string(value) != "0";
        };

        if (envTruthy("TT_FORCE_COLOR") || envTruthy("TT_FORCE_BOLD") ||
            envTruthy("FORCE_COLOR") || envTruthy("CLICOLOR_FORCE")) {
            return true;
        }

        if (std::getenv("NO_COLOR") != nullptr) {
            return false;
        }
        if (!TT_ISATTY(TT_FILENO(stderr))) {
            return false;
        }
        const char* term = std::getenv("TERM");
        return term != nullptr && std::string(term) != "dumb";
    }();
    return kSupported;
}

std::string maybeBold(const std::string& text, bool emphasize) {
    if (!supportsAnsiStyling()) {
        return text;
    }
    if (emphasize) {
        return std::string("\033[1;97m") + text + "\033[0m";
    }
    return std::string("\033[1;37m") + text + "\033[0m";
}

std::string styleLabel(const std::string& label) {
    if (!supportsAnsiStyling()) {
        return label;
    }
    return std::string("\033[1;97m") + label + "\033[0m";
}

void printCommandsOverview(const std::string& activeCommand = "") {
    printSectionHeader("Commands");
    for (const auto& info : kCommandInfos) {
        const bool highlight = !activeCommand.empty() && activeCommand == info.name;
        printOption(maybeBold(info.name, highlight), info.description);
    }
}

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

}  // namespace

void printFullDoc(const std::string& activeCommand);
void printAssemblyDoc(bool emphasizeCommand = false);
void printGermlineSomaticDoc(const std::string& activeCommand = "");
void printGermlineDoc(const std::string& activeCommand = "");
void printLicenseNotice();

void TandemTwister::through_error(std::unordered_map<std::string, std::string>& args, const std::string& message) {
    std::cerr << message << std::endl;
    printUsage(args);
    exit(1);
}

bool TandemTwister::openFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return false;
    }
    file.close();
    return true;
}

void printVersion() {
    printTitle("TandemTwister");
    const std::string indent = "  ";
    std::cerr << indent << styleLabel("Purpose") << ": A tool for genotyping tandem repeats from long reads and aligned genome input" << std::endl;
    std::cerr << indent << styleLabel("Version") << ": 2.0.1" << std::endl;
    std::cerr << indent << styleLabel("Author") << ": Lion Ward Al Raei" << std::endl;
    std::cerr << indent << styleLabel("Contact") << ": Lionward.alraei@gmail.com" << std::endl;
    std::cerr << indent << styleLabel("Institute") << ": Max Planck Institute for Molecular Genetics" << std::endl;
    
    std::cerr << std::endl;
}

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

    int index = 1;
    std::string command;

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

    std::vector<std::string> accepted_arguments = { "-b","--bam","-r", "--ref",  "-m", "--motif_file", "-pad", "--padding", "-t", "--threads", "-h", "--help", "-f", "--fasta_file", "-o", "--output_file", "-s", "--output_file_statistics",
    "-mml", "-min_match_ratio_l", "-mms" , "-min_match_ration_s", "-crs", "-consensus_ratio_str", "-minPF", "--minPts_frac" ,"-crv", "--consensus_ratio_vntr", "-ms", "--match_score", "-mp", "--mismatch_penalty", "-gp", "--gap_penalty", "-seps", "--start_eps_str",
    "-sepv", "--start_eps_vntr", "-minPS", "--minPts_str", "-minPV", "--minPts_vntr", "-nls", "--noise_limit_str", "-nlv", "--noise_limit_vntr", "-qs", "--quality_score",
    "-ci", "--cluster_iter","-s","--sex", "--germline", "--somatic", "--assembly", "-rt", "--reads_type","-tanCon", "--tandem_run_threshold" , "-kpr", "--keepPhasingResults", "--sample", "-sn", "--keepCutReads","-kcr","-cor","--correct","-roz","--removeOutliersZscore", "-rtr", "--refineTrRegions","-minR", "--minReadsInRegion", "-btg", "--bamIsTagged", "--version", "-v", "--verbose", "-h" };

    for (int i = index; i < argc; ++i) {
        if (argv[i][0] == '-') {
            if (std::find(accepted_arguments.begin(), accepted_arguments.end(), argv[i]) != accepted_arguments.end()) {
                const std::string key = argv[i];

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

                if (i + 1 < argc && argv[i + 1][0] != '-' &&
                    std::find(accepted_arguments.begin(), accepted_arguments.end(), argv[i + 1]) == accepted_arguments.end()) {
                    args[key] = argv[i + 1];
                    ++i;
                } else {
                    args[key] = "";
                }
            } else {
                std::cerr << "Error: Unrecognized argument " << argv[i] << std::endl;
                printUsage(args);
                exit(1);
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

    std::string command_line;
    for (int i = 0; i < argc; i++) {
        command_line += argv[i];
        command_line += " ";
    }
    this->commandline = command_line;
    return args;
}

void TandemTwister::open_reference_file(const std::string& reference_file, faidx_t*& fai) {
    fai = fai_load(reference_file.c_str());
    if (fai == NULL) {
        throw std::runtime_error("Error: could not open reference file");
    }
}


void printFullDoc(const std::string& activeCommand) {
    printVersion();

    printSectionHeader("Usage");
    printOption("tandemtwister [global options] <command> [command options]", "Run TandemTwister in the selected analysis mode.");

    printCommandsOverview(activeCommand);

    const bool highlightGermline = activeCommand == "germline";
    const bool highlightSomatic = activeCommand == "somatic";
    const bool highlightAssembly = activeCommand == "assembly";

    printSectionHeader("Getting Started");
    printOption("tandemtwister " + maybeBold("germline", highlightGermline) + " --help", "Display germline-specific options.");
    printOption("tandemtwister " + maybeBold("somatic", highlightSomatic) + " --help", "Display somatic-specific options.");
    printOption("tandemtwister " + maybeBold("assembly", highlightAssembly) + " --help", "Display assembly-specific options.");

    printSectionHeader("Global Options");
    //printOption("-tanCon, --tandem_run_threshold", "Maximum bases for merging tandem repeat runs (default: 2 × motif size).");
    printOption("-v, --verbose", "Verbosity level (0: error, 1: critical, 2: info, 3: debug).");
    printOption("-h, --help | --version", "Display this help or version information.");

    printLicenseNotice();
}

void printAssemblyDoc(bool emphasizeCommand) {
    const bool highlight = emphasizeCommand;
    printSectionHeader("Usage");
    const std::string commandLine = "tandemtwister " + maybeBold("assembly", highlight) + " [options] -b <alignment_input> -m <motif_file> -r <reference_fasta> -o <output_file> -s <sample_sex>";
    printOption(commandLine,
                "Run assembly-based analysis on aligned genome input.");

    printSectionHeader("Required Arguments");
    printOption("Command: " + maybeBold("assembly", highlight) + " [--assembly]", "Selects the assembly analysis workflow (legacy flag supported).");
    printOption("-b, --bam", "Path to the aligned assembly BAM file.");
    printOption("-r, --ref", "Reference FASTA file (.fa / .fna).");
    printOption("-m, --motif_file", "Motif definition file (BED / TSV / CSV).");
    printOption("-o, --output_file", "Output file for region, motif, and haplotype copy numbers.");
    printOption("-s, --sex", "Sample sex (0 | female, 1 | male).");
    printOption("-sn, --sample", "Optional sample identifier.");
    printOption("-rt, --reads_type", "Sequencing platform / read type (default: CCS).");

    printSectionHeader("Alignment Parameters");
    printOption("-mml, --min_match_ratio_l", "Minimum match ratio for long motifs (default: 0.5).");
    // printOption("-ms, --match_score", "Match score (default: 10).");
    // printOption("-mp, --mismatch_penalty", "Mismatch penalty (default: -2).");
    // printOption("-gp, --gap_penalty", "Gap penalty (default: -5).");

    printSectionHeader("General");
    printOption("-tanCon, --tandem_run_threshold", "Maximum bases for merging tandem repeat runs (default: 2 × motif size).");
    printOption("-v, --verbose", "Verbosity level (0: error, 1: critical, 2: info, 3: debug).");
    printOption("-h, --help | --version", "Display this help or version information.");
}

void printGermlineSomaticDoc(const std::string& activeCommand) {
    const bool highlightGermline = activeCommand == "germline";
    const bool highlightSomatic = activeCommand == "somatic";

    printSectionHeader("Required Arguments");
    printOption("Command: " + maybeBold("germline", highlightGermline) + " | " + maybeBold("somatic", highlightSomatic) + " [--germline | --somatic]  ", "Selects the germline or somatic analysis workflow (legacy flags supported).");
    printOption("-b, --bam", "Path to the BAM file containing aligned reads.");
    printOption("-r, --ref", "Reference FASTA file (.fa / .fna).");
    printOption("-m, --motif_file", "Motif definition file (BED / TSV / CSV).");
    printOption("-o, --output_file", "Output file for region, motif, and haplotype copy numbers.");
    printOption("-s, --sex", "Sample sex (0 | female, 1 | male).");
    printOption("-sn, --sample", "Optional sample identifier.");
    printOption("-rt, --reads_type", "Sequencing platform / read type (default: CCS).");

    printSectionHeader("Alignment Parameters");
    printOption("-mml, --min_match_ratio_l", "Minimum match ratio for long motifs (default: 0.5).");
    // printOption("-ms, --match_score", "Match score (default: 10).");
    // printOption("-mp, --mismatch_penalty", "Mismatch penalty (default: -2).");
    // printOption("-gp, --gap_penalty", "Gap penalty (default: -5).");

    printSectionHeader("Read Extraction");
    printOption("-s, --output_file_statistics", "Optional phasing summary output file.");
    printOption("-pad, --padding", "Padding around the STR region when extracting reads (default: 0).");
    printOption("-t, --threads", "Number of threads (default: 1).");
    printOption("-kpr, --keepPhasingResults", "Persist intermediate phasing results (default: false).");
    printOption("-kcr, --keepCutReads", "Retain reads trimmed during preprocessing (default: false).");
    printOption("-minR, --minReadsInRegion", "Minimum spanning reads required per region (default: 2).");
    printOption("-btg, --bamIsTagged", "Treat BAM as pre-tagged/phased (default: false).");
    printOption("-qs, --quality_score", "Minimum read quality score (default: 10, max: 60).");

    printSectionHeader("Correction");
    printOption("-cor, --correct", "Enable genotype correction (CCS default: false; CLR/ONT default: true).");
    printOption("-crs, --consensus_ratio_str", "Minimum consensus ratio for STR correction (default: 0.3).");
    printOption("-crv, --consensus_ratio_vntr", "Minimum consensus ratio for VNTR correction (default: 0.3).");
    printOption("-roz, --removeOutliersZscore", "Remove outliers before phasing (default: false).");
    // printOption("-rtr, --refineTrRegions", "Refine tandem repeat regions (default: true).");

    printSectionHeader("Clustering");
    printOption("-seps, --start_eps_str", "Initial epsilon for STR clustering (default: 0.2).");
    printOption("-sepv, --start_eps_vntr", "Initial epsilon for VNTR clustering (default: 0.2).");
    printOption("-minPF, --minPts_frac", "Minimum read fraction per cluster (default: 0.12).");
    printOption("-nls, --noise_limit_str", "Noise limit for STR clustering (default: 0.2).");
    printOption("-nlv, --noise_limit_vntr", "Noise limit for VNTR clustering (default: 0.35).");
    printOption("-ci, --cluster_iter", "Number of clustering iterations (default: 40).");

    printSectionHeader("General");
    printOption("-v, --verbose", "Verbosity level (0: error, 1: critical, 2: info, 3: debug).");
    printOption("-h, --help | --version", "Display this help or version information.");
}

void printGermlineDoc(const std::string& activeCommand) {
    printCommandsOverview(activeCommand);
    printSectionHeader("Usage");
    const bool highlightGermline = activeCommand == "germline";
    const bool highlightSomatic = activeCommand == "somatic";
    printOption("tandemtwister " + maybeBold("germline", highlightGermline) + " [options] -b <input_bam_file> -m <motif_file> -r <reference_fasta> -o <output_file> -s <sample_sex>",
                "Run germline analysis on long-read alignments.");
    printOption("tandemtwister " + maybeBold("somatic", highlightSomatic) + "  [options] -b <input_bam_file> -m <motif_file> -r <reference_fasta> -o <output_file> -s <sample_sex>",
                "Run somatic analysis on long-read alignments.");
    printGermlineSomaticDoc(activeCommand);
}

void printLicenseNotice() {
    printSectionHeader("License & Warranty");
    std::cerr << "  " << styleLabel("License") << ": BSD 3-Clause" << std::endl;
    std::cerr << "  " << styleLabel("Warranty") << ": Provided \"AS IS\" without warranty of any kind." << std::endl;
    std::cerr << "  " << styleLabel("Usage") << ": For research purposes only; not for diagnostic or clinical use." << std::endl;
}

void TandemTwister::printUsage(std::unordered_map<std::string, std::string> &args) {
    const auto commandIt = args.find("command");

    if (commandIt != args.end()) {
        const std::string& command = commandIt->second;
        printVersion();
        if (command == "germline") {
            printSectionHeader("Germline Analysis");
            std::cerr << "Genotyping tandem repeats from long reads." << std::endl;
            printGermlineDoc("germline");
            printLicenseNotice();
            return;
        }
        if (command == "somatic") {
            printSectionHeader("Somatic Analysis");
            std::cerr << "Genotyping tandem repeats from long reads." << std::endl;
            printGermlineDoc("somatic");
            printLicenseNotice();
            return;
        }
        if (command == "assembly") {
            printSectionHeader("Assembly Analysis");
            std::cerr << "Genotyping tandem repeats from aligned genome input." << std::endl;
            printAssemblyDoc(true);
            printLicenseNotice();
            return;
        }
    }
    
    if (args.find("--germline") != args.end() || args.find("--somatic") != args.end()) {
        printVersion();
        if (args.find("--germline") != args.end()) {
            printSectionHeader("Germline Analysis");
            std::cerr << "Genotyping tandem repeats from long reads." << std::endl;
        } else {
            printSectionHeader("Somatic Analysis");
            std::cerr << "Genotyping tandem repeats from long reads." << std::endl;
        }

        printGermlineDoc(args.find("--germline") != args.end() ? "germline" : "somatic");
        printLicenseNotice();
        return;
    }

    if (args.find("--assembly") != args.end()) {
        printVersion();
        printSectionHeader("Assembly Analysis");
        std::cerr << "Genotyping tandem repeats from aligned genome input." << std::endl;
        printAssemblyDoc(true);
        printLicenseNotice();
        return;
    }

    printFullDoc(commandIt != args.end() ? commandIt->second : std::string());
}

