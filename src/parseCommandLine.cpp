#include "include/TandemTwister.hpp"

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

void printVersion(){
    std::cerr << "TandemTwister: A tool for genotyping tandem repeats from long reads and aligned genome input" << std::endl;
    std::cerr << "Version: 2.0.1" << std::endl;
    std::cerr << "Author: Lion Ward Al Raei " << std::endl;
    std::cerr << "Email: Lionward.alraei@gmail.com " << std::endl;
    std::cerr << "Institue:  (Max Planck institute for molecular genetics)" << std::endl;
    std::cerr << std::endl;
}

std::unordered_map<std::string, std::string> TandemTwister::parseCommandLine(int argc, char** argv) {


    /*
    * Parse the command line arguments

    * @param argc: The number of arguments
    * @param argv: The arguments
    * @return: A map containing the arguments
    * 
    
    */
    std::unordered_map<std::string, std::string> args;
    if (argc == 1) {
        printUsage(args);
        exit(1);
    }

    std::vector<std::string> accepted_arguments = { "-b","--bam","-r", "--ref",  "-m", "--motif_file", "-pad", "--padding", "-t", "--threads", "-h", "--help", "-f", "--fasta_file", "-o", "--output_file", "-s", "--output_file_statistics", 
    "-mml", "-min_match_ratio_l", "-mms" , "-min_match_ration_s", "-crs", "-consensus_ratio_str", "-minPF", "--minPts_frac" ,"-crv", "--consensus_ratio_vntr", "-ms", "--match_score", "-mp", "--mismatch_penalty", "-gp", "--gap_penalty", "-seps", "--start_eps_str",
    "-sepv", "--start_eps_vntr", "-minPS", "--minPts_str", "-minPV", "--minPts_vntr", "-nls", "--noise_limit_str", "-nlv", "--noise_limit_vntr", "-qs", "--quality_score",
    "-ci", "--cluster_iter","-s","--sex", "--germline", "--somatic", "--assembly", "-rt", "--reads_type","-tanCon", "--tandem_run_threshold" , "-kpr", "--keepPhasingResults", "--sample", "-sn", "--keepCutReads","-kcr","-cor","--correct","-roz","--removeOutliersZscore", "-rtr", "--refineTrRegions","-minR", "--minReadsInRegion", "-btg", "--bamIsTagged", "--version", "-v", "--verbose", "-h" };
    for (int i = 1; i < argc; i++){
        // check if the argument starts with a dash
        if (argv[i][0] == '-') {
            // check if the argument is in the accepted arguments
            if (std::find(accepted_arguments.begin(), accepted_arguments.end(), argv[i]) != accepted_arguments.end()) {
                // check if the argument is not the last one
                if (i + 1 < argc) {
                    // check if the next argument does not start with a dash
                    if (std::find(accepted_arguments.begin(), accepted_arguments.end(), argv[i + 1]) == accepted_arguments.end()) {

                        args[argv[i]] = argv[i + 1];
                        i++;
                    }
                    else {
                        args[argv[i]] = "";
                    }
                }
                else {
                    args[argv[i]] = "";
                }
            }
            else {
                std::cerr << "Error: Unrecognized argument " << argv[i] << std::endl;
                printUsage(args);
                exit(1);
            }
        }
        else {
            std::cerr << "Error: Unrecognized argument " << argv[i] << std::endl;
            printUsage(args);
            exit(1);
        }
    }
    // if arguments only 1 
    if (args.size() == 1) {
        if (args.find("-h") != args.end() || args.find("--help") != args.end() || args.find("--assembly") != args.end() || args.find("--germline") != args.end() || args.find("--somatic") != args.end() ) {
            printUsage(args);
            exit(0);
        }
        if (args.find("--version") != args.end()) {
            printVersion();
            exit(0);
        }
    }
    // put the whole command line in string 
    std::string command_line = "";
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


void printFullDoc(){

    printVersion();
    std::cerr << "TandemTwister: A tool for genotyping tandem repeats from long reads and aligned genome input" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage: TandemTwister --germline [options] -b <input_bam_file> -m <motif_file> -r <input_reference_file> -o <output_file>  -s <sample_sex>" << std::endl;
    std::cerr << "Usage: TandemTwister --somatic [options] -b <input_bam_file> -m <motif_file> -r <input_reference_file> -o <output_file>  -s <sample_sex>" << std::endl;
    std::cerr << "Usage: TandemTwister --assembly [options] -b <input_bam_file> -m <motif_file> -r <input_reference_file> -o <output_file>  -s <sample_sex>" << std::endl;

    std::cerr << std::endl;
    std::cerr << "Required input:" << std::endl;
    std::cerr << "--germline/--somatic/--assembly           Type of analysis to perform" << std::endl;
    std::cerr << " -b, --bam                                path to the bam file of the aligned reads to the reference genome" << std::endl;
    std::cerr << " -r, --ref                                Path to the input reference file [.fa/.fna]." << std::endl;
    std::cerr << " -m, --motif_file                         path to the file with reference coordinates and motif sequence [BED/TSV/CSV]" << std::endl;
    std::cerr << " -o, --output_file                        Output file containing region, motif, hap1 and hap2 copy numbers" << std::endl;
    std::cerr << " -s, --sex                                Sample sex (0|female 1|male)" << std::endl;
    std::cerr << " -sn, --sample                            name of sample" << std::endl;
    std::cerr << " -rt, --reads_type                        Type of reads (Default: CCS)" << std::endl;

    std::cerr << "[options]: Dynamic programming alignment Parameters: " << std::endl;
    std::cerr << " -mml, --min_match_ratio_l                Minimum match ratio for long motifs (Default: 0.5)" << std::endl;
    std::cerr << " -mms, --min_match_ratio_s                Minimum match ratio for short motifs (Default: [1.,1.,1,0.55,0.5,0.65,0.55,0.625,0.65,0.5])" << std::endl;
    // std::cerr << " -ms, --match_score                       Match score (Default: 2)" << std::endl;
    // std::cerr << " -mp, --mismatch_penalty                  Mismatch penalty (Default: -1)" << std::endl;
    // std::cerr << " -gp, --gap_penalty                       Gap penalty (Default: -1)" << std::endl;
    std::cerr << std::endl;

    std::cerr << "[options]: (only for germline and somatic analysis)" << std::endl;
    std::cerr << "\t Read extraction parameters:" << std::endl;
    std::cerr << "\t -h, --help                               Show help message" << std::endl;
    std::cerr << "\t -s,  --output_file_statistics            Output file containing phasing information and consensus CN Call for each region" << std::endl;
    std::cerr << "\t -pad, --padding                          padding around the str region to extract reads from (Default: 0)" << std::endl;
    std::cerr << "\t -t, --threads                            Number of threads to use (Default: 1)" << std::endl;
    std::cerr << "\t -kpr, --keepPhasingResults               Keep phasing results (Default: false)" << std::endl;
    std::cerr << "\t -kcr, --keepCutReads                     keep cut reads (Default: false)"  << std::endl;
    std::cerr << "\t -minR, --minReadsInRegion                minimum number of reads that should span the region (Default: 2)" << std::endl;
    std::cerr << "\t -btg, --bamIsTagged                      reads in bam are phased (Default: false)" << std::endl;
    std::cerr << "\t -qs, --quality_score                     minimum quality score for a read to be considered (Default: 20 , Max:60)" << std::endl;
    std::cerr << std::endl;

    std::cerr << "\t Correction parameters:" << std::endl;
    std::cerr << "\t -cor, --correct                          perform genotype calling correction (CCS Default: false, CLR/ONT Default: true)" << std::endl;
    std::cerr << "\t -crs, --consensus_ratio_str              minimum percentage of reads in a cluster that have to call an intrval while correction (STR) (Default: 0.3)" << std::endl;
    std::cerr << "\t -crv, --consensus_ratio_vntr             minimum percentage of reads in a cluster that have to call an intrval while correction (VNTR) (Default: 0.3)" << std::endl;
    std::cerr << "\t -roz, --removeOutliersZscore             remove outliers for phasing (Default: false)" << std::endl;
    //std::cerr << " -rtr, --refineTrRegions                  refine tandem repeat regions (Default: true)" << std::endl;
    std::cerr << std::endl;

    std::cerr << "\t Clustering parameters:" << std::endl;
    std::cerr << "\t -seps, --start_eps_str                   the srart radian for clustering for str regions (Default: 0.7)" << std::endl;
    std::cerr << "\t -sepv, --start_eps_vntr                  the srart radian for clustering for vntr regions (Default: 0.6)" << std::endl;
    std::cerr << "\t -minPF, --minPts_frac                    the minimum fraction of reads that should be in one cluster (Default: 0.12)" << std::endl;
    std::cerr << "\t -nls, --noise_limit_str                  the noise limit for clustering for str regions (Default: 0.2)" << std::endl;
    std::cerr << "\t -nlv, --noise_limit_vntr                 the noise limit for clustering for vntr regions (Default: 0.35)" << std::endl;
    std::cerr << "\t -ci, --cluster_iter                      the number of iteration for clustering (Default: 40)" << std::endl;
    std::cerr << std::endl;

    std::cerr << "\t -v, --verbose                            Verbose mode Verbosity level (0:error, 1:critical, 2:info, 3:debug) (Default: 0)" << std::endl;
    std::cerr << "\t -h, --help                               print help message" << std::endl;
}

void printAssemblyDoc(){

    std::cerr << "--germline/--somatic/--assembly           Type of analysis to perform" << std::endl;
    std::cerr << " -b, --bam                                path to the bam file of the aligned reads to the reference genome" << std::endl;
    std::cerr << " -r, --ref                                Path to the input reference file [.fa/.fna]." << std::endl;
    std::cerr << " -m, --motif_file                         path to the file with reference coordinates and motif sequence [BED/TSV/CSV]" << std::endl;
    std::cerr << " -o, --output_file                        Output file containing region, motif, hap1 and hap2 copy numbers" << std::endl;
    std::cerr << " -s, --sex                                Sample sex (0|female 1|male)" << std::endl;
    std::cerr << " -sn, --sample                            name of sample" << std::endl;
    std::cerr << " -rt, --reads_type                        Type of reads (Default: CCS)" << std::endl;

    std::cerr << "[options]: Dynamic programming alignment Parameters: " << std::endl;
    std::cerr << " -mml, --min_match_ratio_l                Minimum match ratio for long motifs (Default: 0.5)" << std::endl;
    // std::cerr << " -ms, --match_score                       Match score (Default: 10)" << std::endl;
    // std::cerr << " -mp, --mismatch_penalty                  Mismatch penalty (Default: -2)" << std::endl;
    // std::cerr << " -gp, --gap_penalty                       Gap penalty (Default: -5)" << std::endl;
    std::cerr << std::endl;
}

void printGermlineSomaticDoc(){
    std::cerr << "Required input:" << std::endl;
    std::cerr << "--germline/--somatic/--assembly           Type of analysis to perform" << std::endl;
    std::cerr << " -b, --bam                                path to the bam file of the aligned reads to the reference genome" << std::endl;
    std::cerr << " -r, --ref                                Path to the input reference file [.fa/.fna]." << std::endl;
    std::cerr << " -m, --motif_file                         path to the file with reference coordinates and motif sequence [BED/TSV/CSV]" << std::endl;
    std::cerr << " -o, --output_file                        Output file containing region, motif, hap1 and hap2 copy numbers" << std::endl;
    std::cerr << " -s, --sex                                Sample sex (0|female 1|male)" << std::endl;
    std::cerr << " -sn, --sample                            name of sample" << std::endl;
    std::cerr << " -rt, --reads_type                        Type of reads (Default: CCS)" << std::endl;

    std::cerr << "[options]: Dynamic programming alignment Parameters: " << std::endl;
    std::cerr << " -mml, --min_match_ratio_l                Minimum match ratio for long motifs (Default: 0.5)" << std::endl;
    // std::cerr << " -ms, --match_score                       Match score (Default: 10)" << std::endl;
    // std::cerr << " -mp, --mismatch_penalty                  Mismatch penalty (Default: -2)" << std::endl;
    // std::cerr << " -gp, --gap_penalty                       Gap penalty (Default: -5)" << std::endl;
    std::cerr << std::endl;

    std::cerr << "[options]: (only for germline and somatic analysis)" << std::endl;
    std::cerr << "\t Read extraction parameters:" << std::endl;
    std::cerr << "\t -h, --help                               Show help message" << std::endl;
    std::cerr << "\t -s,  --output_file_statistics            Output file containing phasing information and consensus CN Call for each region" << std::endl;
    std::cerr << "\t -pad, --padding                          padding around the str region to extract reads from (Default: 0)" << std::endl;
    std::cerr << "\t -t, --threads                            Number of threads to use (Default: 1)" << std::endl;
    std::cerr << "\t -kpr, --keepPhasingResults               Keep phasing results (Default: false)" << std::endl;
    std::cerr << "\t -kcr, --keepCutReads                     keep cut reads (Default: false)"  << std::endl;
    std::cerr << "\t -minR, --minReadsInRegion                minimum number of reads that should span the region (Default: 2)" << std::endl;
    std::cerr << "\t -btg, --bamIsTagged                      reads in bam are phased (Default: false)" << std::endl;
    std::cerr << "\t -qs, --quality_score                     minimum quality score for a read to be considered (Default: 10 , Max:60)" << std::endl;
    std::cerr << std::endl;

    std::cerr << "\t Correction parameters:" << std::endl;
    std::cerr << "\t -cor, --correct                          perform genotype calling correction (CCS Default: false, CLR/ONT Default: true)" << std::endl;
    std::cerr << "\t -crs, --consensus_ratio_str              minimum percentage of reads in a cluster that have to call an intrval while correction (STR) (Default: 0.3)" << std::endl;
    std::cerr << "\t -crv, --consensus_ratio_vntr             minimum percentage of reads in a cluster that have to call an intrval while correction (VNTR) (Default: 0.3)" << std::endl;
    std::cerr << "\t -roz, --removeOutliersZscore             remove outliers for phasing (Default: false)" << std::endl;
    //std::cerr << " -rtr, --refineTrRegions                  refine tandem repeat regions (Default: true)" << std::endl;
    std::cerr << std::endl;

    std::cerr << "\t Clustering parameters:" << std::endl;
    std::cerr << "\t -seps, --start_eps_str                   the srart radian for clustering for str regions (Default: 0.2)" << std::endl;
    std::cerr << "\t -sepv, --start_eps_vntr                  the srart radian for clustering for vntr regions (Default: 0.2)" << std::endl;
    std::cerr << "\t -minPF, --minPts_frac                    the minimum fraction of reads that should be in one cluster (Default: 0.12)" << std::endl;
    std::cerr << "\t -nls, --noise_limit_str                  the noise limit for clustering for str regions (Default: 0.2)" << std::endl;
    std::cerr << "\t -nlv, --noise_limit_vntr                 the noise limit for clustering for vntr regions (Default: 0.35)" << std::endl;
    std::cerr << "\t -ci, --cluster_iter                      the number of iteration for clustering (Default: 40)" << std::endl;
    std::cerr << std::endl;

    std::cerr << "\t -v, --verbose                            Verbose mode Verbosity level (0:error, 1:critical, 2:info, 3:debug) (Default: 0)" << std::endl;
    std::cerr << "\t -h, --help                               print help message" << std::endl;

}
void printGermlineDoc(){
    std::cerr << "Usage: TandemTwister --germline [options] -b <input_bam_file> -m <motif_file> -r <input_reference_file> -o <output_file>  -s <sample_sex>" << std::endl;
    printGermlineSomaticDoc();
}


void TandemTwister::printUsage(std::unordered_map<std::string, std::string> &args) {
    
    if (args.find("--germline") != args.end() || args.find("--somatic") != args.end()) {
        printVersion();
        if (args.find("--germline") != args.end()) {
        std::cerr << "Germline analysis: Genotyping tandem repeats from long reads" << std::endl;
        }
        else{
            std::cerr << "Somatic analysis: Genotyping tandem repeats from long reads" << std::endl;
        }
        
        printGermlineDoc();
    }
    else if (args.find("--assembly") != args.end()) {
        printVersion();
        std::cerr << "Assembly analysis: Genotyping tandem repeats from aligned genome input" << std::endl;
        printAssemblyDoc();
    }
    else {
        printFullDoc();
    }
    //std::cerr << " -tanCon, --tandem_run_threshold          maximum number of bases for merging the tandem-repeats runs in the reference sequence (Default: 2 * motif_size)" << std::endl;
}

