

#include "include/TandemTwister.hpp"


template<typename T>
T TandemTwister::getArg(const std::unordered_map<std::string, std::string>& args, const std::string& shortArg, const std::string& longArg, T defaultValue) {
    if (args.find(shortArg) != args.end()) {
        return static_cast<T>(std::stod(args.at(shortArg)));
    }
    else if (args.find(longArg) != args.end()) {
        return static_cast<T>(std::stod(args.at(longArg)));
    }
    else {
        return defaultValue;
    }
}

template<>

std::string TandemTwister::getArg<std::string>(const std::unordered_map<std::string, std::string>& args, const std::string& shortArg, const std::string& longArg, std::string defaultValue) {
    if (args.find(shortArg) != args.end()) {
        return args.at(shortArg);
    }
    else if (args.find(longArg) != args.end()) {
        return args.at(longArg);
    }
    else {
        return defaultValue;
    }
}

template<>
bool TandemTwister::getArg<bool>(const std::unordered_map<std::string, std::string>& args, const std::string& shortArg, const std::string& longArg, bool defaultValue) {
    if (args.find(shortArg) != args.end()) {
        return args.at(shortArg) == "true" || args.at(shortArg) == "1";
    }
    else if (args.find(longArg) != args.end()) {
        return args.at(longArg) == "true" || args.at(longArg) == "1";
    }
    else {
        return defaultValue;
    }
}




TandemTwister::TandemTwister(int argc, char *argv[]) {

    args = parseCommandLine(argc, argv);
    if (args.find("-h") != args.end() || args.find("--help") != args.end()) {
        printUsage(args);
        exit(0);
    }
    if (args.find("-b") != args.end() || args.find("--bam") != args.end()) {
        try {
            input_bam = args.at("-b");
        }
        catch (const std::out_of_range& oor) {
            input_bam = args.at("--bam");
        }
    }
    else {
        std::cerr << "Error: Missing input bam file" << std::endl;
        printUsage(args);
        exit(1);
    }

    if (args.find("-r") != args.end() || args.find("--ref") != args.end()) {
        try {
            input_reference = args.at("-r");
        }
        catch (const std::out_of_range& oor) {
            input_reference = args.at("--ref");
        }
    }
    else {
        std::cerr << "Error: Missing input reference file" << std::endl;
        printUsage(args);
        exit(1);
    }
    this->input_reference_fai = input_reference+".fai";

    if (args.find("-m") != args.end() || args.find("--motif_file") != args.end()) {
        try {
            region_file = args.at("-m");
        }
        catch (const std::out_of_range& oor) {
            region_file = args.at("--motif_file");
        }
    }
    else {
        std::cerr << "Error: Missing input motif file" << std::endl;
        printUsage(args);
        exit(1);
    }
    if (args.find("-s") !=   args.end() || args.find("--sex") != args.end()){
        try {
            sample_sex = std::stoi(args.at("-s"));

        }
        catch(const std::out_of_range& oor){
            sample_sex = std::stoi(args.at("--sex"));
        }
        if (sample_sex != 1 && sample_sex != 0){
            std::cerr << "Eroor: sample_sex should be either 0 or 1" << std::endl;
        }
    }
    else {
        std::cerr << "Error: Missing input sample sex" << std::endl;
        printUsage(args);
        exit(1);
    }
    std::string tmp_var = "";
    if (args.find("-o") != args.end() || args.find("--output_file") != args.end()) {
        try {
            tmp_var = args.at("-o");
        }
        catch (const std::out_of_range& oor) {
            tmp_var = args.at("--output_file");
        }
        //output_file = tmp_var +  ".tsv";
        this-> output_file_vcf = tmp_var +".vcf.gz";
        
    }
    else {
        std::cerr << "Error: Missing output file" << std::endl;
        printUsage(args);
        exit(1);
    }
    if (args.find("--somatic") != args.end() || args.find("--germline") != args.end()) {
        if (args.find("--somatic") != args.end() && args.find("--germline") != args.end()) {
            std::cerr << "Error: Please specify either somatic or germline analysis" << std::endl;
            printUsage(args);
            exit(1);
        }
        if (args.find("--somatic") != args.end()) {
            this->analysis_type = "somatic";
        }
        else {
            this->analysis_type = "germline";
        }
    }
    else if (args.find("--assembly") != args.end() && args.find("--germline") == args.end() && args.find("--somatic") == args.end()) {
        this->analysis_type = "assembly";
    }
    else {
        std::cerr << "Error: Missing type of analysis" << std::endl;
        printUsage(args);
        exit(1);
    }


    this->padding = getArg<uint16_t>(args,"-pad","--padding", 0);
    this->sampleName =  getArg<std::string>(args, "-sn", "--sample", "sample");
    this->cons_ratio_str = getArg<float>(args, "-crs", "--consensus_ratio_str", 0.4); 
    if (this->cons_ratio_str < 0 || this->cons_ratio_str > 1) {
        through_error(args,"Error: consensus_ratio_str should be between 0 and 1");
    }
    
    // the rest of the arguments are optional
    this->reads_type = getArg<std::string>(args, "-rt", "--reads_type", "CCS");
    // allow for only CCS , ONT and CLR reads as value for reads_type
    if (reads_type != "CCS" && reads_type != "ONT" && reads_type != "CLR") {
        std::cerr << "Error: reads_type should be either CCS, ONT or CLR" << std::endl;
        printUsage(args);
        exit(1);
    }
    this->verbose = getArg<uint16_t>(args, "-v", "--verbose", 0);


   




    this->match_score = getArg<uint16_t>(args, "-ms", "--match_score", 10);
    this->mismatch_score = getArg<int16_t>(args, "-mp", "--mismatch_penalty", -1);
    this->gap_score = getArg<int16_t>(args, "-gp", "--gap_penalty", -3);
    this->num_threads = getArg<uint32_t>(args, "-t", "--threads", 1);
    this->minReadsInRegion = getArg<uint16_t>(args, "-minR", "--minReadsInRegion", 4);
    this->bamIsTagged = getArg<bool>(args, "-btg", "--bamIsTagged", false);
    this->quality_score = getArg<uint16_t>(args, "-qs", "--quality_score", 10);
    if (this->quality_score > 60 ) {
        std::cerr << "Error: quality score should be less than 60, it is capped at 60" << std::endl;
        this->quality_score = 60;
    }
    if (this->num_threads < 1){
        std::cerr  << "Error: number of threads should be bigger than 0" << std::endl;
    }
 
    uint32_t pos = this->output_file.find_last_of("/"); 
    this->output_path = (pos != std::string::npos) ? this->output_file.substr(0, pos + 1) : "";

    // check if the output directory exists and create it if it does not 
    if (!this->output_path.empty() && !std::filesystem::exists(this->output_path)) {
        if (!std::filesystem::create_directories(this->output_path)) { 
            std::cerr << "Error: could not create directory " << this->output_path << std::endl;
            exit(1);
        }
    }



    //std::ofstream outfile(this->output_file);
    //this->outfile = std::move(outfile);
    
    if (!openFile(this->input_bam) || !openFile(this->input_reference) || !openFile(this->region_file)) {
        std::cerr << "Error: could not open one of the input files" << std::endl;
        exit(1);
    }
    if (this->analysis_type == "somatic" || this->analysis_type == "germline") {
        if (this->analysis_type == "somatic") {
            this->manual_reclustering = false;
            this->somatic = true;
        }
        else {
            this->manual_reclustering = true;
            this->somatic = false;
        }
        this->cons_ratio_vntr = getArg<float>(args, "-crv", "--consensus_ratio_vntr", 0.8);
        if (this->cons_ratio_str < 0 || this->cons_ratio_str > 1) {
            through_error(args,"Error: consensus_ratio_str should be between 0 and 1");
        }
        this->start_eps_str = getArg<float>(args, "-seps", "--start_eps_str", 0.2);
        if (this->start_eps_str < 0 ) {
            through_error(args,"Error: start_eps should be bigger than 0");
        }
        this->start_eps_vntr = getArg<float>(args, "-sepv", "--start_eps_vntr", 0.2);
        if (this->start_eps_str < 0) {
            through_error(args,"Error: start_eps should be bigger than 0");
        }

        this->minPts_vntr = getArg<uint16_t>(args, "-minPV", "--minPts_vntr", 2);
        
        if (noise_limit_str < 0 || noise_limit_str > 1) {
        through_error(args,"Error: noise_limit_str should be between 0 and 1");
        }

        this->noise_limit_vntr = getArg<float>(args, "-nlv", "--noise_limit_vntr", 0.3);
        if (noise_limit_str < 0 || noise_limit_str > 1) {
            through_error(args,"Error: noise_limit_str should be between 0 and 1");
        }
        this->cluster_iter = getArg<uint16_t>(args, "-ci", "--cluster_iter", 40);
        this->cluster_iter = this->cluster_iter == 0 ? 1 : this->cluster_iter;
        this->keep_phasing_results = getArg<bool>(args, "-kpr", "--keepPhasingResults", false);
        this->keep_cut_sequence = getArg<bool>(args, "-kcr", "--keepCutReads", false);
        this->remove_outliers_zscore = getArg<bool>(args, "-roz", "--removeOutliersZscore", false);
        this->refineTrRegions = getArg<bool>(args, "-rtr", "--refineTrRegions", false);
        this->minPts_frac = getArg<float>(args, "-minPF", "--minPts_frac", 0.12);
        this->min_match_ratio_l = getArg<double>(args, "-mml", "--min_match_ratio_l", 0.5);
        this->min_match_ratio_s = getArg<double>(args, "-mms", "--min_match_ratio_s", 0);
        this->correctCalling = getArg<bool>(args, "-cor", "--correct", false);
        this->minReadsInCluster = 0.1; // this is redundant with the minPts_frac (remove later)
        if (reads_type != "CCS"){
            this->motif_threashold_size =4;
            this-> tanCon = getArg<uint16_t>(args, "-tanCon", "--tandem_run_threshold", 1);
            this->noise_limit_str = getArg<float>(args, "-nls", "--noise_limit_str", 0.2);
        }
        else{
            this->motif_threashold_size = 3;
            this-> tanCon = getArg<uint16_t>(args, "-tanCon", "--tandem_run_threshold", 1);
            this->noise_limit_str = getArg<float>(args, "-nls", "--noise_limit_str", 0.4);
        }
        if (this->min_match_ratio_l < 0 || this->min_match_ratio_l > 1) {
            through_error(args,"Error: min_match_ratio should be between 0 and 1");
        }
        
        if (this->keep_phasing_results == true) {
            std::string output_file_phasing = output_path + "haps_phasing_" + output_file.substr(output_file.find_last_of("/") + 1);
            std::ofstream outfile_phasing(output_file_phasing);
            if (!outfile_phasing.is_open()) {
                std::cerr << "Error: could not open the pahsing file " << output_file_phasing << std::endl;
                exit(1);
            }
            this->outfile_phasing = std::move(outfile_phasing);
            this->outfile_phasing << "region\tread_name\tcluster\tCN\tconsensus_CN\n";
        }

        this->TR_regions = regions_not_grouped(region_file);
        this->num_threads = this->num_threads > TR_regions.size() ? TR_regions.size(): this->num_threads;

    }
    else {
        this->min_match_ratio_l = getArg<double>(args, "-mml", "--min_match_ratio_l", 0.5);
        //this->outfile << "region\tmotif\tCN_assembly\tCN_ref\tmotif_ids" << std::endl;
        this->TR_regions_by_chr = regions_grouped_by_chr(region_file,accepted_chromosome_names);
        this->num_threads = this->num_threads > TR_regions_by_chr.size() ? TR_regions_by_chr.size(): this->num_threads;
        this->keep_cut_sequence = getArg<bool>(args, "-kcr", "--keepCutReads", false);
    }

    
    if (this->keep_cut_sequence == true){
        std::string outputFileCutReads = output_path + "cut_reads_" + this->sampleName  + output_file.substr(output_file.find_last_of("/") + 1,output_file.find_last_of(".")) + ".fasta";
        std::ofstream outfile_cutReads(outputFileCutReads);
        this->outfile_cutReads = std::move(outfile_cutReads);
    
    }
    createVcfHeader();

}



// define the destructor
TandemTwister::~TandemTwister() {
    if (this->fp != NULL) {
        sam_close(this->fp);
    }
    if (this->h != NULL) {
        bam_hdr_destroy(this->h);
    }
    if (this->idx != NULL) {
        hts_idx_destroy(this->idx);
    }
    if (this->fai != NULL) {
        fai_destroy(this->fai);
    }

}


