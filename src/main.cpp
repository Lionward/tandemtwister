#include "include/TandemTwister.hpp"
#include <iostream>

void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}


// Set log level based on your custom verbosity scale
void set_verbose_level(int verbosity) {
    switch (verbosity) {
        case 0: spdlog::set_level(spdlog::level::err);     break; // Error only
        case 1: spdlog::set_level(spdlog::level::critical); break; // Critical only
        case 2: spdlog::set_level(spdlog::level::info);     break; // Info and above
        case 3: spdlog::set_level(spdlog::level::debug);    break; // Debug and above
        default: spdlog::set_level(spdlog::level::err);     break; // Default to error
    }
}

// Initialize logging with your preferred settings
void init_logging(int verbosity = 0, const std::string& pattern = "[%H:%M:%S] [%^%l%$] %v") {
    // Create console logger with color
    auto console = spdlog::stdout_color_mt("console");
    spdlog::set_default_logger(console);
    
    // Set custom pattern
    spdlog::set_pattern(pattern);
    
    // Set verbosity level
    set_verbose_level(verbosity);
}


int main(int argc, char *argv[]) {
    auto start_time = std::chrono::high_resolution_clock::now();
    TandemTwister TandemTwister(argc, argv);
    // Print the arguments in a separate section
    spdlog::set_level(spdlog::level::info);
    spdlog::info("========================================");
    spdlog::info("Arguments used:");
    spdlog::info("----------------------------------------");
    spdlog::info("TandemTwister version: {}", TandemTwister.version);
    spdlog::info("Input bam file: {}", TandemTwister.input_bam);
    spdlog::info("Input reference file: {}", TandemTwister.input_reference);
    spdlog::info("Input region file: {}", TandemTwister.region_file);
    // spdlog::info("Output path: {}", TandemTwister.output_path);
    spdlog::info("Output file: {}", TandemTwister.output_file_vcf);
    spdlog::info("Padding: {}", TandemTwister.padding);
    spdlog::info("Number of threads: {}", TandemTwister.num_threads);
    spdlog::info("Analysis type: {}", TandemTwister.analysis_type);
    spdlog::info("Sample name: {}", TandemTwister.sampleName);
    spdlog::info("Sample sex: {}", (TandemTwister.sample_sex == 0 ? "female" : "male"));

    if (TandemTwister.analysis_type == "somatic" || TandemTwister.analysis_type == "germline") {
        spdlog::info("reads_type: {}", TandemTwister.reads_type);
        spdlog::info("Input reads are {}", TandemTwister.bamIsTagged == true ? "phased" : "not phased");
        spdlog::info("keep phasing results: {}", TandemTwister.keep_phasing_results == 0 ? "false" : "true");
        spdlog::info("keep cut reads: {}", TandemTwister.keep_cut_sequence == 0 ? "false" : "true");
        spdlog::info("remove outliers z-score: {}", TandemTwister.remove_outliers_zscore == 0 ? "false" : "true");
        spdlog::info("min_match_ratio_l: {}", TandemTwister.min_match_ratio_l);
        if (TandemTwister.min_match_ratio_s != 0){
            spdlog::info("min_match_ratio_s: {}", TandemTwister.min_match_ratio_s);
        }
        else {
            spdlog::info("min_match_ratio_s: [1.,1.,1,0.55,0.5,0.65,0.55,0.625,0.65,0.5]");
        }
        spdlog::info("start_eps_str: {}", TandemTwister.start_eps_str);
        spdlog::info("start_eps_vntr: {}", TandemTwister.start_eps_vntr);
        spdlog::info("minPts_frac: {}", TandemTwister.minPts_frac);
        spdlog::info("noise_limit_str: {}", TandemTwister.noise_limit_str);
        spdlog::info("noise_limit_vntr: {}", TandemTwister.noise_limit_vntr);
        spdlog::info("cluster_iter: {}", TandemTwister.cluster_iter);
        //spdlog::info("minPts_str: {}", TandemTwister.minPts_str);
        //spdlog::info("minPts_vntr: {}", TandemTwister.minPts_vntr);

        spdlog::info("correct genotype calling: {}", TandemTwister.correctCalling == 0 ? "false" : "true");
        if (TandemTwister.correctCalling) {
            spdlog::info("cons_ratio_str: {}", TandemTwister.cons_ratio_str);
            spdlog::info("cons_ratio_vntr: {}", TandemTwister.cons_ratio_vntr);
        }
        spdlog::info("Minimum number of reads in region: {}", TandemTwister.minReadsInRegion);
        spdlog::info("========================================");
        spdlog::info("Number of regions to process: {}", TandemTwister.TR_regions.size());
    }
    else{
        spdlog::info("min_match_ratio_l: {}", TandemTwister.min_match_ratio_l);

        if (TandemTwister.min_match_ratio_s != 0){
            spdlog::info("min_match_ratio_s: {}", TandemTwister.min_match_ratio_s);
        }
        else {
            spdlog::info("min_match_ratio_s: [1.,1.,1,0.55,0.5,0.65,0.55,0.625,0.65,0.5]");
        }
    }
    spdlog::info("========================================");

    init_logging(TandemTwister.verbose);
    //std::cout << "refine TanRegions: " << (TandemTwister.refineTrRegions == 0 ? "false" : "true") << std::endl;
    if (TandemTwister.refineTrRegions){
        std::cout << "tandem run threshold: " << TandemTwister.tanCon << std::endl;
    }
    // double rss_gb_before = getRSS();
    if (TandemTwister.analysis_type == "somatic" || TandemTwister.analysis_type == "germline") {
        TandemTwister.processRegionsForLongReadsInput();
    }
    else if (TandemTwister.analysis_type == "assembly") {
        TandemTwister.processRegionsForAssemblyInput();
    }
    else {
        std::cerr << "Error: bam_type should be either reads or assembly" << std::endl;
        exit(1);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;
    double vm, rss;
    process_mem_usage(vm, rss);
    // double rss_gb_after = getRSS();
    // double rss_gb_diff = rss_gb_after - rss_gb_before;
    spdlog::set_level(spdlog::level::info);
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
    // Print CPU time and RSS
    spdlog::info("========================================");
    spdlog::info("CPU time: {:.2f} sec, VM usage: {:.2f} GB, RSS usage: {:.2f} GB", elapsed_time.count(), vm / (1024 * 1024), rss / (1024 * 1024));
    spdlog::info("========================================");


    return 0;
}
