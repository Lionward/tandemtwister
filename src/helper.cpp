#include "include/TandemTwister.hpp"


std::vector<uint16_t> TandemTwister::get_motif_ids(const std::vector<Interval>& path) {
    /**
     * @brief Get the motif ids.
     * 
     * This function gets the motif ids.
     * 
     * @param path The path.
     * @return The motif ids.
     * 
     */
    std::vector<uint16_t> motif_ids = {};
    for (const Interval& interval : path) {
        motif_ids.push_back(interval.motif_id);
    }
    return motif_ids;
}

std::string TandemTwister::get_motif_ids_str(const std::vector<uint16_t>& motif_ids) {
    /**
     * @brief Get the motif ids as a string.
     * 
     * This function gets the motif ids as a string.
     * 
     * @param motif_ids The motif ids.
     * @return The motif ids as a string.
     * 
     */
    if (motif_ids.empty()) {
        return "";
    }
    std::string motif_ids_str = std::accumulate(motif_ids.begin(), motif_ids.end(), std::string(),
        [](const std::string& a, uint16_t b) {
            return a.empty() ? std::to_string(b) : a + "_" + std::to_string(b);
        });
    return motif_ids_str;
}



std::string TandemTwister::generateRegionGenotypeResult(
    const std::vector<std::string>& motifs,
    const std::vector<std::vector<uint16_t>>& motif_ids_alleles,
    const std::string &region,
    std::vector<uint32_t> copy_numbers,
    uint16_t CN_reference_all_runs,
    uint16_t num_reads_covering_TR
    ) {


    /**
     * @brief Generate the region genotype result.
     * 
     * This function generates the region genotype result.
     * 
     * @param motifs The motifs.
     * @param motif_ids_alleles The motif ids for the alleles.
     * @param region The region.
     * @param copy_numbers The copy numbers of the alleles.
     * @param CN_reference_all_runs The copy number of the reference.
     * @param num_reads_covering_TR The number of reads covering the tandem repeat.
     * @return The region genotype result.
     * 
     */
    std::string motif_list = std::accumulate(motifs.begin(), motifs.end(), std::string(),
        [](const std::string& a, const std::string& b) {
            return a.empty() ? b : a + "," + b;
        });

    std::ostringstream region_genotype_result;
    std::string cns_string = std::accumulate(copy_numbers.begin(), copy_numbers.end(), std::string(),
        [](const std::string& a, uint32_t b) {
            return a.empty() ? std::to_string(b) : a + "\t" + std::to_string(b);
        });
    region_genotype_result << region << "\t"
                           << motif_list << "\t"
                           << cns_string << "\t"
                           << CN_reference_all_runs << "\t"
                           << num_reads_covering_TR << "\t";


    std::vector<std::string> motif_ids_list = {};
    for (const auto& motif_ids : motif_ids_alleles) {
        motif_ids_list.push_back(get_motif_ids_str(motif_ids));
    }
    for (const auto& motif_ids : motif_ids_list) {
        region_genotype_result << motif_ids << "\t";
    }
    region_genotype_result << "\n";


    return region_genotype_result.str();
}

std::vector<uint16_t> TandemTwister::computeLCP(const std::vector<std::string>& motifs) {
    /**
    * @brief  Compute the Longest Common Prefix (LCP) between motifs.

    * This function computes the Longest Common Prefix (LCP) between motifs. It compares each motif
    * with the previous motif and calculates the length of the common prefix between them.

    *
    * @param motifs  A vector of motifs.
    * @return A vector of integers representing the LCP between motifs.
    */

    std::vector<uint16_t> lcp(motifs.size(), 0);
    for (size_t i = 1; i < motifs.size(); ++i) {
        uint16_t j = 0;
        while (j < motifs[i - 1].size() && j < motifs[i].size() && motifs[i - 1][j] == motifs[i][j]) {
            ++j;
        }
        lcp[i] = j;
    }
    return lcp;
}


std::tuple<std::string, std::string, std::string> TandemTwister::parse_region(const std::vector<std::tuple<std::string,std::vector<std::string>>>& chunk ,size_t j) {
    /**
     * @brief Parse the region.
     * 
     * This function parses the region.
     * 
     * @param chunk The chunk.
     * @param j The index.
     * @return The region.
     * 
     */
    std::string region = std::get<0>(chunk[j]);
    std::string chr, start, end;
    std::stringstream ss(region);
    std::getline(ss, chr, ':');
    std::getline(ss, start, '-');
    std::getline(ss, end);
    return std::make_tuple(chr, start, end);
}

bool TandemTwister::isChrXorY(const std::string& chr) {
    return "chrX" == chr || "chrY" == chr;
}


void TandemTwister::open_bam_file(const std::string& bam_file, samFile*& fp, bam_hdr_t*& h, hts_idx_t*& idx) {
    /**
     * @brief Open the BAM file.
     * 
     * This function opens the BAM file.
     * 
     */
    fp = sam_open(bam_file.c_str(), "r");
    if (fp == NULL) {
        spdlog::error("Failed to open BAM file {}", bam_file);
        throw std::runtime_error("Error opening BAM file.");
    }

    h = sam_hdr_read(fp);
    if (h == NULL) {
        sam_close(fp);
        throw std::runtime_error("Error reading header for " + bam_file);
    }

    idx = sam_index_load(fp, bam_file.c_str());
    if (idx == NULL) {
        sam_close(fp);
        throw std::runtime_error("Failed to load index for " + bam_file);
    }
}


std::vector<std::string> split(const std::string& str, char delimiter) {
    /**
     * @brief Split a string by a delimiter.
     * 
     * This function splits a string by a delimiter.
     * 
     * @param str The string.
     * @param delimiter The delimiter.
     * @return A vector of strings.
     * 
     */
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}


std::vector<std::tuple<std::string, std::vector<std::string>>> TandemTwister::regions_not_grouped(std::string& regionFile) {
    /**
     * @brief Get the regions not grouped by chromosome, sorted for balanced parallel processing.
     * 
     * This function gets the regions sorted by region length to ensure balanced distribution
     * across multiple threads, preventing workload imbalance.
     * 
     * @param regionFile The region file.
     * @return A vector of tuples containing the region and the motifs, sorted for balanced processing.
     * 
     */
    // open the file
    std::ifstream region_file(regionFile);
    if (!region_file.is_open()) {
        spdlog::error("Error: Unable to open the region file");
        exit(1);
    }
    
    std::vector<std::tuple<std::string, std::vector<std::string>, int, std::string, int, int>> TR_regions_with_metadata;
    std::string line;
    
    // First pass: read all regions and calculate their lengths
    while (std::getline(region_file, line)) {
        std::istringstream iss(line);
        std::string chr, start, end, motifs;
        
        if (!(iss >> chr >> start >> end >> motifs)) {
            spdlog::error("Error: Incorrect format in line: {}", line);
            continue;
        }
        
        // Calculate region length for balancing
        int start_int = std::stoi(start);
        int end_int = std::stoi(end);
        int region_length = end_int - start_int;
        
        // Split the motifs by comma
        std::vector<std::string> motifList = split(motifs, ',');
        
        // Sort motifs by length from longest to shortest
        std::sort(motifList.begin(), motifList.end(), 
                 [](const std::string& a, const std::string& b) { 
                     return a.size() > b.size(); 
                 });
        
        std::string region = chr + ":" + start + "-" + end;
        
        TR_regions_with_metadata.push_back(std::make_tuple(region, motifList, region_length, chr, start_int, end_int));
    }
    region_file.close();
    
    // Sort regions by length for block disturbution
    std::sort(TR_regions_with_metadata.begin(), TR_regions_with_metadata.end(),
             [](const auto& a, const auto& b) {
                 return std::get<2>(a) > std::get<2>(b);
             });
    
    // Ensure we don't have more threads than regions
    if (this->num_threads > TR_regions_with_metadata.size()) {
        this->num_threads = TR_regions_with_metadata.size();
    }

    //tem block vector
    std::vector<std::vector<std::tuple<std::string, std::vector<std::string>, std::string, int, int>>> temp_blocks(this->num_threads);
    
    // Distribute the regions in a way to ensure that the blocks contains the same number of big regions.
    // this will optimize the run time
    for (size_t i = 0; i < TR_regions_with_metadata.size(); ++i) {
        auto& region_data = TR_regions_with_metadata[i];
        size_t block_index = i % this->num_threads;
        temp_blocks[block_index].push_back(
            std::make_tuple(
                std::get<0>(region_data), // region string
                std::get<1>(region_data), // motifs
                std::get<3>(region_data), // chr
                std::get<4>(region_data), // start
                std::get<5>(region_data)  // end
            )
        );
    }
    
    // Sort each block by chromosome, then start, then end
    for (auto& block : temp_blocks) {
        std::sort(block.begin(), block.end(),
                 [](const auto& a, const auto& b) {
                     if (std::get<2>(a) != std::get<2>(b)) {
                         return std::get<2>(a) < std::get<2>(b);
                     }
                     if (std::get<3>(a) != std::get<3>(b)) {
                         return std::get<3>(a) < std::get<3>(b);
                     }
                     return std::get<4>(a) < std::get<4>(b);
                 });
    }
    
    // Flatten the blocks 
    std::vector<std::tuple<std::string, std::vector<std::string>>> TR_regions;
    for (size_t block = 0; block < this->num_threads; ++block) {
        for (auto& region : temp_blocks[block]) {
            TR_regions.push_back(std::make_tuple(std::get<0>(region), std::get<1>(region)));
        }
    }
    
    return TR_regions;
}

std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>> TandemTwister::regions_grouped_by_chr(std::string  &regionFile,std::vector<std::string>& accepted_chromosome_names) {
    
    /*
    * @brief Get the regions grouped by chromosome.
    *
    * This function gets the regions grouped by chromosome.
    * 
    * @param regionFile The region file.
    * @param accepted_chromosome_names The accepted chromosome names.
    * @return A map containing the regions grouped by chromosome.
    *  
    */

    std::ifstream region_file(regionFile);
    if (!region_file.is_open()) {
        spdlog::error("Error: Unable to open the region file");
        exit(1);
    }

    std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>> TR_regions;
    std::string line;
    while (std::getline(region_file, line)) {
        std::istringstream iss(line);
        std::string chr;
        std::string start;
        std::string end;
        std::string motifs;


        if (!(iss >> chr >> start >> end >> motifs)) {
            spdlog::error("Error: Incorrect format in line: {}", line);
            continue;
        }
        std::vector<std::string> motifList = split(motifs, ',');
        std::sort(motifList.begin(), motifList.end(), [](const std::string& a, const std::string& b) { return a.size() > b.size(); });

        // if the chromosome is not in the accepted chromosome names, skip it
        if (std::find(accepted_chromosome_names.begin(), accepted_chromosome_names.end(), chr) == accepted_chromosome_names.end()) {
            continue;
        }
        // if the sample sex is female and the chromosome is Y, skip it
        if (this->sample_sex ==0 && chr == "chrY") {
            continue;
        }
        // if the sample
        if (TR_regions.find(chr) == TR_regions.end()) {
            TR_regions[chr] = std::vector<std::tuple<std::string, std::vector<std::string>>>();
        }

        TR_regions[chr].push_back(std::make_tuple(chr + ":" + start + "-" + end, motifList));
        

        
    }
    // get the unique keys in TR_regions
    for (const auto& key : TR_regions) {
        this->input_chromosome_names.push_back(key.first);
    }
    // sort the input_chromosome_names by chromosome name
    std::sort(this->input_chromosome_names.begin(), this->input_chromosome_names.end(), [](const std::string& a, const std::string& b) { return a < b; });
    
    // close the file
    region_file.close();
    return TR_regions;
}



