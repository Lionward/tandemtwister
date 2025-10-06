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
        std::cerr << "Failed to open BAM file " << bam_file << std::endl;
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
     * @brief Get the regions not grouped by chromosome.
     * 
     * This function gets the regions not grouped by chromosome.
     * 
     * @param regionFile The region file.
     * @return A vector of tuples containing the region and the motifs.
     * 
     */
    // open the file
    std::ifstream region_file(regionFile);
    if (!region_file.is_open()) {
        std::cerr << "Error: Unable to open the region file" << std::endl;
        exit(1);
    }
    std::vector<std::tuple<std::string, std::vector<std::string>>> TR_regions;
    std::string line;
    while (std::getline(region_file, line)) {
        std::istringstream iss(line);
        std::string chr, start, end, motifs;
        // check if start and end are the same or the difference is one or negative
        if (!(iss >> chr >> start >> end >> motifs)) {
            std::cerr << "Error: Incorrect format in line: " << line << std::endl;
            continue;
        }
        // Split the motifs by comma
        std::vector<std::string> motifList = split(motifs, ',');
        // sort the motifs by length from the longest to the shortest
        std::sort(motifList.begin(), motifList.end(), [](const std::string& a, const std::string& b) { return a.size() > b.size(); });
        std::string region = chr + ":" + start + "-" + end;        
        TR_regions.push_back(std::make_tuple(region, motifList));
    }
    // close the file
    region_file.close();

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
        std::cerr << "Error: Unable to open the region file" << std::endl;
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
            std::cerr << "Error: Incorrect format in line: " << line << std::endl;
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



