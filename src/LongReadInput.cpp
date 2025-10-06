// candidate regions to compare to :
// chr1:109,564-109,614
// chr1:108444-108469,chr1:109564-109614 

#include "include/TandemTwister.hpp"
#include <thread>
// import logging module
#include <iostream>
// time 
#include <chrono>


bool TandemTwister::IsMotifRun(const std::vector<Interval> & intervals ,const uint16_t motif_size){
    /**
    * Check if the intervals are motif runs
    * @param intervals: The intervals to check
    * @param motif_size: The size of the motif
    * @return: True if the intervals are motif runs, false otherwise
    */
    if (intervals.empty()) {
        return false;
    }

    uint32_t maxMotifSpan = motif_size * intervals.size();
    uint32_t intervalSpan = intervals.back().end - intervals.front().start +1 ;
    uint32_t intervalSpan2 = 0;
    if (intervalSpan > maxMotifSpan){
        intervalSpan2 = maxMotifSpan;
    }
    else {
        intervalSpan2 = intervalSpan;
    }
    
    if (intervalSpan2 < 4 && motif_size == 1){ // for homopolymers we  only consider motif runs that are at least 5 bases long         
        return false;
    }
    else if (intervalSpan2 < 8 && motif_size > 1){
        return false;
    }

    else {
        uint32_t actualmotifspan = 0;
        const float threshold = 0.7f;
        if (intervals[0].start == intervals[0].end){
            actualmotifspan = intervals.size();
        }
        else{
            actualmotifspan = std::accumulate(intervals.begin(), intervals.end(), 0, 
                [](uint32_t sum, const Interval& interval) {
                    return sum + interval.end - interval.start + 1;
                });
        }
        return (static_cast<float>(actualmotifspan) / static_cast<float>(intervalSpan)) > threshold;
    }
}


std::vector<std::vector<Interval>> TandemTwister::findAllTandemRuns(const std::vector<Interval>& intervals, const uint16_t motif_size) {
    /**
     * Find all tandem runs in the intervals
     * @param intervals: The intervals to search for tandem runs
     * @param motif_size: The size of the motif
     * @return: A vector of vectors containing the tandem runs
     */
    if (intervals.empty()) {
        return std::vector<std::vector<Interval>>();
    }

    uint32_t threshold = this->tanCon * motif_size +1;

    uint16_t currentRunStart = 0;
    std::vector<std::vector<Interval>> allRuns;
    allRuns.reserve(intervals.size()); 
    // Split intervals into runs based on the threshold the we have 
    // example : if the sequence has this run --------ATATATATATATATAT--ATATATATA-----ATATATATAT------
    // the intervals will be split into 2 runs, the first run will be the first 16 intervals, since the distance between the 2nd run and the first is only 2, which is less than the threshold so they weill be merged into one run
    // the second run will be the last 10 intervals 
    for (uint32_t i = 1; i < intervals.size(); ++i) {
        const Interval& currentInterval = intervals[i];
        const Interval& previousInterval = intervals[i-1];
        if (std::abs(static_cast<int>(currentInterval.start) - static_cast<int>(previousInterval.end)) > threshold) {
            allRuns.push_back(std::vector<Interval>(intervals.begin() + currentRunStart, intervals.begin() + i));
            currentRunStart = i;
        }
    }
    allRuns.push_back(std::vector<Interval>(intervals.begin() + currentRunStart, intervals.end()));
    
    // Filter runs in one pass
    // 1- remove empty runs
    // 2- remove runs that are not motif runs
    // 3- remove runs that are not motif runs and have only one interval ex: like the last AT in this sequence.. --------ATATATATATATATAT--ATATATATA-----ATATATATAT------------AT--
    allRuns.erase(std::remove_if(allRuns.begin(), allRuns.end(), [this, motif_size,allRuns](const std::vector<Interval>& run) { 
        return run.empty() || run.size() == 1  || !IsMotifRun(run , motif_size) ; 
    }), allRuns.end());


    if (allRuns.empty()) {
        return std::vector<std::vector<Interval>>{intervals};
    }
    


    return allRuns;
}


std::vector<Interval> TandemTwister::genotype_reference(std::string &regionWithPadding, std::vector<std::string> & motifs ,std::string & seq, int len, faidx_t* fai, const std::vector<uint16_t> & lcp) {   
    /**
     * Genotype the reference sequence
     * @param regionWithPadding: The region to genotype
     * @param motifs: The motifs to search for in the sequence
     * @param seq: The sequence to search for the motifs
     * @param len: The length of the sequence
     * @param fai: The FASTA index
     * @param lcp: The longest common prefix array of the motifs to avoid redundant calculations
     * @return: A vector of intervals that meet the threshold for a match
     */
    seq = std::string(fai_fetch(fai, regionWithPadding.c_str(), &len));  
    if (seq.empty()) {
        std::cerr << "Failed to fetch sequence for region " << regionWithPadding << std::endl;
        return std::vector<Interval>();
    }
    //std::unordered_map<char, uint16_t> nucleotides_occurences;
    std::vector<Interval> path_reference = {};
   
    path_reference = findBestPath(seq, motifs, this->match_score, this->mismatch_score, this->gap_score, lcp);
    return path_reference;
}




std::string  TandemTwister::generate_phasing_results(std::vector<std::pair<std::string,std::vector<Interval>>> & first_cluster, std::vector<std::pair<std::string,std::vector<Interval>>> & second_cluster ,std::vector<std::pair<std::string,std::vector<Interval>>> &noise_cluster , unsigned int CN_first, unsigned int CN_second, const std::string& region) {
    /**
     * Generate the phasing results
     * @param first_cluster: The first cluster
     * @param second_cluster: The second cluster
     * @param noise_cluster: The noise cluster
     * @param CN_first: The copy number of the first cluster
     * @param CN_second: The copy number of the second cluster
     * @param region: The region
     * @return: The phasing results
     */

    std::string phasing_result = "";
    for (const auto& cluster : first_cluster) {
        phasing_result += region  + "\t" + cluster.first + "\t" + "0" + "\t" + std::to_string(cluster.second.size())  + "\t" + std::to_string(CN_first) +  "\n";
    }
    for (const auto& cluster : second_cluster) {
        phasing_result += region  + "\t" + cluster.first + "\t" + "1" + "\t" + std::to_string(cluster.second.size())  + "\t" + std::to_string(CN_second) + "\n";
    }
    for (const auto& cluster : noise_cluster) {
        phasing_result += region  + "\t" + cluster.first + "\t" + "-1" + "\t" + std::to_string(cluster.second.size())  + "\t" + "-1" + "\n";
    }
    return phasing_result;
}


std::vector<std::pair<uint32_t, uint32_t>> TandemTwister::cut_read_in_TR_region(const  uint32_t reference_start, uint32_t* cigar_data, size_t cigar_length,const std::vector<std::pair<uint32_t, uint32_t>>& tandem_runs) {
    /**
     * Cut the read in the tandem repeat region
     * @param reference_start: The reference start position
     * @param cigar_data: The cigar data
     * @param cigar_length: The length of the cigar data
     * @param tandem_runs: The tandem runs to process the read only in the regions of high purity 
     * @return: A vector of pairs containing the cut positions
    */

    uint32_t current_ref_position = reference_start;
    uint32_t current_read_position = 0;
    uint32_t start_read = 0;
    uint32_t end_read = 0;
    int32_t remaining_length = 0;
    std::vector<bool> reached_start(tandem_runs.size(), false);
    std::vector<std::pair<uint32_t, uint32_t>> tandem_runs_cut_positions(tandem_runs.size());
    tandem_runs_cut_positions.resize(tandem_runs.size());
    uint16_t tandem_run_idx = 0;
    uint16_t cigartype = 0;
    uint32_t cigartype_length =0;
    for (int i : std::views::iota(static_cast<int>(0), static_cast<int>(cigar_length))) {
        cigartype = bam_cigar_op(cigar_data[i]);
        cigartype_length = bam_cigar_oplen(cigar_data[i]);

        // if the cigar type is a deletion then we add the current ref position to the breakpoints with its size 
        switch (cigartype) {
            case 0:
            case 7:
            case 8: // match or mismatch or sequence match
                current_ref_position += cigartype_length;
                current_read_position += cigartype_length;
                break;
            case 1:
            case 4:
            case 3: // insertion or soft clippingv or hard clipping
                current_read_position += cigartype_length;
                break;
            case 2:
                current_ref_position += cigartype_length;
                break;
            default:
                break;
        }
        remaining_length = 0;

        if (reached_start[tandem_run_idx] != true && current_ref_position >=  tandem_runs[tandem_run_idx].first-1) {
            if (cigartype != 2  && cigartype != 1) { 
                remaining_length = current_ref_position -  tandem_runs[tandem_run_idx].first;
            }
            else if (cigartype == 1){
                remaining_length = current_ref_position - tandem_runs[tandem_run_idx].first + cigartype_length ;
            }
            remaining_length = remaining_length < 0 ? 0 : remaining_length;
            start_read = current_read_position - remaining_length;
            reached_start[tandem_run_idx] = true;
            tandem_runs_cut_positions[tandem_run_idx].first = start_read;
        }
        remaining_length = 0;
        if (current_ref_position >=  tandem_runs[tandem_run_idx].second) {
            if (cigartype != 2) { 
                remaining_length = current_ref_position - tandem_runs[tandem_run_idx].second;
            }
            end_read = current_read_position - remaining_length;
            // reached_end = true;
            tandem_runs_cut_positions[tandem_run_idx].second = end_read;
            ++tandem_run_idx;
            // consider the case if the current ref postion is also over the next tandem run
            while(true){
                if (tandem_run_idx == tandem_runs.size()) {
                    break;
                }
                
                remaining_length = 0;
                if (reached_start[tandem_run_idx] != true &&  current_ref_position >=  tandem_runs[tandem_run_idx].first -1) {

                    if (cigartype != 2  && cigartype != 1 ) { //  && cigartype != 6 && cigartype != 5
                        remaining_length = current_ref_position -  tandem_runs[tandem_run_idx].first ;
                    }
                    else if (cigartype == 1){
                        remaining_length = current_ref_position   - tandem_runs[tandem_run_idx].first  + cigartype_length ;
                    }
                    remaining_length = remaining_length < 0 ? 0 : remaining_length;
   
                    start_read = current_read_position - remaining_length;

                    tandem_runs_cut_positions[tandem_run_idx].first = start_read;
                    reached_start[tandem_run_idx] = true;
                }
                remaining_length = 0;
                if (current_ref_position >=  tandem_runs[tandem_run_idx].second  ) {
                    if (cigartype != 2 && cigartype != 6 && cigartype != 5) {
                        remaining_length = current_ref_position - tandem_runs[tandem_run_idx].second;
                    }
                    
                    end_read = current_read_position - remaining_length;
                    tandem_runs_cut_positions[tandem_run_idx].second = end_read;
                    ++tandem_run_idx;
                }
                else {
                    break;
                }
            }
            if (tandem_run_idx == tandem_runs.size()) {
                break;
            }
        }

    }

    return tandem_runs_cut_positions;

}




std::tuple<uint32_t, uint32_t, uint16_t> TandemTwister::findMostOccurringCopyNumber(const std::vector<std::pair<std::string, std::vector<Interval>>>& cluster) {
    /**
     * Find the most occurring copy number in the cluster
     * @param cluster: The cluster to search for the most occurring copy number
     * @return: A tuple containing the most occurring copy number, the frequency of the copy number, and the index of the most occurring copy number
     */

    if (cluster.empty()) {
        return std::make_tuple(0, 0, 0);
    }
    std::unordered_map<uint32_t, uint32_t> frequencyMap;

    for (const auto& tuple : cluster) {
        uint32_t copyNumber = tuple.second.size();

        ++frequencyMap[copyNumber];
    }


    // if the size of the cluster is 2 then check if the CN of both is the same then return the index of the first read, if not return the index of the maximum CN
    if (cluster.size() == 2) {
        if (cluster[0].second.size() == cluster[1].second.size()) {
            return std::make_tuple(cluster[0].second.size(), 2, 0);
        }
        else {
            return std::make_tuple(cluster[0].second.size() > cluster[1].second.size() ? cluster[0].second.size() : cluster[1].second.size(), 1, cluster[0].second.size() > cluster[1].second.size() ? 0 : 1);
        }
    }

    // if there are more than 2 copy numbers with the same frequency then return the index of the highest one


    auto maxElement = std::max_element(frequencyMap.begin(), frequencyMap.end(),
        [](const std::pair<uint32_t, uint32_t>& a, const std::pair<uint32_t, uint32_t>& b) {
            return a.second < b.second;
        });
    // get the second most occuring copy number
    auto secondMaxElement = std::max_element(frequencyMap.begin(), frequencyMap.end(),
        [maxElement](const std::pair<uint32_t, uint32_t>& a, const std::pair<uint32_t, uint32_t>& b) {
            return a.second < b.second && a.first != maxElement->first;
        });
    
    // if the second most occuring copy number has the same frequency as the most occuring copy number then return the index of the highest copy number
    if (secondMaxElement != frequencyMap.end() && secondMaxElement->second == maxElement->second) {
        auto maxCopyNumber = std::max(maxElement->first, secondMaxElement->first);
        for (size_t i = 0; i < cluster.size(); ++i) {
            if (cluster[i].second.size() == maxCopyNumber) {
                return std::make_tuple(maxCopyNumber, maxElement->second, i);
            }
        }
    }


  
    uint32_t mostOccurringCopyNumber = maxElement->first;
    uint32_t frequency = maxElement->second;
    
    uint16_t index = 0;
    for (size_t i = 0; i < cluster.size(); ++i) {
        if (cluster[i].second.size() == mostOccurringCopyNumber) {
            index = i;
            break;
        }
    }

    return std::make_tuple(mostOccurringCopyNumber, frequency, index);
}
std::tuple<uint32_t, uint32_t, uint16_t> TandemTwister::findMedianCopyNumber(const std::vector<std::pair<std::string, std::vector<Interval>>>& cluster) {
    /**
     * Find the median copy number in the cluster
     * @param cluster: The cluster to search for the median copy number
     * @return: A tuple containing the median copy number, the frequency of the median copy number, and the index of the median copy number
     */

    if (cluster.empty()) {
        return std::make_tuple(0, 0, 0);
    }

    std::vector<uint32_t> copyNumbers;
    for (const auto& tuple : cluster) {
        copyNumbers.push_back(tuple.second.size());
    }

    std::sort(copyNumbers.begin(), copyNumbers.end());

    uint32_t medianCopyNumber = copyNumbers[copyNumbers.size() / 2];
    uint32_t frequency = std::count(copyNumbers.begin(), copyNumbers.end(), medianCopyNumber);

    uint16_t index = 0;
    for (size_t i = 0; i < cluster.size(); ++i) {
        if (cluster[i].second.size() == medianCopyNumber) {
            index = i;
            break;
        }
    }

    return std::make_tuple(medianCopyNumber, frequency, index);
}


void TandemTwister::print_dashes_representation(const std::string& read_name, const std::string& cut_read, const std::vector<Interval>& path) {
    /**
     * this function is a helper function to print the dashes representation of the read
     * @param read_name: The name of the read
     * @param cut_read: The cut read
     * @param path: The path to print the dashes representation
     */
    if (path.empty()) {
        return;
    }
    std::string read_seq(cut_read.length(), '-');
    for (const auto& interval : path) {
        read_seq.replace(interval.start, interval.end - interval.start, cut_read.substr(interval.start-1, interval.end - interval.start+1));
    }
    std::cout << read_name << "\t" << read_seq << std::endl;
}

void TandemTwister::cluster_by_features(const std::string& chr, 
                                        std::vector<std::tuple<std::string, std::vector<Interval>, uint8_t>>& reads_intervals, 
                                        std::vector<std::vector<std::pair<std::string, std::vector<Interval>>>> & clusters, 
                                        std::vector<std::pair<std::string, std::vector<Interval>>>& noise_cluster, 
                                        arma::mat &region_features, 
                                        uint16_t motif_size,
                                        bool &haplody) {
    /**
     * Cluster reads by features (e.g. cn , motif occurences , length, purity, score)
     * @param chr: The chromosome
     * @param reads_intervals: The reads intervals
     * @param clusters: The clusters
     * @param noise_cluster: The noise cluster
     * @param region_features: The region features
     * @param motif_size: The size of the motif
     * @param haplody: The haplody flag
     */

    if (!((this->sample_sex == 1) && (isChrXorY(chr)))) {
        cluster_reads(region_features, reads_intervals, clusters, noise_cluster, motif_size);
    } else {
        haplody = true;
        cluster_reads(region_features, reads_intervals, clusters, noise_cluster, motif_size);
        clusters.resize(2);

        // replace the second cluster with the first cluster if the copy number of the second cluster is more than the first cluster (this )
        uint cn_1_avg = 0;
        uint cn_2_avg = 0;
        for (const auto& read : clusters[0]) {
            cn_1_avg += read.second.size();
        }
        cn_1_avg /= clusters[0].size();

        for (const auto& read : clusters[1]) {
            cn_2_avg += read.second.size();
        }
        cn_2_avg /= clusters[1].size();

        if (cn_1_avg > cn_2_avg) {
            clusters[1] = clusters[0];
        }
        else {
            clusters[0] = clusters[1];
        }

        clusters[1] = clusters[0];
    }
}


void TandemTwister::get_hp_tag(uint8_t *&hp_tag, bam1_t *reads, uint32_t &hp_value)
{
    if (this->bamIsTagged)
    {
        hp_tag = bam_aux_get(reads, "HP");
        if (hp_tag)
        {
            hp_value = bam_aux2i(hp_tag);
        }
    }
    else
    {
        hp_value = 0;
    }
}

regionResult TandemTwister::processLongReads(samFile* &fp, hts_itr_t* &iter,const size_t start_pos, const size_t end_pos, const std::string& region, 
                                    std::vector<std::string>& motifs, uint16_t motif_size, ::string &chr,
                                    std::vector<std::pair<uint32_t, uint32_t>> &tandem_runs,
                                    std::string &region_phasing_result, std::string &cut_reads_fasta, const std::vector<uint16_t>& lcp) {
    /**
     * Process long reads
     * @param fp: The file pointer
     * @param iter: The iterator
     * @param start_pos: The start position
     * @param end_pos: The end position
     * @param region: The region
     * @param motifs: The motifs
     * @param chr: The chromosome
     * @param tandem_runs: The tandem runs
     * @param region_phasing_result: The region phasing result
     * @param num_reads: The number of reads
     * @param cut_reads_fasta: The cut reads in FASTA format
     * @param lcp: The longest common prefix array
     */

    // define the variables that will be outputted
    arma::mat region_features; 
    std::vector<std::vector<std::pair<std::string, std::vector<Interval>>>> clusters;
    std::vector<std::pair<std::string, std::vector<Interval>>>  noise_cluster;
    std::vector<std::string> alt_seqs;
    std::vector<uint16_t> indices_consensus_read = {};
    std::vector<uint32_t> copy_numbers = {};
    std::vector<uint32_t> cn_occurences_in_cluster = {};
    uint16_t num_reads_covering_TR = 0;


    regionResult region_results{clusters, noise_cluster, alt_seqs, indices_consensus_read, copy_numbers, cn_occurences_in_cluster, num_reads_covering_TR};
    bam1_t* reads = bam_init1();
    uint8_t* seq_data;
    char* qseq;
    std::string read_name="";
    uint32_t* cigar_data;
    size_t read_reference_start;
    std::vector<Interval> path;
    // vector of tuples containing the read name, the path of the read, and the hp tag value
    std::vector<std::tuple<std::string,std::vector<Interval>,uint8_t >> reads_intervals = {};

    // map containing the read name and the read sequence
    std::unordered_map<std::string, std::string> reads_sequences = {};
    
    // map containing the read sequence and the path of the read (avoid processing identical reads)
    std::unordered_map<std::string, std::vector<Interval>> processed_reads;

    // map containing the motif and the number of occurences
    std::unordered_map<std::string, uint16_t> motif_occurences;

    // variables to store the read features
    float cut_read_length = 0.0f;
    float CN_total = 0.0f;
    float path_purity = 0.0f;
    float path_total_score = 0.0f;
    uint32_t start = 0;
    uint32_t end = 0;
    uint32_t length = 0;
    uint16_t region_features_idx = 0;
    uint8_t* hp_tag = 0;
    uint32_t hp_value = 0;

    // reserve the memory for the processed reads
    processed_reads.reserve(200);

    bool first_cluter = false;
    bool second_cluster = false;

    // will be added later to the command

    spdlog::info("Processing region {}", region);
    
    // iterate over the reads
    while (sam_itr_next(fp, iter, reads) >= 0) {

        if (reads->core.flag & BAM_FSECONDARY || reads->core.flag & BAM_FSUPPLEMENTARY || reads->core.flag & BAM_FDUP || reads->core.flag & BAM_FQCFAIL || reads->core.pos >  static_cast<hts_pos_t>(start_pos) || bam_endpos(reads) <  static_cast<hts_pos_t>(end_pos) || (static_cast<int>(reads->core.qual) < this->quality_score)) {
            continue;
        }

        if (num_reads_covering_TR > 200){
            break;
        }

        ++num_reads_covering_TR;

        // get the read sequence
        seq_data = bam_get_seq(reads);

        // get the read name
        read_name = std::string(bam_get_qname(reads));

        // get the hp tag value
        get_hp_tag(hp_tag, reads, hp_value);

    
        if (hp_value == 1){
            first_cluter = true;
        }
        else if (hp_value == 2){
            second_cluster = true;
        }
        // get the cigar data
        cigar_data = bam_get_cigar(reads);  

        // get the read reference start
        read_reference_start = reads->core.pos;

        // cut the read in the tandem repeat region
        std::vector<std::pair<uint32_t, uint32_t>>  tandemRepeatsCutPositions = cut_read_in_TR_region(read_reference_start , cigar_data,reads->core.n_cigar, tandem_runs);
        
        // if the read is not in the tandem repeat region then continue (this could happen if the read is not covering the whole tandem repeat region)
        if (tandemRepeatsCutPositions.empty()) {
            continue;
        }

        // get the sequence of the read
        qseq  = new char[tandemRepeatsCutPositions.back().second +5];
        
        for (uint32_t i = 0; i < tandemRepeatsCutPositions.back().second+4; ++i) {
            qseq[i] = seq_nt16_str[bam_seqi(seq_data, i)]; // convert the sequence to string
        }
        qseq[tandemRepeatsCutPositions.back().second+4] = '\0'; // null terminate the string

        region_features_idx = 0;
        start = tandemRepeatsCutPositions[0].first;
        end = tandemRepeatsCutPositions.back().second;
        length = end - start +1;
        start = start > 0 ? start : 1;
        std::string cut_read(qseq + start-1,length);

        // add the cut read to the cut reads in FASTA format (this is used for debugging) -kcr option in commanline
        if (this->keep_cut_sequence){
            cut_reads_fasta += region + ";" + read_name + "\n" + cut_read + "\n";
        }

        // store the read sequence
        reads_sequences[read_name] = cut_read;
        

        cut_read.clear();
        std::pair<std::string,std::vector<Interval>> path_dict;
        std::vector<Interval> path_total = {};
        
        for (const auto &pos : tandemRepeatsCutPositions) {
            start = pos.first;
            end = pos.second;
            length = end - start + 1;  
            std::string cut_read_run(qseq + start -1, length);
            path.clear();

            // check if the read is already processed
            if (processed_reads.find(cut_read_run) != processed_reads.end()) {
                path = processed_reads[cut_read_run];
            }
            else {
                path = findBestPath(cut_read_run, motifs, this->match_score, this->mismatch_score, this->gap_score,lcp);
                processed_reads[cut_read_run] = path;
            }
            
            // caluclate the occurences of the motifs in the read
            for (const auto& interval : path) {
                ++motif_occurences[motifs[interval.motif_id]];
            }
            uint32_t  path_size = path.size();

            // adjust the path intervals to the total path
            if (region_features_idx> 0 && !path_total.empty()){
                auto temp = start - tandemRepeatsCutPositions[region_features_idx-1].second  + path_total.back().end ;
                for (auto &interval : path) {
                    interval.start += temp;
                    interval.end += temp;
                }
            }

            // if the path is empty then continue
            if (path_size == 0)
            {
                ++region_features_idx;
                continue; 
            }

            // calculate the read features
            cut_read_length += cut_read_run.length();

            // path puriity is how much of a read is covered by the path
            path_purity = std::accumulate(path.begin(), path.end(), 0, 
                [](uint32_t sum, const Interval& interval) {
                    return sum + interval.end - interval.start + 1;
                });
            
            // path total score is the sum of the scores of the path intervals
            path_total_score = std::accumulate(path.begin(), path.end(), 0, 
                [](uint32_t sum, const Interval& interval) {
                    return sum + interval.score;
                });

            // CN total is the sum of cn for each subregion
            CN_total += path_size;
            
            path_total.insert(path_total.end(), path.begin(), path.end());


            std::string().swap(cut_read_run);
            ++region_features_idx;
            length = 0;
        }


        // store the read intervals
        reads_intervals.emplace_back(std::make_tuple(read_name, path_total,hp_value));

        // store the read features
        arma::rowvec read_features = {cut_read_length, CN_total, path_purity, path_total_score};

        // add the motif occurences to the read features
        for (const auto& motif : motifs) {
            read_features.insert_cols(read_features.n_cols, arma::rowvec({static_cast<double>(motif_occurences[motif])}));
        }

        // insert the read features to the reads features matrix
        region_features.insert_rows(region_features.n_rows, read_features);


        cut_read_length = 0.0f;
        path_purity = 0.0f;
        path_total_score = 0.0f;
        CN_total = 0.0f;
        delete[] qseq;
        motif_occurences.clear();
        tandemRepeatsCutPositions.clear();
        path_total.clear();

    }
    spdlog::info("Number of reads in region {}: {}", region, reads_intervals.size());
    bam_destroy1(reads);
    bool haplody = false;


    // this parameter is controlled by -minR  parameters which is the minimum number of reads in the region for the region to be processed
    if (region_features.n_rows < this->minReadsInRegion) {
        spdlog::warn("Not enough reads in region {}: {}", region, region_features.n_rows);
        clusters.resize(2);
        alt_seqs.resize(2, "");
        copy_numbers.resize(2, 0);
        cn_occurences_in_cluster.resize(2, 0);
        indices_consensus_read.resize(2, 0);        
        region_results = regionResult{clusters, noise_cluster, alt_seqs, indices_consensus_read, copy_numbers, cn_occurences_in_cluster, num_reads_covering_TR};
        return region_results;
    }


    // if the bam is tagged then cluster the reads based on the hp tag (for now other reads that are not tagged will be ignored)
    if (this->bamIsTagged){
        clusters.resize(2);
        if (first_cluter && second_cluster){
            
            for (const auto& read : reads_intervals) {
                if (std::get<2>(read) == 1){
                    clusters[0].push_back(std::make_pair(std::get<0>(read),std::get<1>(read)));
                }
                else if (std::get<2>(read) == 2){
                    clusters[1].push_back(std::make_pair(std::get<0>(read),std::get<1>(read)));
                }
                else {
                    noise_cluster.push_back(std::make_pair(std::get<0>(read),std::get<1>(read)));
                }
            }
            if (clusters[0].empty() && clusters[1].empty()) {
                spdlog::warn("No reads in region {}: {}", region, region_features.n_rows);
                clusters.clear();
                first_cluter = false;
                second_cluster = false;
            }
        }
        else if (first_cluter || second_cluster) {
            clusters[0].reserve(reads_intervals.size());
            for (const auto& read : reads_intervals) {
                if (std::get<2>(read) == 1 || std::get<2>(read) == 2){
                    clusters[0].push_back(std::make_pair(std::get<0>(read),std::get<1>(read)));
                }
                else {
                    noise_cluster.push_back(std::make_pair(std::get<0>(read),std::get<1>(read)));
                }
            }
        }
        else {

            cluster_by_features(chr, reads_intervals, clusters, noise_cluster, region_features, motif_size, haplody);
        }
    }
    else{
        cluster_by_features(chr, reads_intervals, clusters, noise_cluster, region_features, motif_size, haplody);
    }

    // do correction for the copy numbers for motif sizes less than 10
    if  (this->correctCalling && motif_size < 10) {
        spdlog::debug("Correcting copy numbers");
        float consensus_ratio = (motif_size <= this->STR_max_length) ? this->cons_ratio_str : this->cons_ratio_vntr;
        for (auto& cluster : clusters) {
            spdlog::debug("Correcting copy numbers for cluster with size {}", cluster.size());
            if (cluster.size() < 2) {
                spdlog::debug("Cluster size is less than 2");
                if (cluster.size() == 0) {
                    copy_numbers.push_back(0);
                    alt_seqs.push_back("");
                }
                else {
                    copy_numbers.push_back(cluster[0].second.size());
                    alt_seqs.push_back(reads_sequences[cluster[0].first]);
                }
                cn_occurences_in_cluster.push_back(1);
                indices_consensus_read.push_back(0);
                continue;
            }
            uint16_t consensus_range = cluster.size() <= 2 ? 2 : std::ceil(cluster.size() * consensus_ratio);
            uint16_t cn = correct_intervals(cluster, consensus_range, motif_size, reads_sequences);

            // if the corrected copy number is 0 it means the main read has no intervals that pass the consensus range
            if (cn == 0) {
                spdlog::debug("No intervals that pass the consensus range after correction");
                copy_numbers.push_back(0);
                alt_seqs.push_back("");
                cn_occurences_in_cluster.push_back(0);
                indices_consensus_read.push_back(0);
                continue;
            }

            // after removing the intervals that didn't pass the consensus range we need to cut the sequence from the start of the first interval to the end of the last interval (if intermediate intervals are removed they will appear as interrupted sequences; consider removing this part of the sequence if the intermediate intervals are removed)
            std::string cut_seq = reads_sequences[cluster[0].first].substr(cluster[0].second[0].start - 1, cluster[0].second.back().end - cluster[0].second[0].start + 1);
            
            // reset the intervals starts and ends by subtracting the start of the first interval (this is to make the intervals start from 0) after removing the intervals that didn't pass the consensus range
            auto start_interval = cluster[0].second[0].start;
            for (size_t j = 0; j < cluster[0].second.size(); ++j) {

                cluster[0].second[j].start -= start_interval - 1;
                cluster[0].second[j].end -= start_interval - 1;

            }
            
            // add the corrected cluster info to the cluster
            copy_numbers.push_back(cn);
            alt_seqs.push_back(cut_seq);
            indices_consensus_read.push_back(0);
            cn_occurences_in_cluster.push_back(cluster.size());

        }
    } 
    else {
        for (const auto& cluster : clusters) {
            if (cluster.size() > 0) {
                uint16_t cn = 0;
                uint16_t occurences = 0;
                uint16_t index_cn = 0;
                std::tie(cn, occurences, index_cn) = findMostOccurringCopyNumber(cluster);
                copy_numbers.push_back(cn);
                alt_seqs.push_back(reads_sequences[cluster[index_cn].first]);
                indices_consensus_read.push_back(index_cn);
                cn_occurences_in_cluster.push_back(occurences);
            }
            else {
                copy_numbers.push_back(0);
                alt_seqs.push_back("");
                cn_occurences_in_cluster.push_back(0);
                indices_consensus_read.push_back(0);

            }
        }
    }


    if (copy_numbers[0] == 0 && copy_numbers[1] == 0) {
        region_results = regionResult{clusters, noise_cluster, alt_seqs, indices_consensus_read, copy_numbers, cn_occurences_in_cluster, num_reads_covering_TR};
        return region_results;
    }

    if (clusters.size() == 2) {
        const float homopolymerThreshold = 0.2f;
        const uint16_t cnDifferenceThreshold = 2;
        float homozygousThreshold = 0.7f;
        const uint16_t minReadsForCorrection = 7;

      
        bool isHomopolymerCorrection = (motif_size == 1  
                                        && std::abs(static_cast<int>(copy_numbers[0]) - static_cast<int>(copy_numbers[1])) < cnDifferenceThreshold 
                                        && clusters[1].size() < std::round(homopolymerThreshold * (region_features.n_rows - noise_cluster.size())));
 
        bool isHomozygousCorrection = ((copy_numbers[1] == 0 
                                      && clusters[0].size() > std::round(homozygousThreshold * (region_features.n_rows - noise_cluster.size()))) 
                                      && region_features.n_rows >= minReadsForCorrection);
        

        if (isHomopolymerCorrection || isHomozygousCorrection){ 
            spdlog::debug("isHomopolymerCorrection: {} isHomozygousCorrection: {}", isHomopolymerCorrection, isHomozygousCorrection);
            copy_numbers[1] = copy_numbers[0];
            alt_seqs[1] = alt_seqs[0];
            clusters[1] = clusters[0];
            indices_consensus_read[1] = indices_consensus_read[0];
            cn_occurences_in_cluster[1] = cn_occurences_in_cluster[0];
        }
    }


    spdlog::debug("Copy numbers: {} {}", copy_numbers[0], copy_numbers[1]);
    

    if (this->keep_phasing_results){
        region_phasing_result += generate_phasing_results(clusters[0], clusters[1], noise_cluster, copy_numbers[0], copy_numbers[1],region); 
    }
    region_results = regionResult{clusters, noise_cluster, alt_seqs, indices_consensus_read, copy_numbers, cn_occurences_in_cluster, num_reads_covering_TR};
    reads_intervals.clear();
    reads_sequences.clear();
    return region_results;
}



std::tuple<std::vector<std::string> ,std::vector<vcfRecordInfoReads>> TandemTwister::process_chunk_reads(size_t start_dict, size_t end_dict, samFile* fp, bam_hdr_t* h, hts_idx_t* idx, faidx_t* fai,std::string & cut_reads_fasta) {
    /**
     * Process chunk reads
     * @param start_dict: The start dictionary
     * @param end_dict: The end dictionary
     * @param fp: The file pointer
     * @param h: The header
     * @param idx: The index
     * @param fai: The FASTA index
     * @param cut_reads_fasta: The cut reads in FASTA format
     * @return: A tuple containing the phasing results, the genotype results, and the genotype records
     */

    std::vector<Interval> path_reference;
    std::string chr, start, end;
    std::string region;
    uint32_t  start_pos = 0; // this is the start position of the region to extract reads from
    uint32_t end_pos= 0; // this is the end position of the region to extract reads from
    uint32_t ref_seq_len = 0; // length of the reference sequence
    std::string  seq = "";
    std::string region_phasing_result; // this is the phasing result for the region
    // std::string region_genotype_result; // this is the genotype result for the region in Bed format
    std::vector<std::string> chunk_genotype_results; // this is a vector of genotype results for the regions for each process
    std::vector<vcfRecordInfoReads> chunk_genotype_records; // this is a vector of genotype results for the regions for each process
    std::vector<std::string> chunk_phasing_results;
    std::string TR_type ="";
    hts_itr_t* iter = nullptr;
    uint16_t motif_size = 0;

    for (size_t j = start_dict; j < end_dict; ++j){
        std::vector<std::string> motifs = std::get<1>(this->TR_regions[j]);
        std::sort(motifs.begin(), motifs.end());
        const std::vector<uint16_t> lcp = computeLCP(motifs);
        // get the max length of motifs as motif size
        motif_size = std::max_element(motifs.begin(), motifs.end(), 
            [](const std::string& a, const std::string& b) {
            return a.size() < b.size();
            })->size();
        
        if ( motif_size < 7)
        {
            TR_type = "STR";
        }
        else
        {
            TR_type = "VNTR";
        }

        std::tie(chr, start, end) = parse_region(this->TR_regions,j);
        start_pos = (std::stoi(start) - this->padding > 0) ? std::stoi(start) - this->padding : 1; 
        end_pos = std::stoi(end) + this->padding;
        if (start_pos == end_pos || end_pos < start_pos ) {
            if (this->verbose == 2){
                spdlog::warn("WARNING: The start and end positions are the same in region " + chr + ":" + start + "-" + end + " skipping the region...");
            }
            continue;
        }
        uint16_t CN_reference_all_runs=0;
        std::vector<std::pair<uint32_t, uint32_t>> tandem_runs_positions = {}; // this is a vector of pairs of start and end positions of the tandem runs
        uint16_t buffer = 0; // this is a buffer for the start position of the tandem runs, if the start position is 0 or 1 then the buffer will be 1, otherwise it will be 2, to make sure to include ins/del in the first position of the tandem run
        std::string updated_region = chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos);
        std::vector<Interval> path_reference  = genotype_reference(updated_region,  std::get<1>(this->TR_regions[j]), seq, ref_seq_len, fai, lcp);
        std::vector<std::vector<Interval>> tandem_runs = {};

        if (path_reference.size() == 0) {
            spdlog::warn("unexpected behavior: no motif was found in the reference sequence for region, consider changing the mml threshold " + chr + ":" + start + "-" + end);
            continue;
        }

        if (motif_size <= this->motif_threashold_size && this->refineTrRegions) { 
            tandem_runs = findAllTandemRuns(path_reference,motif_size);
            tandem_runs_positions.reserve(tandem_runs.size());
        }
        else {
            tandem_runs = {path_reference};
            if (!this->refineTrRegions){
                CN_reference_all_runs = path_reference.size();
            }
            tandem_runs_positions = {std::make_pair(start_pos, end_pos)};
        }
        if (this->refineTrRegions){
            tandem_runs_positions.clear();
            buffer = tandem_runs[0][0].start < 2 ? 1 : 2;
            for (auto& tandem_run : tandem_runs) {
                // print_dashes_representation(std::get<0>(this->TR_regions[j]).c_str(), seq, tandem_run);
                CN_reference_all_runs += tandem_run.size();
                tandem_runs_positions.emplace_back(std::make_pair(tandem_run[0].start + start_pos - buffer, tandem_run.back().end + start_pos -1 ));
            }
        }

        tandem_runs.clear();
        iter = sam_itr_querys(idx,h, updated_region.c_str());
        if (iter == NULL) {
                std::cerr << "Failed to fetch iterator for region " << std::get<0>(this->TR_regions[j]) << std::endl;
        }


        std::string result_region = "";
        std::string region = std::get<0>(this->TR_regions[j]);
        vcfRecordInfoReads recordInfo("", {}, {}, {}, {}, {},{},{}, "", "", "", 0, 0, 0, {});
        std::vector<std::string> ref_cord = {};
        std::string reference_seq = seq;
        std::vector<std::vector<uint16_t>> motif_ids_alleles = {};
        std::vector<std::vector<std::string>> motif_cords_alleles = {};
        std::vector<std::uint16_t> cluster_sizes = {};
        seq.clear();
        regionResult region_results = processLongReads(fp, iter, start_pos, end_pos, region,  motifs,motif_size, chr,tandem_runs_positions,region_phasing_result,cut_reads_fasta,lcp);//
        std::vector<uint16_t> motif_ids_reference = get_motif_ids(path_reference);
        for(const auto &interval : path_reference){
                ref_cord.push_back("(" + std::to_string(interval.start) + "-" + std::to_string(interval.end) + ")");
        }
  
        if (region_results.clusters[0].empty() && region_results.clusters[1].empty()) {
            recordInfo = vcfRecordInfoReads(region, motifs, {{},{}}, motif_ids_reference, {{},{}}, ref_cord, {0,0},{"",""}, reference_seq, TR_type, chr, static_cast<uint32_t>(start_pos), 0, 0, {0,0});
        }
        else {
            ref_seq_len = 0;
            chunk_phasing_results.push_back(region_phasing_result);
            region_phasing_result.clear();
            for (uint16_t i = 0; i < region_results.clusters.size(); ++i) {
                auto &cluster = region_results.clusters[i];
                if (cluster.empty()) {
                    cluster_sizes.push_back(0);
                    motif_ids_alleles.push_back({});
                    motif_cords_alleles.push_back({});
                    continue;
                }
                cluster_sizes.push_back(cluster.size());
                std::vector<uint16_t> motif_ids_allele = {};
                std::vector<std::string> motif_ids_allele_cord = {};
                for (const auto &interval : cluster[region_results.indices_consensus_read[i]].second) {
                    motif_ids_allele.push_back(interval.motif_id);
                    motif_ids_allele_cord.push_back("(" + std::to_string(interval.start) + "-" + std::to_string(interval.end) + ")");
                }
                motif_ids_alleles.push_back(motif_ids_allele);
                motif_cords_alleles.push_back(motif_ids_allele_cord);
            }
            if (motif_ids_alleles.size() == 1) {
                motif_ids_alleles.push_back(motif_ids_alleles[0]);
                motif_cords_alleles.push_back(motif_cords_alleles[0]);
            }

            recordInfo = vcfRecordInfoReads(region, motifs,motif_ids_alleles,  motif_ids_reference, motif_cords_alleles, ref_cord, region_results.copy_numbers, region_results.alt_seqs,  reference_seq, TR_type, chr, static_cast<uint32_t>(start_pos), CN_reference_all_runs, region_results.num_reads_covering_TR,cluster_sizes);
        }

        hts_itr_destroy(iter);
        chunk_genotype_records.push_back(recordInfo);
        tandem_runs_positions.clear();
        CN_reference_all_runs = 0;

    }
    return std::make_tuple(chunk_phasing_results,chunk_genotype_records);
}



void TandemTwister::processRegionsForLongReadsInput() {
    /**
     * Process regions for long reads input
    */ 
    std::vector<pid_t> pids;

    uint32_t num_processes = std::min(this->num_threads, static_cast<uint32_t>(this->TR_regions.size()));
    uint32_t num_of_regions = this->TR_regions.size() / num_processes;
    for (size_t  i = 0; i < num_processes; ++i) {
        pid_t pid = fork();
        if (pid == 0) {

            uint32_t start = i * num_of_regions;
            uint32_t end = (i == num_processes - 1) ? this->TR_regions.size() : (i + 1) * num_of_regions;
            std::vector <std::string> cut_reads_chunk;
            std::tuple<std::vector<std::string> ,std::vector<vcfRecordInfoReads>>  results_chunk;
            samFile* fp = NULL;
            bam_hdr_t* h = NULL;
            hts_idx_t* idx = NULL;
            faidx_t* fai = NULL;
            open_bam_file(this->input_bam, fp, h, idx);
            open_reference_file(this->input_reference, fai);

            try{
                std::string cut_reads_fasta = "";
                results_chunk = process_chunk_reads(start, end,fp,h,idx,fai,cut_reads_fasta);
                // open a temp file to write the cut reads to
                if (this->keep_cut_sequence){
                    std::string cut_reads_file = this->output_path + "cut_reads_" + this->sampleName + "_" + this->reads_type + "_" + std::to_string(i) + ".fasta";
                    std::ofstream outfile_cutReads(cut_reads_file);
                    if (!outfile_cutReads.is_open()) {
                        std::cerr << "Failed to open file " << cut_reads_file << std::endl;
                        exit(1);
                    }
                    outfile_cutReads << cut_reads_fasta;
                    outfile_cutReads.flush();
                    outfile_cutReads.close();
                }
                
                if (this->keep_phasing_results){
                    std::string phasing_file = this->output_path + "phasing_" + this->sampleName + "_" + this->reads_type + "_" + std::to_string(i) + ".tsv";
                    std::ofstream outfile_phasing(phasing_file);
                    if (!outfile_phasing.is_open()) {
                        std::cerr << "Failed to open file " << phasing_file << std::endl;
                        exit(1);
                    }
                    for (const std::string &phasing_result : std::get<0>(results_chunk)) {

                        outfile_phasing << phasing_result;
               
                    }
                    outfile_phasing.flush();
                    outfile_phasing.close();
                }

                std::string vcf_file = this->output_path  + this->sampleName + "_" + this->reads_type + "_" + std::to_string(i) + ".vcf";
                // print the results_chunk
              
                writeRecordsToVcf(std::get<1>(results_chunk),vcf_file);
                outfile.close();
                hts_close(fp);
                bam_hdr_destroy(h);
                hts_idx_destroy(idx);
                fai_destroy(fai);
                exit(0);
            }
            catch (const std::runtime_error& e) {
                std::cerr << e.what() << std::endl;
                exit(1);
            }          
        }
        else if (pid < 0) {
            std::cerr << "Error: could not fork" << std::endl;
            return;
        }
        else {
            pids.push_back(pid);
        }
    }
    for (auto pid : pids) {
        int status;
        waitpid(pid, &status, 0);
    }
    

    std::string cut_reads_fasta = "";
    std::string genotype_result = "";
    std::string phasing_result = "";
    std::string vcf_result = "";

    // open the main vcf file
    htsFile* outfile_vcf = bcf_open(this->output_file_vcf.c_str(), "wz");
    // write the this->vcf_header to the vcf file

    if (outfile_vcf == NULL) {
        std::cerr << "Failed to open file " << this->output_file_vcf << std::endl;
        exit(1);
    }
    if (bcf_hdr_write(outfile_vcf, this->vcf_header) != 0) {
        fprintf(stderr, "Error writing VCF header.\n");
        return;
    }

    
    if (outfile_vcf == NULL) {
        std::cerr << "Failed to open file " << vcf_result << std::endl;
        exit(1);
    }
    
    
    for (size_t i = 0; i < num_processes; ++i) {
        std::string cut_reads_file_chunk = this->output_path + "cut_reads_" + this->sampleName + "_"+ this->reads_type + "_" + std::to_string(i) + ".fasta";
        std::ifstream infile_cutReads(cut_reads_file_chunk);
        if (this->keep_cut_sequence){
            if (!infile_cutReads.is_open()) {
                std::cerr << "Failed to open file " << cut_reads_file_chunk << std::endl;
                exit(1);
            }
        } 

        std::string phasing_file_chunk = this->output_path + "phasing_" + this->sampleName + "_"+ this->reads_type + "_" + std::to_string(i) + ".tsv";
        std::ifstream infile_phasing(phasing_file_chunk);
        if (this->keep_phasing_results){

        if (!infile_phasing.is_open()) {
            std::cerr << "Failed to open file " << phasing_file_chunk << std::endl;
            exit(1);
        }

        }

        
        std::string vcf_file_chunk = this->output_path + this->sampleName + "_" + this->reads_type + "_" + std::to_string(i) + ".vcf";
        htsFile* infile_vcf = hts_open(vcf_file_chunk.c_str(), "r");
       
        if (infile_vcf == NULL) {
            std::cerr << "Failed to open file " << vcf_file_chunk << std::endl;
            exit(1);
        }
        
        bcf_hdr_t *hdr  = bcf_hdr_read(infile_vcf);
        if (!hdr) {
            std::cerr << "Error: Unable to read header of the VCF file " << vcf_file_chunk << std::endl;

            hts_close(infile_vcf);
            return ;
        }
        
        // verify the headers match
        if (bcf_hdr_nsamples(this->vcf_header) != bcf_hdr_nsamples(hdr)) {
            std::cerr << "Warning: Input and output headers may be incompatible" << std::endl;
        }

        bcf1_t *rec = bcf_init();
        while (bcf_read(infile_vcf, hdr, rec) == 0) {
            bcf_unpack(rec, BCF_UN_ALL);
            // print the record
            if (bcf_write(outfile_vcf,  this->vcf_header, rec) < 0) {
                //std::cerr << "Error writing VCF record to file." << std::endl;
                // skip the record
                //std::cerr << "Skipping the record" << std::endl;
                continue;
            }
        }
        bcf_destroy(rec);
        hts_close(infile_vcf);

        std::stringstream buffer;
        buffer << infile_cutReads.rdbuf();
        cut_reads_fasta += buffer.str();
        std::stringstream buffer2;
        // buffer2 << infile.rdbuf();
        genotype_result += buffer2.str();
        std::stringstream buffer3;
        buffer3 << infile_phasing.rdbuf();
        phasing_result += buffer3.str();
        std::stringstream buffer4;
        infile_cutReads.close();
        // infile.close();
        infile_phasing.close();
        std::remove(cut_reads_file_chunk.c_str());
        //std::remove(genotype_file_chunk.c_str());
        std::remove(phasing_file_chunk.c_str());
        std::remove(vcf_file_chunk.c_str());
        bcf_hdr_destroy(hdr);
    }
    
    // gzip the vcf file
    if (this->keep_cut_sequence){
        this->outfile_cutReads << cut_reads_fasta;
        outfile_cutReads.flush();
        outfile_cutReads.close();
    }
    this->outfile << genotype_result;
    this->outfile_phasing << phasing_result;

    this->outfile_phasing.flush();
    this->outfile.flush();

    this->outfile.close();
    this->outfile_phasing.close();
    
    bcf_hdr_destroy(this->vcf_header);
    bcf_close(outfile_vcf);

    

    // index the vcf file
    std::string vcf_file_index = this->output_file_vcf + ".tbi";
    if (bcf_index_build(this->output_file_vcf.c_str(), 0) != 0) {
        std::cerr << "Failed to index the vcf file " << this->output_file_vcf << std::endl;
        exit(1);
    }
   
}   
