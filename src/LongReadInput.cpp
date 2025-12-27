#include "include/TandemTwister.hpp"
#include <thread>
#include <mutex>
#include <iostream>
#include <chrono>
#include <limits>
#include <array>


namespace {
    /*
    * This namespace is used to colorize the motifs for debugging
    * @param: none
    * @return: none
    */
#ifdef _WIN32
    constexpr const char* ANSI_RESET = "";
    std::string colorizeSegment(const std::string& content, const char* /*color_code*/) {
        return content;
    }
    const char* motifColorFor(uint32_t /*motif_id*/) {
        return "";
    }
#else
    constexpr const char* ANSI_RESET = "\033[0m";
    std::string colorizeSegment(const std::string& content, const char* color_code) {
        if (color_code == nullptr || *color_code == '\0') {
            return content;
        }
        return std::string(color_code) + content + ANSI_RESET;
    }
    const char* motifColorFor(uint32_t motif_id) {
        static const std::array<const char*, 8> palette = {
            "\033[38;5;45m",
            "\033[38;5;214m",
            "\033[38;5;99m",
            "\033[38;5;46m",
            "\033[38;5;160m",
            "\033[38;5;226m",
            "\033[38;5;27m",
            "\033[38;5;201m"
        };
        return palette[motif_id % palette.size()];
    }
#endif

    void visualizeMotifRunsDebug(const std::vector<Interval>& intervals,
                                 const std::vector<std::vector<Interval>>& keptRuns,
                                 const std::vector<std::vector<Interval>>& removedRuns,
                                 const std::vector<std::string>& motifs,
                                 const std::string& sequence) {
        /*
        * This function is used to visualize the motif runs for debugging
        * @param intervals: The intervals to visualize
        * @param keptRuns: The kept runs
        * @param removedRuns: The removed runs
        * @param motifs: The motifs to visualize
        * @param sequence: The sequence to visualize
        * @return: none
        */
        if (sequence.empty()) {
            spdlog::debug("visualizeMotifRunsDebug: sequence is empty");
            return;
        }

        struct Annotation {
            size_t start;
            size_t end;
            uint32_t motif_id;
            bool kept;
        };

        std::vector<Annotation> annotations;
        annotations.reserve(intervals.size());
        std::string statusLine(sequence.size(), '.');

        auto recordRuns = [&](const std::vector<std::vector<Interval>>& runs, bool kept) {
            const char mark = kept ? '+' : 'X';
            for (const auto& run : runs) {
                for (const auto& interval : run) {
                    if (interval.end < interval.start) {
                        continue;
                    }
                    size_t startIdx = (interval.start > 0) ? static_cast<size_t>(interval.start - 1) : 0;
                    if (startIdx >= sequence.size()) {
                        continue;
                    }
                    size_t endIdx = static_cast<size_t>(interval.end);
                    endIdx = std::min(endIdx, sequence.size());
                    if (endIdx <= startIdx) {
                        continue;
                    }
                    annotations.push_back({startIdx, endIdx, interval.motif_id, kept});
                    for (size_t pos = startIdx; pos < endIdx; ++pos) {
                        statusLine[pos] = mark;
                    }
                }
            }
        };

        recordRuns(keptRuns, true);
        recordRuns(removedRuns, false);

        std::sort(annotations.begin(), annotations.end(), [](const Annotation& lhs, const Annotation& rhs) {
            return lhs.start < rhs.start;
        });

        std::string coloredSequence;
        coloredSequence.reserve(sequence.size() * 2);
        size_t cursor = 0;

        for (const auto& annotation : annotations) {
            if (cursor < annotation.start) {
                coloredSequence.append(sequence.substr(cursor, annotation.start - cursor));
            }
            const char* color_code = motifColorFor(annotation.motif_id);
            coloredSequence += colorizeSegment(sequence.substr(annotation.start, annotation.end - annotation.start), color_code);
            cursor = annotation.end;
        }
        if (cursor < sequence.size()) {
            coloredSequence.append(sequence.substr(cursor));
        }

        std::string coloredMotifsList = "{";
        for (size_t i = 0; i < motifs.size(); ++i) {
            const char* motif_color = motifColorFor(static_cast<uint32_t>(i));
            coloredMotifsList += colorizeSegment(motifs[i], motif_color);
            if (i + 1 < motifs.size()) {
                coloredMotifsList += ", ";
            }
        }
        coloredMotifsList += "}";

        spdlog::debug("findAllTandemRuns: motifs {}", coloredMotifsList);
        spdlog::debug("findAllTandemRuns: sequence view {}", coloredSequence);
        spdlog::debug("findAllTandemRuns: status line   {}", statusLine);
        spdlog::debug("findAllTandemRuns: legend '+' kept, 'X' removed, '.' untouched");
    }
}


bool TandemTwister::IsMotifRun(const std::vector<Interval> & intervals, const std::vector<uint16_t>& motif_lengths){
    /**
    * Check if the intervals are motif runs
    * @param intervals: The intervals to check
    * @param motif_lengths: The size of the motifs indexed by motif id
    * @return: True if the intervals are motif runs, false otherwise
    */
    if (intervals.empty()) {
        spdlog::debug("IsMotifRun: rejecting empty interval run");
        return false;
    }

    if (spdlog::should_log(spdlog::level::debug)) {
        spdlog::debug("IsMotifRun: evaluating run [{} intervals] span {}-{}", intervals.size(), intervals.front().start, intervals.back().end);
    }

    auto motifLengthFor = [&motif_lengths](const Interval& interval) -> uint16_t {
        if (interval.motif_id < motif_lengths.size()) {
            return motif_lengths[interval.motif_id];
        }
        if (interval.end >= interval.start) {
            uint32_t span = interval.end - interval.start + 1;
            return span > std::numeric_limits<uint16_t>::max() ? std::numeric_limits<uint16_t>::max() : static_cast<uint16_t>(span);
        }
        return 0;
    };

    uint32_t maxMotifSpan = 0;
    uint16_t maxMotifLenInRun = 0;
    bool hasKnownLength = false;
    for (const auto& interval : intervals) {
        uint16_t motifLen = motifLengthFor(interval);
        if (motifLen > 0) {
            hasKnownLength = true;
        }
        maxMotifSpan += motifLen;
        maxMotifLenInRun = std::max<uint16_t>(maxMotifLenInRun, motifLen);
    }

    if (!hasKnownLength) {
        spdlog::debug("IsMotifRun: rejecting run {}-{} due to missing motif length mapping", intervals.front().start, intervals.back().end);
        return false;
    }

    uint32_t intervalSpan = (intervals.back().end >= intervals.front().start)
        ? intervals.back().end - intervals.front().start + 1
        : 0;
    if (intervalSpan == 0) {
        spdlog::debug("IsMotifRun: rejecting run with zero span ({}-{})", intervals.front().start, intervals.back().end);
        return false;
    }

    uint32_t intervalSpan2 = std::min(maxMotifSpan, intervalSpan);
    
    bool isHomopolymer = maxMotifLenInRun <= 1;
    if (isHomopolymer && intervalSpan2 < 4){ // for homopolymers we only consider motif runs that are at least 4 bases long
        spdlog::debug("IsMotifRun: rejecting homopolymer run {}-{} span {} <4", intervals.front().start, intervals.back().end, intervalSpan2);
        return false;
    }
    if (!isHomopolymer && intervalSpan2 < 8){
        spdlog::debug("IsMotifRun: rejecting mixed-motif run {}-{} span {} <8", intervals.front().start, intervals.back().end, intervalSpan2);
        return false;
    }

    uint32_t actualmotifspan = 0;
    const float threshold = 0.7f;
    if (intervals[0].start == intervals[0].end){
        actualmotifspan = intervals.size();
    }
    else{
        actualmotifspan = std::accumulate(intervals.begin(), intervals.end(), 0u, 
            [](uint32_t sum, const Interval& interval) {
                return sum + interval.end - interval.start + 1;
            });
    }
    float purity = static_cast<float>(actualmotifspan) / static_cast<float>(intervalSpan);
    bool isMotif = purity > threshold;
    spdlog::debug("IsMotifRun: run {}-{} purity {:.3f} (threshold {:.2f}) -> {}", intervals.front().start, intervals.back().end, purity, threshold, isMotif ? "accept" : "reject");
    return isMotif;
}


std::vector<std::vector<Interval>> TandemTwister::findAllTandemRuns(const std::vector<Interval>& intervals, const std::vector<std::string>& motifs, const std::string& sequence) {
    /**
     * Find all tandem runs in the intervals
     * @param intervals: The intervals to search for tandem runs
     * @param motifs: The motifs observed in the region
     * @return: A vector of vectors containing the tandem runs
     */
    if (intervals.empty()) {
        spdlog::debug("findAllTandemRuns: received empty interval list");
        return std::vector<std::vector<Interval>>();
    }

    std::vector<uint16_t> motif_lengths;
    motif_lengths.reserve(motifs.size());
    for (const auto& motif : motifs) {
        motif_lengths.push_back(static_cast<uint16_t>(motif.size()));
    }

    if (spdlog::should_log(spdlog::level::debug)) {
        uint16_t minLen = motif_lengths.empty() ? 0 : *std::min_element(motif_lengths.begin(), motif_lengths.end());
        uint16_t maxLen = motif_lengths.empty() ? 0 : *std::max_element(motif_lengths.begin(), motif_lengths.end());
        spdlog::debug("findAllTandemRuns: evaluating {} intervals across {} motifs (len range {}-{})", intervals.size(), motifs.size(), minLen, maxLen);
    }

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
        uint16_t currentMotifLength = std::max<uint16_t>(
            (currentInterval.motif_id < motif_lengths.size()) ? motif_lengths[currentInterval.motif_id] : 1,
            static_cast<uint16_t>(1));
        uint16_t previousMotifLength = std::max<uint16_t>(
            (previousInterval.motif_id < motif_lengths.size()) ? motif_lengths[previousInterval.motif_id] : 1,
            static_cast<uint16_t>(1));
        uint32_t threshold = this->tanCon * std::max<uint16_t>(currentMotifLength, previousMotifLength) + 1;
        int gap = static_cast<int>(currentInterval.start) - static_cast<int>(previousInterval.end);
        int absGap = std::abs(gap);
        if (absGap > static_cast<int>(threshold)) {
            spdlog::debug("findAllTandemRuns: splitting at idx {} gap {} (thr {}) prev motif {} curr motif {}", i, absGap, threshold, previousInterval.motif_id, currentInterval.motif_id);
            allRuns.push_back(std::vector<Interval>(intervals.begin() + currentRunStart, intervals.begin() + i));
            currentRunStart = i;
        }
    }
    allRuns.push_back(std::vector<Interval>(intervals.begin() + currentRunStart, intervals.end()));
    
    spdlog::debug("findAllTandemRuns: {} candidate runs before filtering", allRuns.size());

    std::vector<std::vector<Interval>> removedRuns;
    std::vector<std::vector<Interval>> keptRuns;
    removedRuns.reserve(allRuns.size());
    keptRuns.reserve(allRuns.size());

    for (auto& run : allRuns) {
        if (run.empty() || run.size() == 1 || !IsMotifRun(run , motif_lengths)) {
            removedRuns.push_back(run);
        } else {
            keptRuns.push_back(run);
        }
    }

    if (spdlog::should_log(spdlog::level::debug)) {
        visualizeMotifRunsDebug(intervals, keptRuns, removedRuns, motifs, sequence);
    }

    allRuns = std::move(keptRuns);

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
        spdlog::error("Failed to fetch sequence for region {}", regionWithPadding);
        return std::vector<Interval>();
    }
    //std::unordered_map<char, uint16_t> nucleotides_occurences;
    std::vector<Interval> path_reference = {};
   
    path_reference = findBestPath(seq, motifs, this->match_score, this->mismatch_score, this->gap_score, lcp);
    return path_reference;
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




std::tuple<uint32_t, uint32_t, uint16_t> TandemTwister::findMostOccurringCopyNumber(const std::vector<std::tuple<std::string, std::vector<Interval>, std::string>>& cluster) {
    /**
     * Find the most occurring copy number in the cluster
     * @param cluster: The cluster to search for the most occurring copy number
     * @return: A tuple containing the most occurring copy number, the frequency of the copy number, and the index of the most occurring copy number
     */

    if (cluster.empty()) {
        return std::make_tuple(0, 0, 0);
    }
    std::unordered_map<std::string, uint32_t> frequencyMap;



    std::string temp_feature = "";
    for (const auto& tuple : cluster) {
        temp_feature = std::get<2>(tuple);
        ++frequencyMap[temp_feature];
    }


    // Visualize the frequencyMap for debugging/inspection
    spdlog::debug("Copy number frequency map visualization:");
    for (const auto& entry : frequencyMap) {
        spdlog::debug("Copy number: {}, Frequency: {}", entry.first, entry.second);
    }


    // if the size of the cluster is 2 then check if the CN of both is the same then return the index of the first read, if not return the index of the maximum CN
    if (cluster.size() == 2) {
        if (std::get<1> (cluster[0]).size() == std::get<1>(cluster[1]).size()) {
            return std::make_tuple(std::get<1>(cluster[0]).size(), 2, 0);
        }
        else {
            return std::make_tuple(std::get<1>(cluster[0]).size() > std::get<1>(cluster[1]).size() ? std::get<1>(cluster[0]).size() : std::get<1>(cluster[1]).size(), 1, std::get<1>(cluster[0]).size() > std::get<1>(cluster[1]).size() ? 0 : 1);
        }
    }

    // if there are more than 2 copy numbers with the same frequency then return the index of the highest one
    auto maxElement = std::max_element(frequencyMap.begin(), frequencyMap.end(),
        [](const std::pair<std::string, uint32_t>& a, const std::pair<std::string, uint32_t>& b) {
            return a.second < b.second;
        });
    spdlog::debug("maxElement: copyNumber = {}, frequency = {}", maxElement->first, maxElement->second);

    // get the second most occuring copy number
    auto secondMaxElement = frequencyMap.end();
    uint32_t secondMaxFreq = 0;
    for (const auto& entry : frequencyMap) {
        if (entry.first != maxElement->first && entry.second > secondMaxFreq) {
            secondMaxFreq = entry.second;
            secondMaxElement = frequencyMap.find(entry.first);
        }
    }

    if (secondMaxElement != frequencyMap.end()) {
        spdlog::debug("Second max element: copyNumber = {}, frequency = {}", secondMaxElement->first, secondMaxElement->second);
    } else {
        spdlog::debug("Second max element not found.");
    }
    auto sumCopyNumberFromFeatureString = [](const std::string& featureStr) -> uint32_t {
        uint32_t sum = 0;
        std::istringstream iss(featureStr);
        std::string token;
        while (iss >> token) {
            try {
                sum += static_cast<uint32_t>(std::stoul(token));
            } catch (const std::exception& e) {
                spdlog::warn("Failed to convert part of feature string '{}' to uint32_t: {}", token, e.what());
            }
        }
        return sum;
    };

    auto maxElement_cn = sumCopyNumberFromFeatureString(maxElement->first);
    
    std::string mostOccurringCopyNumber_str =  maxElement->first;
    uint32_t mostOccurringCopyNumber_int = maxElement_cn;
    uint32_t frequency = maxElement->second;
    uint16_t index = 0;

    // if the second most occuring copy number has the same frequency as the most occuring copy number then return the index of the highest copy number
    if (secondMaxElement != frequencyMap.end() && secondMaxElement->second == maxElement->second && frequencyMap.size() > 1) {
        auto secondMaxElement_cn = sumCopyNumberFromFeatureString(secondMaxElement->first);
        spdlog::debug("Second max element: copyNumber = {}, frequency = {}", secondMaxElement->first, secondMaxElement->second);
        auto maxCopyNumber = std::max(maxElement_cn, secondMaxElement_cn);
        if (maxCopyNumber != maxElement_cn){
            mostOccurringCopyNumber_str = secondMaxElement->first;
            frequency = secondMaxElement ->second;
            mostOccurringCopyNumber_int = secondMaxElement_cn;
        }
    }

    spdlog::debug("mostOccurringCopyNumber compostion is:  ({})", mostOccurringCopyNumber_str);

    for (size_t i = 0; i < cluster.size(); ++i) {
        if (std::get<2>(cluster[i]) == mostOccurringCopyNumber_str) {
            index = i;
            break;
        }
    }
    spdlog::debug("Index of the mostOccusring copy number is {}: " ,index );
    return std::make_tuple(mostOccurringCopyNumber_int, frequency, index);
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



void TandemTwister::cluster_by_features(const std::string& chr, 
                                        std::vector<std::tuple<std::string, std::vector<Interval>, uint8_t>>& reads_intervals, 
                                        std::vector<std::vector<std::tuple<std::string, std::vector<Interval>, std::string>>> & clusters, 
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
            cn_1_avg += std::get<1>(read).size();
        }
        cn_1_avg /= clusters[0].size();

        for (const auto& read : clusters[1]) {
            cn_2_avg += std::get<1>(read).size();;
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

std::string get_read_features_str(const arma::rowvec& read_features) {
    /*
    * This function is used to get the read features as a string
    * @param read_features: The read features
    * @return: The read features as a string
    *   
    */
    std::string read_features_str = "";
    for (uint16_t i = 4; i < read_features.n_cols; ++i) {
        read_features_str += fmt::format("{} ", read_features(i));
    }
    return read_features_str;
}


regionResult TandemTwister::processLongReads(samFile* &fp, hts_itr_t* &iter,const size_t start_pos, const size_t end_pos, const std::string& region, 
                                    std::vector<std::string>& motifs, uint16_t motif_size, ::string &chr,
                                    std::vector<std::pair<uint32_t, uint32_t>> &tandem_runs,
                                    std::string &cut_reads_fasta, const std::vector<uint16_t>& lcp) {
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
     * @param num_reads: The number of reads
     * @param cut_reads_fasta: The cut reads in FASTA format
     * @param lcp: The longest common prefix array
     */

    // define the variables that will be outputted
    arma::mat region_features; 
    std::vector<std::vector<std::tuple<std::string, std::vector<Interval>, std::string>>> clusters;
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

        if (path_total.empty()) {
            // remove the read from the reads_sequences
            reads_sequences.erase(read_name);
            continue;
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
                // Use the index of the read in reads_intervals as the row index
                size_t read_index = &read - &reads_intervals[0];
                auto read_features = region_features.row(read_index);
                std::string read_features_str = get_read_features_str(read_features);
                if (std::get<2>(read) == 1){
                    // skip the first 4 features (cut_read_length, CN_total, path_purity, path_total_score)
                    clusters[0].push_back(std::make_tuple(std::get<0>(read),std::get<1>(read),read_features_str));
                }
                else if (std::get<2>(read) == 2){
                    clusters[1].push_back(std::make_tuple(std::get<0>(read),std::get<1>(read),read_features_str));
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
                size_t read_index = &read - &reads_intervals[0];
                auto read_features = region_features.row(read_index);
                std::string read_features_str = get_read_features_str(read_features);
                if (std::get<2>(read) == 1 || std::get<2>(read) == 2){
                    clusters[0].push_back(std::make_tuple(std::get<0>(read), std::get<1>(read), read_features_str));
                }
                else {
                    noise_cluster.push_back(std::make_pair(std::get<0>(read), std::get<1>(read)));
                }
            }
        }
        else {
            spdlog::debug("Clustering by features as no haplotypes tags were found in the region.");
            cluster_by_features(chr, reads_intervals, clusters, noise_cluster, region_features, motif_size, haplody);
        }
    }

    else{
        spdlog::debug("Clustering by features (input bam is not tagged).");
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
                    copy_numbers.push_back(std::get<1>(cluster[0]).size());
                    alt_seqs.push_back(reads_sequences[std::get<0>(cluster[0])]);
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
            std::string cut_seq = reads_sequences[std::get<0>(cluster[0])].substr(std::get<1>(cluster[0])[0].start - 1, std::get<1>(cluster[0]).back().end - std::get<1>(cluster[0])[0].start + 1);
            
            // reset the intervals starts and ends by subtracting the start of the first interval (this is to make the intervals start from 0) after removing the intervals that didn't pass the consensus range
            auto start_interval = std::get<1>(cluster[0])[0].start;
            for (size_t j = 0; j < std::get<1>(cluster[0]).size(); ++j) {

                std::get<1>(cluster[0])[j].start -= start_interval - 1;
                std::get<1>(cluster[0])[j].end -= start_interval - 1;

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
                alt_seqs.push_back(reads_sequences[std::get<0>(cluster[index_cn])]);
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
    


    region_results = regionResult{clusters, noise_cluster, alt_seqs, indices_consensus_read, copy_numbers, cn_occurences_in_cluster, num_reads_covering_TR};
    reads_intervals.clear();
    reads_sequences.clear();
    return region_results;
}



std::vector<vcfRecordInfoReads> TandemTwister::process_chunk_reads(size_t start_dict, size_t end_dict, samFile* fp, bam_hdr_t* h, hts_idx_t* idx, faidx_t* fai,std::string & cut_reads_fasta) {
    /**
     * Process chunk reads
     * @param start_dict: The start dictionary
     * @param end_dict: The end dictionary
     * @param fp: The file pointer
     * @param h: The header
     * @param idx: The index
     * @param fai: The FASTA index
     * @param cut_reads_fasta: The cut reads in FASTA format
     * @return: A vector of genotype records
     */

    std::vector<Interval> path_reference;
    std::string chr, start, end;
    std::string region;
    uint32_t  start_pos = 0; // this is the start position of the region to extract reads from
    uint32_t end_pos= 0; // this is the end position of the region to extract reads from
    uint32_t ref_seq_len = 0; // length of the reference sequence
    std::string  seq = "";
    // std::string region_genotype_result; // this is the genotype result for the region in Bed format
    std::vector<std::string> chunk_genotype_results; // this is a vector of genotype results for the regions for each process
    std::vector<vcfRecordInfoReads> chunk_genotype_records; // this is a vector of genotype results for the regions for each process
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
            spdlog::debug("Finding all tandem runs");
            tandem_runs = findAllTandemRuns(path_reference, motifs, seq);
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
                CN_reference_all_runs += tandem_run.size();
                tandem_runs_positions.emplace_back(std::make_pair(tandem_run[0].start + start_pos - buffer, tandem_run.back().end + start_pos -1 ));
            }
        }

        tandem_runs.clear();
        iter = sam_itr_querys(idx,h, updated_region.c_str());
        if (iter == NULL) {
                spdlog::error("Failed to fetch iterator for region {}", std::get<0>(this->TR_regions[j]));
        }

        // initialize the result region
        std::string result_region = "";
        // get the region
        std::string region = std::get<0>(this->TR_regions[j]);
        // initialize the recordInfo
        vcfRecordInfoReads recordInfo("", {}, {}, {}, {}, {},{},{}, "", "", "", 0, 0, 0, {});
        // initialize the reference coordinates
        std::vector<std::string> ref_cord = {};
        // initialize the reference sequence
        std::string reference_seq = seq;
        // initialize the motif ids for the alleles
        std::vector<std::vector<uint16_t>> motif_ids_alleles = {};
        // initialize the motif coordinates for the alleles
        std::vector<std::vector<std::string>> motif_cords_alleles = {};
        // initialize the cluster sizes
        std::vector<std::uint16_t> cluster_sizes = {};
        // clear the sequence
        seq.clear();
        // process the long reads
        regionResult region_results = processLongReads(fp, iter, start_pos, end_pos, region,  motifs,motif_size, chr,tandem_runs_positions,cut_reads_fasta,lcp);//
        // get the motif ids for the reference
        std::vector<uint16_t> motif_ids_reference = get_motif_ids(path_reference);
        // get the reference coordinates
        for(const auto &interval : path_reference){
                ref_cord.push_back("(" + std::to_string(interval.start) + "-" + std::to_string(interval.end) + ")");
        }
        // if the clusters are empty then set the recordInfo to the empty record
        if (region_results.clusters[0].empty() && region_results.clusters[1].empty()) {
            recordInfo = vcfRecordInfoReads(region, motifs, {{},{}}, motif_ids_reference, {{},{}}, ref_cord, {0,0},{"",""}, reference_seq, TR_type, chr, static_cast<uint32_t>(start_pos), 0, 0, {0,0});
        }
        else {
            ref_seq_len = 0;
            // iterate over the clusters
            for (uint16_t i = 0; i < region_results.clusters.size(); ++i) {
                auto &cluster = region_results.clusters[i];
                if (cluster.empty()) {
                    // if the cluster is empty then set the cluster sizes to 0
                    cluster_sizes.push_back(0);
                    motif_ids_alleles.push_back({});
                    motif_cords_alleles.push_back({});
                    continue;
                }
                // add the cluster size to the cluster sizes
                cluster_sizes.push_back(cluster.size());
                std::vector<uint16_t> motif_ids_allele = {};
                std::vector<std::string> motif_ids_allele_cord = {};
                // iterate over the intervals in the cluster to get the motif ids and coordinates
                for (const auto &interval : std::get<1>(cluster[region_results.indices_consensus_read[i]])) {
                    motif_ids_allele.push_back(interval.motif_id);
                    motif_ids_allele_cord.push_back("(" + std::to_string(interval.start) + "-" + std::to_string(interval.end) + ")");
                }
                // update the motif ids and coordinates for the alleles
                motif_ids_alleles.push_back(motif_ids_allele);
                motif_cords_alleles.push_back(motif_ids_allele_cord);
            }
            if (motif_ids_alleles.size() == 1) {
                // if there is only one cluster then set the motif ids and coordinates  for the second cluster to the first cluster
                motif_ids_alleles.push_back(motif_ids_alleles[0]);
                motif_cords_alleles.push_back(motif_cords_alleles[0]);
            }
            // return the recordInfo
            recordInfo = vcfRecordInfoReads(region, motifs,motif_ids_alleles,  motif_ids_reference, motif_cords_alleles, ref_cord, region_results.copy_numbers, region_results.alt_seqs,  reference_seq, TR_type, chr, static_cast<uint32_t>(start_pos), CN_reference_all_runs, region_results.num_reads_covering_TR,cluster_sizes);
        }

        hts_itr_destroy(iter);
        // add the recordInfo to the chunk genotype records
        chunk_genotype_records.push_back(recordInfo);
        // clear the tandem runs positions
        tandem_runs_positions.clear();
        // clear the CN reference all runs
        CN_reference_all_runs = 0;

    }
    return chunk_genotype_records;
}



void TandemTwister::processRegionsForLongReadsInput() {
    /**
     * Process regions for long reads input using multithreading
    */ 
    // Thread-safe data structures to collect results
    std::mutex results_mutex;
    std::vector<vcfRecordInfoReads> all_genotype_records;
    std::string all_cut_reads_fasta = "";
    
    // Track peak memory usage across all threads
    std::mutex memory_mutex;
    double peak_rss = 0.0;
    double peak_vm = 0.0;

    uint32_t num_threads = std::min(this->num_threads, static_cast<uint32_t>(this->TR_regions.size()));
    uint32_t num_of_regions = this->TR_regions.size() / num_threads;
    // Lambda function for thread worker
    auto thread_worker = [&](size_t thread_id) {
        uint32_t start = thread_id * num_of_regions;
        uint32_t end = (thread_id == num_threads - 1) ? this->TR_regions.size() : (thread_id + 1) * num_of_regions;
        
        samFile* fp = NULL;
        bam_hdr_t* h = NULL;
        hts_idx_t* idx = NULL;
        faidx_t* fai = NULL;
        
        try {
            open_bam_file(this->input_bam, fp, h, idx);
            open_reference_file(this->input_reference, fai);

            std::string cut_reads_fasta = "";
            std::vector<vcfRecordInfoReads> chunk_genotype_records = process_chunk_reads(start, end, fp, h, idx, fai, cut_reads_fasta);

            // Collect results thread-safely (file writing done in main thread)
            {
                std::lock_guard<std::mutex> lock(results_mutex);
                all_genotype_records.insert(all_genotype_records.end(), 
                                           chunk_genotype_records.begin(), 
                                           chunk_genotype_records.end());
                all_cut_reads_fasta += cut_reads_fasta;
            }
            
            // Measure and update peak memory usage (thread-safe)
            {
                unsigned long vsize;
                long rss_val;
                {
                    std::string ignore;
                    std::ifstream ifs("/proc/self/stat", std::ios_base::in);
                    ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                            >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                            >> ignore >> ignore >> vsize >> rss_val;
                }
                long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
                double vm = vsize / 1024.0;
                double rss = rss_val * page_size_kb;
                
                // Update peak memory usage thread-safely
                {
                    std::lock_guard<std::mutex> lock(memory_mutex);
                    if (rss > peak_rss) peak_rss = rss;
                    if (vm > peak_vm) peak_vm = vm;
                }
            }
       
            // Close file handles
            hts_close(fp);
            bam_hdr_destroy(h);
            hts_idx_destroy(idx);
            fai_destroy(fai);
        }
        catch (const std::runtime_error& e) {
            spdlog::error("Error in thread {}: {}", thread_id, e.what());
            // Clean up on error
            if (fp) hts_close(fp);
            if (h) bam_hdr_destroy(h);
            if (idx) hts_idx_destroy(idx);
            if (fai) fai_destroy(fai);
        }
    };
    
    // Launch threads
    std::vector<std::thread> threads;
    for (size_t i = 0; i < num_threads; ++i) {
        threads.emplace_back(thread_worker, i);
    }
    
    // Wait for all threads to complete
    for (auto& t : threads) {
        t.join();
    }
    
    // Log peak memory usage if verbose
    if (this->verbose >= 2 && peak_rss > 0) {
        spdlog::info("Peak memory usage during processing - VM: {:.2f} GB, RSS: {:.2f} GB", 
                     peak_vm / (1024 * 1024), peak_rss / (1024 * 1024));
    }

    std::string temp_vcf = this->output_path + "temp_" + this->sampleName + "_" + this->reads_type + ".vcf";
    writeRecordsToVcf(all_genotype_records, temp_vcf);
    
    // Open the main vcf file
    htsFile* outfile_vcf = bcf_open(this->output_file_vcf.c_str(), "wz");
    if (outfile_vcf == NULL) {
        spdlog::error("Failed to open file {}", this->output_file_vcf);
        return;
    }
    
    if (bcf_hdr_write(outfile_vcf, this->vcf_header) != 0) {
        fprintf(stderr, "Error writing VCF header.\n");
        bcf_close(outfile_vcf);
        return;
    }
    
    // Read records from temporary VCF file for sorting
    std::vector<std::tuple<std::string, int, int, bcf1_t*>> records_with_pos; // chr, pos, rid, record
    htsFile* temp_vcf_file = hts_open(temp_vcf.c_str(), "r");
    if (temp_vcf_file != NULL) {
        bcf_hdr_t *hdr = bcf_hdr_read(temp_vcf_file);
        if (hdr) {
            bcf1_t *rec = bcf_init();
            while (bcf_read(temp_vcf_file, hdr, rec) == 0) {
                bcf_unpack(rec, BCF_UN_ALL);
                
                // create a copy of the record
                bcf1_t *rec_copy = bcf_dup(rec);
                bcf_unpack(rec_copy, BCF_UN_ALL);
                
                // get the chromosome name
                std::string chrom = bcf_hdr_id2name(hdr, rec_copy->rid);
                // get the reference ID
                int rid = rec_copy->rid;
                
                // add the record to the records with position
                records_with_pos.emplace_back(chrom, rec_copy->pos, rid, rec_copy);
            }
            bcf_destroy(rec);
            bcf_hdr_destroy(hdr);
        }
        hts_close(temp_vcf_file);
    }
    
    // Remove temporary file
    std::remove(temp_vcf.c_str());
    // sort all records by chromosome and position
    std::sort(records_with_pos.begin(), records_with_pos.end(),
    [](const auto& a, const auto& b) {
        // compare by reference ID (rid) - this ensures correct chromosome order
        if (std::get<2>(a) != std::get<2>(b)) {
            return std::get<2>(a) < std::get<2>(b);
        }
        // compare by position
        return std::get<1>(a) < std::get<1>(b);
    });

    
    
    // write sorted records to final VCF
    for (auto& record_data : records_with_pos) {
        bcf1_t* rec = std::get<3>(record_data);
        if (bcf_write(outfile_vcf, this->vcf_header, rec) < 0) {
            spdlog::error("Error writing VCF record to file.");
        }
        bcf_destroy(rec);
    }
    // Write cut reads if requested
    if (this->keep_cut_sequence) {
        this->outfile_cutReads << all_cut_reads_fasta;
        this->outfile_cutReads.flush();
        this->outfile_cutReads.close();
    }
    
    this->outfile.flush();
    this->outfile.close();
    
    // Close the vcf file
    bcf_close(outfile_vcf);
    
    // Index the vcf file
    std::string vcf_file_index = this->output_file_vcf + ".tbi";
    if (bcf_index_build(this->output_file_vcf.c_str(), 0) != 0) {
        spdlog::error("Failed to index the vcf file {}", this->output_file_vcf);
    }
   
}   
                               