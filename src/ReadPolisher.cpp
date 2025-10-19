
#include "include/TandemTwister.hpp"

std::vector<Interval> updateIntervals(const std::string& alignedSeq,const std::string& seq, const std::vector<Interval>& intervals) {
    std::vector<Interval> intervals_aligned;
    intervals_aligned.reserve(intervals.size());
    std::vector<unsigned int> dashCounts(seq.size(), 0);
    unsigned int dashCount = 0;
    size_t seqIndex = 0;
  
    for (size_t i = 0; i < alignedSeq.size(); ++i) {
        if (alignedSeq[i] == '-') {
            dashCount++;
        } else if (seqIndex < seq.size()) {
            dashCounts[seqIndex] = dashCount;
            seqIndex++;
        }
      
    }
    // Adjust intervals based on the precomputed dash counts
    for (const Interval &interval : intervals) {
        Interval interval_aligned;
        if (interval.start >= dashCounts.size() || interval.end >= dashCounts.size()) {
            continue;
        }

        interval_aligned.start = interval.start + dashCounts[interval.start];
        interval_aligned.end = interval.end + dashCounts[interval.end];
        intervals_aligned.push_back(interval_aligned);
    }
    return intervals_aligned;
}


std::pair<std::vector<Interval> ,std::vector<Interval>> TandemTwister::needlemanWunsch(const std::string& seq1, const std::string& seq2, const AlignmentScores& scores,const std::vector<Interval> &ref_intevals, const std::vector<Interval> &read_intervals) {
    if (seq1.empty() || seq2.empty()) {
        return make_pair(ref_intevals, read_intervals);
    }
    int m = seq1.length();
    int n = seq2.length();
    // Initialize the matrix
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> backtrace(m + 1, std::vector<int>(n + 1, 0));

    // Fill the matrix and backtrace matrix
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match = dp[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? scores.match : scores.mismatch);
            int delete1 = dp[i - 1][j] - (j ==n ? 0 : scores.gap_extend) + scores.gap_open;
            int delete2 = dp[i][j - 1] - (i ==m ? 0 : scores.gap_extend) + scores.gap_open;

            // Update backtrace matrix
            if (match >= delete1 && match >= delete2) {
                backtrace[i][j] = 1;  // Diagonal
                dp[i][j] = match;
            } else if (delete1 >= delete2) {
                backtrace[i][j] = 2;  // Up
                dp[i][j] = delete1;
            } else {
                backtrace[i][j] = 3;  // Left
                dp[i][j] = delete2;
            }
        }
    }

    if (dp[m][n] < 0) { // if the score is negative then return the original intervals
        return make_pair(ref_intevals, read_intervals);
    }
    // Traceback to reconstruct the alignment
    int i = m;
    int j = n;
    std::string alignedSeq1;
    std::string alignedSeq2;
    alignedSeq1.reserve(m);
    alignedSeq2.reserve(n);

    // back trace the matrix to get the aligned sequences
    while ((i > 0  && j > 0)) {
        if (backtrace[i][j] == 1) {
            alignedSeq1 += seq1[i - 1];
            alignedSeq2 += seq2[j - 1];
            --i;
            --j;
        } else if (backtrace[i][j] == 2) {
            alignedSeq1 += seq1[i - 1];
            alignedSeq2 += '-';
            --i;
        } else {
            alignedSeq1 += '-';
            alignedSeq2 += seq2[j - 1];
            --j;
        }
    }
    // finish the backtrace if we reached the end of one of the sequences
    while (i > 0) {
        alignedSeq1 += seq1[i - 1];
        alignedSeq2 += '-';
        --i;
    }
    while (j > 0) {
        alignedSeq1 += '-';
        alignedSeq2 += seq2[j - 1];
        --j;
    }

    if (alignedSeq1.size() != alignedSeq2.size()) {
        spdlog::error("Error: aligned sequences are not the same size.");
        exit(1);
    }
    
    std::reverse(alignedSeq1.begin(), alignedSeq1.end());
    std::reverse(alignedSeq2.begin(), alignedSeq2.end());



    std::vector<Interval> ref_intervals_aligned;
    std::vector<Interval> read_intervals_aligned;
    // make a space before the next debug message

    // // update the intervals based on the alignment
    spdlog::debug("alignedSeq1: {}", alignedSeq1);
    spdlog::debug("alignedSeq2: {}", alignedSeq2);

    // std::cout << std::endl;
    ref_intervals_aligned = updateIntervals(alignedSeq1, seq1, ref_intevals);
    read_intervals_aligned = updateIntervals(alignedSeq2, seq2, read_intervals);

    // Visualize the corrected intervals compared to the reference

    return make_pair(ref_intervals_aligned, read_intervals_aligned);
}


bool TandemTwister::areIntervalOverlapping(Interval interval1, Interval interval2,const uint16_t motif_size){ //  unsigned motif_length
    // this function checks if the intervals are overlapping (for now it checks for exact overlap, 
    //maybe change it to check if the intervals are at most x base pair away from each other)
    // check if the interval are overlapping or at least 20% of the motif size  away from each other to the left or to the right
    if (std::max(interval1.start, interval2.start) <= std::min(interval1.end, interval2.end)) {
        return true;
    }

    if(motif_size==0){
        if (std::abs(static_cast<int>(interval1.end) - static_cast<int>(interval2.start)) <=2) {
            return true;
        }

        // Additional check for reversed order
        if (std::abs(static_cast<int>(interval2.end) - static_cast<int>(interval1.start)) <=2) {
            return true;
        }
    }
    return false;

}


void TandemTwister::print_intervals_in_cluster(std::vector<std::pair<std::string,std::vector<Interval>>>  &reads_cluster){

        for (auto const &cluster : reads_cluster)
        {
            std::cout << cluster.first << "\t";
            for (auto const &interval : cluster.second)
            { 
                std::cout << interval.start << " " << interval.end << "\t";
            }
            std::cout << std::endl;
        }
         std::cout << std::endl;
}


bool TandemTwister::IsIntervalOnLeft(Interval &interval1, Interval &interval2){
    if (interval2.start < interval1.start){
        return true;
    }
    return false;
}

// int16_t findMostOccurringCopyNumberIndex(const std::vector<std::pair<std::string, std::vector<Interval>>>& cluster) {
//     if (cluster.empty()) {
//         return -1;
//     }

//     std::unordered_map<int16_t, std::pair<uint16_t, uint16_t>> frequencyMap;
//     int16_t mostOccurringIndex = -1;
//     uint16_t maxFrequency = 0;

//     for (uint16_t i = 0; i < cluster.size(); ++i) {
//         uint16_t copyNumber = cluster[i].second.size();
//         ++frequencyMap[copyNumber].first;
//         if (frequencyMap[copyNumber].first == 1) {
//             frequencyMap[copyNumber].second = i;
//         }

//         if (frequencyMap[copyNumber].first > maxFrequency) {
//             maxFrequency = frequencyMap[copyNumber].first;
//             mostOccurringIndex = frequencyMap[copyNumber].second;
//         }
//     }

//     return mostOccurringIndex;
// }



int16_t findMostOccurringCopyNumberIndex(const std::vector<std::pair<std::string, std::vector<Interval>>>& cluster) {
    if (cluster.empty()) {
        return -1;
    }

    // if the size of the cluster is 2 then check if the CN of both is the same then return the index of the first read, if not return the index of the minimum CN
    if (cluster.size() == 2){
        if (cluster[0].second.size() == cluster[1].second.size()){
            return 0;
        }
        else{
            return cluster[0].second.size() > cluster[1].second.size() ? 0 : 1;
        }
    }
    std::unordered_map<int16_t, std::pair<uint16_t, std::vector<uint16_t>>> frequencyMap;
    std::vector<int16_t> mostFrequentCopyNumbers;
    uint16_t maxFrequency = 0;

    for (uint16_t i = 0; i < cluster.size(); ++i) {
        uint16_t copyNumber = cluster[i].second.size();
        auto& freqData = frequencyMap[copyNumber];
        ++freqData.first;
        freqData.second.push_back(i);

        if (freqData.first > maxFrequency) {
            maxFrequency = freqData.first;
            mostFrequentCopyNumbers.clear();
            mostFrequentCopyNumbers.push_back(copyNumber);
        } else if (freqData.first == maxFrequency) {
            if (std::find(mostFrequentCopyNumbers.begin(), mostFrequentCopyNumbers.end(), copyNumber) == mostFrequentCopyNumbers.end()) {
                mostFrequentCopyNumbers.push_back(copyNumber);
            }
        }
    }
    // sort the mostFrequentCopyNumbers form the biggest to the smallest
    std::sort(mostFrequentCopyNumbers.begin(), mostFrequentCopyNumbers.end(), std::greater<int16_t>());

    // if there is only one most frequent copy number then return the index of the first read with that copy number
    if (mostFrequentCopyNumbers.size() == 1) {
        return frequencyMap[mostFrequentCopyNumbers[0]].second[0];
    } else {
        // in this case we have more than one copy number that occurs the most, so return the index of the median copy number
        std::sort(mostFrequentCopyNumbers.begin(), mostFrequentCopyNumbers.end());
        return frequencyMap[mostFrequentCopyNumbers[mostFrequentCopyNumbers.size() / 2]].second[0];
    }

}


void TandemTwister::print_0_1_matrix(std::vector<std::vector<uint32_t>> clustered_intervals){
    for (auto const &cluster : clustered_intervals)
    {
        for (auto const &interval : cluster)
        {
            std::cout << interval << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}




uint32_t TandemTwister::correct_intervals(std::vector<std::pair<std::string,std::vector<Interval>>>  &reads_cluster, const float consensus_range, const uint16_t motif_size, const std::unordered_map<std::string, std::string>& sequences){
    
    std::sort(reads_cluster.begin(), reads_cluster.end(), [](auto  a, auto b) {
        return a.second.size() > b.second.size();
    });

    // if there is only one read in the cluster then return the number of intervals in the read
    if (reads_cluster.size() == 1){
        return reads_cluster[0].second.size();
    }

    // if the most occuring copy number is 0 the return 0
    if (reads_cluster[0].second.size() == 0){
        return 0;
    }


    // make the main read as the medien copy number (the main read here is the one with the most intervals)
    std::string main_read = sequences.at(reads_cluster[0].first);
    // this vector will hold the number of reads that overlap with each interval in the main read
    std::vector<uint32_t> intervalsConsesus(reads_cluster[0].second.size(), 1); 
    // this vector will hold the number of reads that is on the left of each interval in the main read
    std::vector<uint32_t> left_intervals(reads_cluster[0].second.size(),0);


    std::string read = "";

    for(uint16_t i = 1; i < reads_cluster.size(); ++i){  
        read = sequences.at(reads_cluster[i].first);
        if (read.size() == 0){

            continue;
        }

        auto [ref_intervals, read_intervals] = needlemanWunsch(main_read, read, AlignmentScores{2, -7, -7, -1}, reads_cluster[0].second, reads_cluster[i].second);        
        uint32_t reads_intervals_idx = 0;
        uint32_t ref_intervals_idx = 0;
        Interval ref_interval = Interval();

        while (ref_intervals_idx < ref_intervals.size()){
            if (reads_intervals_idx >= read_intervals.size() ){
                break;
            }
            ref_interval= ref_intervals[ref_intervals_idx];
  
        
            if (areIntervalOverlapping(ref_interval, read_intervals[reads_intervals_idx],motif_size)){
                intervalsConsesus[ref_intervals_idx] += 1;
                ++ref_intervals_idx;
                ++reads_intervals_idx;
            }
            else{
                // check if the Interval is on the left of the read_interval
                if (IsIntervalOnLeft(ref_interval, read_intervals[reads_intervals_idx])){

                    left_intervals[ref_intervals_idx] += 1;
                    ++reads_intervals_idx;
                }
                else{
                    
                    ++ref_intervals_idx;
                }

            }
            
        }
    }
    uint32_t copy_number_left = 0;
    uint32_t copy_number = 0;
    
    for (auto &interval: left_intervals){
        if (interval >= consensus_range){
            ++copy_number_left;
        }
    }

    for (uint32_t i = 0; i < intervalsConsesus.size(); ++i){
        if (intervalsConsesus[i] >= consensus_range){
            ++copy_number;
        }
    }
    uint32_t copy_number_total = copy_number+copy_number_left;

    spdlog::debug("consensus_range: {}", consensus_range);
    spdlog::debug("copy number total: {}", copy_number_total);
    spdlog::debug("copy number left: {}", copy_number_left);



    
    // remove the intervals in the main read that are not passing the threshold
    if (copy_number_total == 0){
        return 0;
    }

    spdlog::debug("removing intervals that are not passing the threshold");
    spdlog::debug("number of main read intervals before removing: {}", reads_cluster[0].second.size());
    // remove the intervals that are not passing the threshold from the main read
    reads_cluster[0].second.erase(
        std::remove_if(
            reads_cluster[0].second.begin(),
            reads_cluster[0].second.end(),
            [&](const Interval& interval) {
                size_t index = &interval - &reads_cluster[0].second[0];
                return intervalsConsesus[index] < consensus_range;
            }
        ),
        reads_cluster[0].second.end()
    );

    spdlog::debug("number of main read intervals after removing: {}", reads_cluster[0].second.size());
    return copy_number_total;

}