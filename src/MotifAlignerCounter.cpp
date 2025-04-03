#include "include/TandemTwister.hpp"
#include <iomanip>
using IntervalScorePair = std::pair<Interval, std::tuple<int32_t, int32_t,int32_t>>;


std::tuple<Interval,bool,uint32_t> back_track(const std::vector<int32_t>& matrix_s, const std::vector<std::pair<uint32_t, uint32_t>> & matrix_p, uint32_t j,uint32_t rows, uint32_t cols, int32_t threshold) { 
    /**
    * @brief Backtracks through a scoring matrix to find an interval and determine if it meets a threshold.
    *
    * This function traverses a scoring matrix and a corresponding pair matrix to identify an interval
    * that meets a specified score threshold. It starts from the bottom-right of the matrix and moves
    * backwards until it reaches the top or the left edge of the matrix.
    *
    * @param matrix_s A vector of integers representing the scoring matrix.
    * @param matrix_p A vector of pairs representing the backtracking pointers for the scoring matrix.
    * @param j The column index to start backtracking from.
    * @param rows The number of rows in the scoring matrix.
    * @param cols The number of columns in the scoring matrix.
    * @param threshold The score threshold to determine if the interval is a match.
    * @return A pair containing the interval and a boolean indicating if the interval meets the threshold.
    *         The interval is represented by a struct with start, end, and motif_id fields.
    *         The boolean is true if the interval's score meets or exceeds the threshold, false otherwise.
    */
    auto index = [&](int32_t i, int32_t j) { return i * cols + j; };
    uint32_t intervalEnd = j;
    uint32_t i = rows - 1;
    Interval interval;
    interval.start = 0;
    interval.end = j;
    interval.motif_id = 0;
    std::pair<uint32_t, uint32_t> lastTrack;
    bool isMatch = false;
    int32_t current_score =  matrix_s[index(i,j)];
    int32_t best_score = current_score;
    while (j > 0 && i > 0) {
        if (current_score  >= threshold){
            isMatch = true;
        }
        lastTrack = matrix_p[index(i, j)];
        i = lastTrack.first;
        j = lastTrack.second;
        current_score =  matrix_s[index(i,j)];
        if (current_score > best_score){
            best_score = current_score;
        }
    }
    interval.start = lastTrack.second + 1;
    interval.end = intervalEnd;
    return std::make_tuple(interval, isMatch, best_score);
}




std::pair<std::vector<int32_t>, std::vector<std::pair<uint32_t, uint32_t>>> SWA(
    const std::string& seq1,
    const std::string& seq2,
    uint16_t match_score,
    int16_t mismatch_penalty,
    int16_t gap_penalty,
    int lcp_length,
    std::vector<int32_t> &scoring_matrix,
    std::vector<std::pair<uint32_t, uint32_t>> &matrix_pointers
) {    /**
    * @brief this function implements the Smith-Waterman algorithm for local alignment.

    * This function implements the Smith-Waterman algorithm for local alignment. It uses dynamic programming
    * to calculate the alignment score matrix and backtracking pointers for two sequences. The function also
    * keeps track of the number of occurrences of each nucleotide in the second sequence.

    *
    * @param seq1  The first sequence to align.
    * @param seq2 The second sequence to align.
    * @param match_score The score for a match.
    * @param mismatch_penalty The penalty for a mismatch.
    * @param gap_penalty The penalty for a gap.
    * @param lcp_length The length of the longest common prefix between the two sequences.
    * @param scoring_matrix A vector of the previous scoring matrix.
    * @param matrix_pointers A vector of the previous backtracking pointers.
    * @return A pair containing the new scoring matrix and the backtracking pointers.
    */

    uint32_t rows = seq1.length() + 1;
    uint32_t cols = seq2.length() + 1;
    auto index = [&](uint32_t i, uint32_t j) { return i * cols + j; };

    scoring_matrix.resize(rows * cols, 0);
    matrix_pointers.resize(rows * cols, {0, 0});

    for (uint32_t i = lcp_length + 1; i <= seq1.length(); ++i) {
        char char1 = seq1[i - 1];
        #pragma omp smid (i)
        for (uint32_t j = 1; j <= seq2.length(); ++j) {
            char char2 = seq2[j - 1];
            bool isMatch = (char1 == char2);
            int32_t match = scoring_matrix[index(i - 1, j - 1)] + (isMatch ? match_score : mismatch_penalty);
            int32_t delete_ = scoring_matrix[index(i - 1, j)] + gap_penalty;
            int32_t insert = scoring_matrix[index(i, j - 1)] + gap_penalty;
            int32_t score = std::max({match, delete_, insert, 0});

            if (score == delete_ ) {
                matrix_pointers[index(i, j)] = {i - 1, j};
            } else if (score == match) {
                matrix_pointers[index(i, j)] = {i - 1, j - 1};
            } else if (score == insert) {
                matrix_pointers[index(i, j)] = {i, j - 1};
            } else {
                matrix_pointers[index(i, j)] = {i - 1, j - 1}; // No backtracking if score is 0
            }

            scoring_matrix[index(i, j)] = score;
        }
    }


return {std::cref(scoring_matrix), std::cref(matrix_pointers)};
}


std::vector<Interval> findBestIntervalSet(std::vector<std::vector<Interval>>& allIntervals) {

    /**
    * @brief  find the best set of non-overlapping intervals.
    *
    * @param allIntervals  A vector of vectors of intervals.
    */

    std::unordered_map<uint16_t, std::vector<Interval>> intervalsByEnd;
    for (const auto& motifIntervals : allIntervals) {
        for (const auto& interval : motifIntervals) {
            intervalsByEnd[interval.end].push_back(interval);
        }
    }

    // Find the maximum end position
    uint16_t maxEnd = 0;
    for (const auto& entry : intervalsByEnd) {
        maxEnd = std::max(maxEnd, entry.first);
    }

    // DP table to store maximum scores
    std::vector<uint32_t> maxScores(maxEnd + 1, 0);

    // Vector to store selected intervals for each position
    std::vector<std::vector<Interval>> selectedIntervals(maxEnd + 1);

    // Iterate through each position
    for (uint16_t i = 1; i <= maxEnd; ++i) {
        // Initialize with the previous maximum score and intervals
        maxScores[i] = maxScores[i - 1];
        selectedIntervals[i] = selectedIntervals[i - 1];

        // Check if there are intervals ending at the current position
        if (intervalsByEnd.find(i) != intervalsByEnd.end()) {
            // Iterate through all intervals ending at position i
            for (const auto& interval : intervalsByEnd[i]) {
                // Calculate the potential new score
                uint32_t newScore = maxScores[interval.start - 1] + interval.score;

                // Update the DP table and selected intervals if the new score is higher
                if (newScore > maxScores[i]) {
                    maxScores[i] = newScore;
                    selectedIntervals[i] = selectedIntervals[interval.start - 1];
                    selectedIntervals[i].push_back(interval);
                }
            }
        }
    }


    // Return the selected intervals for the last position
    return selectedIntervals[maxEnd];
}




std::vector<Interval>  findMotifIntervals(std::string &seq, std::string& motif, double allowError, uint16_t matchScore, int16_t mismatchPenalty, int16_t gapPenalty, uint16_t motif_id,int lcp_length,std::vector<int32_t>& scoring_matrix, std::vector<std::pair<uint32_t, uint32_t>>& matrix_pointers) {
    /**
    * @brief  Find the intervals for a motif in a sequence.

    * This function finds the intervals for a motif in a sequence. It uses the Smith-Waterman algorithm
    * to calculate the scoring matrix and backtracking pointers for the sequence and motif. It then
    * backtracks through the scoring matrix to find the intervals that meet the threshold for a match.

    *
    * @param seq  The sequence to search for the motif.
    * @param motif The motif to search for in the sequence.
    * @param allowError The error rate allowed for the motif.
    * @param matchScore The score for a match.
    * @param mismatchPenalty The penalty for a mismatch.
    * @param gapPenalty The penalty for a gap.
    * @param motif_id The ID of the motif.
    * @param lcp_length The length of the longest common prefix between the sequence and the motif.
    * @param scoring_matrix The scoring matrix for the sequence and motif.
    * @param matrix_pointers The backtracking pointers for the scoring matrix.
    * @return A vector of intervals that meet the threshold for a match.
    */

   
    float motifSize = motif.size(); 
    float seq_size = seq.size();
    if (seq_size == 0 || seq_size < motifSize * allowError) {
        return {};
    }
 
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    std::transform(motif.begin(), motif.end(), motif.begin(), ::toupper);

    auto [matrix_s, matrix_p] = SWA(motif, seq, matchScore, mismatchPenalty, gapPenalty, lcp_length, scoring_matrix, matrix_pointers);

    std::vector<Interval> intervals;
    intervals.reserve(seq.size());
    uint32_t rows = motif.size()+1;
    uint32_t cols = seq.size()+1;
    std::vector<int32_t> lastRow(cols);
    std::copy(matrix_s.begin() + (rows - 1) * cols, matrix_s.end(), lastRow.begin());
    
    allowError = motifSize <= 2 ? 1 : allowError;
    
    //Interval tempInterval;
    if (motifSize == 3 && motif.find('N') != std::string::npos) {
        allowError = 0.5;
    }
    auto threshold = std::ceil(allowError * motif.size() * matchScore);
    for (uint32_t i = 0; i < lastRow.size(); i++) {
        int32_t & score = lastRow[i];
        int32_t idx = i;
        intervals.push_back({static_cast<uint32_t>(idx), static_cast<uint32_t>(idx), motif_id,0}); 
        if (idx == 0 && score !=0 ) { 
            continue;
        }
        auto [tempInterval, isMatch,best_score] = back_track(matrix_s, matrix_p, i,rows,cols, threshold);
        if (!isMatch ) {
            continue;
        }
        Interval& interval = intervals.back();
        interval = tempInterval;
        interval.motif_id = motif_id;
        interval.score = best_score;
    }

    return std::cref(intervals);
}


std::vector<Interval> TandemTwister::findBestPath(std::string &seq, std::vector<std::string> &motifs, uint16_t matchScore, int16_t mismatchPenalty, int16_t gapPenalty,const std::vector<uint16_t> &lcp ) {
    
    /**
    * @brief  Find the best path for a sequence.

    * This function finds the best path for a sequence by searching for multiple motifs in the sequence.
    * It uses the findMotifIntervals function to find the intervals for each motif and then flattens and
    * sorts the intervals and scores. It then finds the columns that meet the threshold for a match.

    *
    * @param seq  The sequence to search for the motifs.
    * @param motifs The motifs to search for in the sequence.
    * @param matchScore The score for a match.
    * @param mismatchPenalty The penalty for a mismatch.
    * @param gapPenalty The penalty for a gap.
    * @param nucleotides_occurences A map to store the number of occurrences of each nucleotide in the sequence.
    * @return A vector of intervals that meet the threshold for a match.
    */


    
    std::vector<std::vector<Interval>> allIntervals;
    allIntervals.reserve(motifs.size());
    std::vector<Interval> intervals;
    std::vector<int32_t> thresholds;
    double allowError =0.0;
    uint16_t motif_size = 0; 



    // declare the scoring matrix and matrix pointers
    std::vector<int32_t> scoring_matrix;
    std::vector<std::pair<uint32_t, uint32_t>> matrix_pointers;
    std::sort(motifs.begin(), motifs.end());


    for (uint16_t i = 0; i < motifs.size(); i++) {
        motif_size = motifs[i].size();
        // check if motif is 3 and N is part of the motif 
        if (this->min_match_ratio_s  == 0) {
            allowError = motif_size < 11 ? this->motif_match_ratio[motif_size] : this->min_match_ratio_l;
        }
        else {
            allowError = motif_size < 11 ? this->min_match_ratio_s : this->min_match_ratio_l;
        }
        if (motifs[i].find('N') != std::string::npos) {
            // reduce the allowEroor by the number of Ns in the motif
            uint16_t Ns = std::count(motifs[i].begin(), motifs[i].end(), 'N');
            double NsRatio = Ns / motif_size;
            allowError = allowError - NsRatio - 0.1;
       
        }
        int32_t threshold = std::ceil(allowError * motifs[i].size() * matchScore);
        thresholds.push_back(threshold);

        std::vector<Interval> intervals = findMotifIntervals(seq, motifs[i], allowError, matchScore, mismatchPenalty, gapPenalty, i, lcp[i], scoring_matrix, matrix_pointers);

        allIntervals.push_back(intervals); 
    }

    return findBestIntervalSet(allIntervals);
}