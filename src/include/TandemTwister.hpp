/*


Copyright <2024> <Lion Ward Al Raei>

this code was made available only for reviewing purposes. it is not allowed to be used or distributed without the permission of the authors.

*/


#ifndef TandemTwister_hpp
#define TandemTwister_hpp

#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <fstream>
#include <filesystem>
#include <sys/wait.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <iostream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <chrono>
#include <unistd.h>
#include <mlpack/core.hpp>
#include <iomanip>
#include <unordered_set>
#include <map>
#include <ranges>
#include <numeric>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
using namespace std;

struct Interval {
    uint32_t start;
    uint32_t end;
    uint32_t motif_id;
    uint16_t score;
    //float purity;

    bool operator<(const Interval& rhs) const {
        return start < rhs.start || (start == rhs.start && end < rhs.end);
    }
    bool operator==(const Interval& rhs) const {
        return start == rhs.start && end == rhs.end && motif_id == rhs.motif_id;
    }
    bool operator!=(const Interval& rhs) const {
        return start != rhs.start || end != rhs.end;
    }
    bool operator>(const Interval& rhs) const {
        return start > rhs.start || (start == rhs.start && end > rhs.end);
    }
    bool operator<=(const Interval& rhs) const {
        return start <= rhs.start || (start == rhs.start && end <= rhs.end ) ;
    }
};

struct regionResult {
    std::vector<std::vector<std::tuple<std::string, std::vector<Interval>, std::string>>> clusters = {};
    std::vector<std::pair<std::string, std::vector<Interval>>>  noise_cluster = {};
    std::vector<std::string> alt_seqs = {};
    std::vector<uint16_t> indices_consensus_read = {};
    std::vector<uint32_t> copy_numbers = {};
    std::vector<uint32_t> cn_occurences_in_cluster = {};
    uint16_t num_reads_covering_TR = 0;
    regionResult(
        std::vector<std::vector<std::tuple<std::string, std::vector<Interval>, std::string>>> clusters,
        std::vector<std::pair<std::string, std::vector<Interval>>> noise_cluster,
        std::vector<std::string> alt_seqs,
        std::vector<uint16_t> indices_consensus_read,
        std::vector<uint32_t> copy_numbers,
        std::vector<uint32_t> cn_occurences_in_cluster,
        uint16_t num_reads_covering_TR
    ) : clusters(clusters), noise_cluster(noise_cluster), alt_seqs(alt_seqs), indices_consensus_read(indices_consensus_read),
        copy_numbers(copy_numbers), cn_occurences_in_cluster(cn_occurences_in_cluster), num_reads_covering_TR(num_reads_covering_TR) {}
};

struct ReadInfo
{
    std::string readSeq;
    std::string readName;
    std::string readState;
    unsigned int readPos_to_ref;
    bool operator==(const ReadInfo &other) const
    {
        return readSeq == other.readSeq && readName == other.readName && readState == other.readState;
    }
    
};



struct left_interval{
    Interval interval;
    uint16_t read_idx; // the index of the read in the cluster
    uint16_t idx; // the index of the position in the left_intervals
    std::vector<uint16_t> positions_for_0; // positions in the cluster that should be 0 
};


struct AlignmentScores {
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
};

struct vcfRecordInfoReads {
    std::string region;
    std::vector<std::string> motifs;
    std::vector<::vector<uint16_t>> motif_occurrences_alleles;
    std::vector<uint16_t> motif_occurrences_REF;
    std::vector<std::vector<std::string>> motifs_intervals_alleles;
    std::vector<std::string> motifs_intervals_REF;
    std::vector<uint32_t> alleles_CN;
    std::vector<std::string> alleles_seq;
    std::string reference_seq;
    std::string TR_type;
    std::string chr;
    uint32_t start_pos;
    uint32_t ref_CN;
    uint32_t numReadsInRegion;
    std::vector<std::uint16_t> cluster_sizes;
    vcfRecordInfoReads(const std::string &region,
                const std::vector<std::string> &motifs,
                const std::vector<std::vector<uint16_t>> &motif_occurrences_alleles,
                const std::vector<uint16_t> &motif_occurrences_REF,
                const std::vector<std::vector<std::string>> &motifs_intervals_alleles,
                const std::vector<std::string>& motifs_intervals_REF,
                const std::vector<uint32_t>& alleles_CN,
                const std::vector<std::string>& alleles_seq,
                const std::string& reference_seq,
                const std::string& TR_type,
                const std::string& chr,
                const uint32_t &start_pos,
                const uint32_t &ref_CN,
                const uint32_t &numReadsInRegion,
                const std::vector<std::uint16_t> &cluster_sizes
               )
        : region(region),
        motifs(motifs),
        motif_occurrences_alleles(motif_occurrences_alleles),
        motif_occurrences_REF(motif_occurrences_REF),
        motifs_intervals_alleles(motifs_intervals_alleles),
        motifs_intervals_REF(motifs_intervals_REF),
        alleles_CN(alleles_CN),
        alleles_seq(alleles_seq),
        reference_seq(reference_seq),
        TR_type(TR_type),
        chr(chr),
        start_pos(start_pos),
        ref_CN(ref_CN),
        numReadsInRegion(numReadsInRegion),
        cluster_sizes(cluster_sizes) {}
};





struct vcfRecordInfoAssembly {
    std::string region;
    std::vector<std::string> motifs;
    std::vector<uint16_t> motif_occurrences_allele;
    std::vector<uint16_t> motif_occurrences_ref;
    std::vector<std::string> motifs_intervals_allele;
    std::vector<std::string> motifs_intervals_ref;
    std::string TR_type;
    std::string chr;
    uint32_t start_pos;
    uint16_t allele_CN;
    uint16_t ref_CN;
    uint16_t alleleSpanStart;
    uint16_t alleleSpanEnd ;
    std::string allele_seq;
    std::string ref_seq;
    vcfRecordInfoAssembly(const std::string& region,
                const std::vector<std::string>& motifs,
                const std::vector<uint16_t>& motif_occurrences_allele,
                const std::vector<uint16_t>& motif_occurrences_ref,
                const std::vector<std::string>& motifs_span_allele,
                const std::vector<std::string>& motifs_span_ref,
                const std::string& TR_type,
                const std::string& chr,
                uint32_t start_pos,
                uint16_t allele_CN,
                uint16_t ref_CN,
                uint16_t alleleSpanStart,
                uint16_t alleleSpanEnd,
                std::string allele_seq,
                std::string ref_seq)

        : region(region),
        motifs(motifs),
        motif_occurrences_allele(motif_occurrences_allele),
        motif_occurrences_ref(motif_occurrences_ref),
        motifs_intervals_allele(motifs_span_allele),
        motifs_intervals_ref(motifs_span_ref),
        TR_type(TR_type),
        chr(chr),
        start_pos(start_pos),
        allele_CN(allele_CN),
        ref_CN(ref_CN),
        alleleSpanStart(alleleSpanStart),
        alleleSpanEnd(alleleSpanEnd),
        allele_seq(allele_seq),
        ref_seq(ref_seq) {}
};


class TandemTwister {
public:
    //input arguments Files
    std::string region_file="";
    std::string input_bam="";
    std::string output_path="";
    std::string output_file="";
    std::string input_reference="";
    std::string input_reference_fai="";
    std::string output_file_vcf="";
    std::string sampleName ="";
    bcf_hdr_t* vcf_header;
    htsFile *vcf_file;
    //reads extraction parameters
    uint16_t padding = 0;
    //genotyping parameters
    double min_match_ratio_l=0.5;
    double min_match_ratio_s=0;
    uint16_t match_score=0;
    int16_t mismatch_score=0;
    int16_t gap_score=0;
    //clustering parameters
    float start_eps_str=0;
    float start_eps_vntr=0;
    float minPts_frac=0;
    float minPts_str =0;
    float minPts_vntr =0;
    float noise_limit_str = 0;
    float noise_limit_vntr= 0;
    uint16_t cluster_iter=0;
    //consensus parameters
    float cons_ratio_str=0;
    float cons_ratio_vntr=0;
    //other parameters
    uint32_t num_threads=1;
    uint8_t sample_sex=0;
    const uint16_t STR_max_length = 7; // maximum length of a str motif
    std::string reads_type = "CCS"; // type of reads CLR , CCS or ONTa
    bool correctCalling = false; // if 1 then correct the reads copy number otherwise use the most occuring copy number
    uint16_t  tanCon = 2; // the number of consecutive intervals to be considered a tandem run "2*motif_size"
    std::ofstream outfile;
    std::ofstream outfile_phasing;
    std::ofstream outfile_cutReads;
    std::unordered_map<std::string, std::string> args;
    bool keep_phasing_results = false;
    bool keep_cut_sequence = false;
    bool remove_outliers_zscore = false;
    bool refineTrRegions = false;
    bool bamIsTagged = true;
    bool manual_reclustering = true;
    bool somatic = true;
    float minReadsInCluster = 0.3;
    uint16_t motif_threashold_size = 0;
    uint16_t minReadsInRegion = 2;
    uint16_t quality_score = 10;
    std::string analysis_type = "";
    std::string version = "v2.0.1";
    std::string commandline ="";
    uint16_t verbose = 0;
    std::vector<float> motif_match_ratio = {0,1.,1.,1,0.55,0.5,0.65,0.55,0.625,0.65,0.5};
    std::map<std::string, std::map<unsigned int , unsigned int>> lookup_table_ONT;
    std::map<std::string, char*> lookup_table_ONT_sequences;
    std::vector<std::string> accepted_chromosome_names = {"chr1" , "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                                                    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
                                                    "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                                                    "chr20", "chr21", "chr22", "chrX", "chrY"};

    std::vector<std::string> input_chromosome_names = {};
    std::unordered_map<std::string, std::string> reference_sequences;
    std::vector<std::tuple<std::string,std::vector<std::string>>> TR_regions;
    std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>> TR_regions_by_chr;
    std::unordered_map<std::string, uint16_t> contig_lines;
    // constructor
    TandemTwister(int argc, char *argv[]);
    // print usage
    void printUsage(std::unordered_map<std::string, std::string>& args);
    ~TandemTwister();
    void open_bam_file(const std::string& bam_file, samFile* &fp, bam_hdr_t* &h, hts_idx_t* &idx);
    void open_reference_file(const std::string& reference_file, faidx_t* &fai);
    void processRegionsForLongReadsInput();
    void processRegionsForAssemblyInput();
    template<typename T>
    T getArg(const std::unordered_map<std::string, std::string>& args, const std::string& shortArg, const std::string& longArg, T defaultValue);
    std::string generateRegionGenotypeResult(
    const std::vector<std::string>& motifs,
    const std::vector<std::vector<uint16_t>>& motif_ids_alleles,
    const std::string &region,
    std::vector<uint32_t> copy_numbers,
    uint16_t CN_reference_all_runs,
    uint16_t num_reads_covering_TR
    );
    private:
        samFile* fp = NULL;
        bam_hdr_t* h = NULL;
        hts_idx_t* idx = NULL;
        faidx_t* fai = NULL;
        // command line arguments functions
        void through_error(std::unordered_map<std::string, std::string>& args, const std::string& error_message);
        bool isChrXorY(const std::string& chr);
        bool openFile(const std::string& filename);
        std::tuple<std::string, std::string, std::string> parse_region(const std::vector<std::tuple<std::string,std::vector<std::string>>>& chunk ,size_t j);
        
        
        std::unordered_map<std::string, std::string> parseCommandLine(int argc, char** argv);

        // helper functions
        std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>> regions_grouped_by_chr(std::string& regionFile,std::vector<std::string>& accepted_chromosome_names);
        std::vector<std::tuple<std::string, std::vector<std::string>>> regions_not_grouped(std::string& regionFile);
        void print_dashes_representation(const std::string& read_name, const std::string& cut_read, const std::vector<Interval>& path);
        std::vector<uint16_t> get_motif_ids(const std::vector<Interval>& path);
        std::string get_motif_ids_str(const std::vector<uint16_t>& motif_ids);
        // reference genotyper functions
        std::vector<Interval> genotype_reference(std::string &regionWithPadding, std::vector<std::string> & motif,  std::string &  seq, int len, faidx_t* fai, const std::vector<uint16_t> &lcp);
        void createVcfHeader();
        void writeRecordsToVcf(std::vector<vcfRecordInfoReads> & recordsInfo,std::string  vcf_file);
        void writeRecordsToVcfAssembly(std::vector<vcfRecordInfoAssembly> & recordsInfo,std::string  vcf_file);
        std::vector<uint16_t> computeLCP(const std::vector<std::string>& motifs);
        // assembly input functions 
        std::unordered_map< unsigned int , unsigned int> process_contig(size_t reference_start, uint32_t* cigar_data, size_t n_cigar, size_t stop_pos);
        std::vector<std::string> process_chunk_assembly(const std::vector<std::tuple<std::string,std::vector<std::string>>>& chunk,std::string chromosome,
                                                        hts_idx_t* idx, bam_hdr_t* h, samFile* fp, faidx_t *fai,std::vector<vcfRecordInfoAssembly> &vcfRecordInfoAssembly,std::string &cut_reads_fasta);
        std::vector<std::vector<std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>>>> createIterators(unsigned int num_processes, unsigned int chunk_size);
        // ONT reads functions
        unsigned int getLastRegionPosInRead(unsigned int lastReadPos,const std::vector<std::tuple<std::string,std::vector<std::string>>>& chunk,unsigned int chunk_idx);

        // // Long reads  functions
        regionResult processLongReads(samFile* &fp, hts_itr_t* &iter,const size_t start_pos, const size_t end_pos, const std::string& region, 
                                    std::vector<std::string>& motifs,uint16_t motif_size, std::string &chr,
                                    std::vector<std::pair<uint32_t, uint32_t>> &tandem_runs,
                                    std::string &region_phasing_result, std::string &cut_reads_fasta, const std::vector<uint16_t>& lcp);
        void get_hp_tag(uint8_t *&hp_tag, bam1_t *reads, uint32_t &hp_value);

        //void process_chunk_ONT(const std::vector<std::tuple<std::string,std::string>>& chunk,  hts_idx_t* idx, bam_hdr_t* h, samFile* fp, faidx_t *fai);
        std::tuple<uint32_t, uint32_t, uint16_t> findMostOccurringCopyNumber(const std::vector<std::tuple<std::string, std::vector<Interval>, std::string>>& cluster);
        std::tuple<uint32_t, uint32_t, uint16_t> findMedianCopyNumber(const std::vector<std::pair<std::string, std::vector<Interval>>>& cluster);

        std::string generate_phasing_results(std::vector<std::tuple<std::string,std::vector<Interval>, std::string>> & first_cluster, std::vector<std::tuple<std::string,std::vector<Interval>, std::string>> & second_cluster,
                                            std::vector<std::pair<std::string,std::vector<Interval>>> &noise_cluster , unsigned int CN_first, unsigned int CN_second, const std::string& region);
        std::vector<std::pair<uint32_t, uint32_t>> cut_read_in_TR_region(const  uint32_t reference_start,  uint32_t* cigar_data, size_t cigar_length,
                                                                        const std::vector<std::pair<uint32_t,uint32_t> > &tandem_runs);
        //void update_feature_vector(std::vector<vecf>& region_features, char* &cut_read, const std::vector<Interval>& path);
        void update_feature_vector(arma::mat& region_features, char* &cut_read, const std::vector<Interval>& path);
        std::tuple<std::vector<std::string>, std::vector<vcfRecordInfoReads> > process_chunk_reads(size_t start_dict, size_t end_dict, samFile* fp, bam_hdr_t* h, 
                                                                                                                           hts_idx_t* idx, faidx_t* fai,std::string & cut_reads_fasta);
        
        void cluster_by_features(const std::string& chr, 
                                        std::vector<std::tuple<std::string, std::vector<Interval>, uint8_t>>& reads_intervals, 
                                        std::vector<std::vector<std::tuple<std::string, std::vector<Interval>, std::string>>>  & clusters, 
                                        std::vector<std::pair<std::string, std::vector<Interval>>>& noise_cluster, 
                                        arma::mat &region_features, 
                                        uint16_t motif_size,
                                        bool &haplody);
        // Motif genotyper functions
        std::vector<Interval> findBestPath(std::string &seq, std::vector<std::string> &motifs, 
                                            uint16_t matchScore, int16_t mismatchPenalty, int16_t gapPenalty, const std::vector<uint16_t> &lcp);
     
        // Reads Intervals Corrrection functions 
        uint32_t correct_intervals(std::vector<std::tuple<std::string,std::vector<Interval>,std::string>>  &reads_cluster, const float consensus_range, const uint16_t motif_size, const std::unordered_map<std::string,std::string> &reads_sequences); 
        bool IsIntervalOnLeft(Interval &interval1, Interval &interval2);

        void print_intervals_in_cluster(std::vector<std::pair<std::string,std::vector<Interval>>>  &reads_cluster);
        bool areIntervalOverlapping(Interval interval1, Interval interval2,const uint16_t motif_size);
   
        std::pair<std::vector<Interval> ,std::vector<Interval>>  needlemanWunsch(const std::string& seq1, const std::string& seq2, const AlignmentScores& scores,
                                                                                const std::vector<Interval> &ref_intervals, const std::vector<Interval> &read_intervals);
        // std::pair<std::vector<Interval>, std::vector<Interval>> smithWaterman(const std::string& seq1, const std::string& seq2, const AlignmentScores& scores, const std::vector<Interval>& ref_intervals, const std::vector<Interval>& read_intervals);

        void print_0_1_matrix(std::vector<std::vector<uint32_t>> clustered_intervals);
        std::vector<Interval> findLongestTandemRun(const std::vector<Interval>& intervals, const uint16_t motif_size);
        std::vector<std::vector<Interval>> findAllTandemRuns(const std::vector<Interval>& intervals, const std::vector<std::string>& motifs, const std::string& sequence);
        bool IsMotifRun(const std::vector<Interval> & intervals, const std::vector<uint16_t>& motif_lengths); //, const unsigned int motif_size
        
        // phasing functions (consider putting them in a different class)
        std::tuple<std::vector<std::vector<uint16_t>>, std::vector<uint16_t>> runDBSCAN(const arma::mat &features_matrix_stdandarized,float eps, uint16_t minPts);
        void cluster_reads(arma::mat &region_features,std::vector<std::tuple<std::string,std::vector<Interval>,uint8_t>> &reads_intervals,
                                                            std::vector<std::vector<std::tuple<std::string,std::vector<Interval>,std::string>>> &final_clusters ,
                                                            std::vector<std::pair<std::string,std::vector<Interval>>> & noise_cluster,
                                                            uint16_t motif_length);
        arma::mat standardizeMatrix(arma::mat mmatrix);
        std::vector<uint16_t> removeOutlierReadsZScore(const arma::mat &data, float zThreshold, uint16_t columnIdx);
        float calculateFeatureMeanInCluster(const std::vector<uint16_t>& cluster, const arma::mat &region_features, uint16_t dataIndex);


};

#endif /* TandemTwister_hpp */