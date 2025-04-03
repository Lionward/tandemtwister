
#include <iostream>
#include <sstream>
#include <fstream>
#include "include/TandemTwister.hpp"

void TandemTwister::createVcfHeader(){ 

    this->vcf_header = bcf_hdr_init("w"); 
    // Add contig information to the header
    std::ifstream reference_file(this->input_reference_fai);
    if (!reference_file.is_open()) {
        std::cerr << "Error opening reference file: " << this->input_reference_fai << std::endl;
        return;
    }
    std::string line;
    uint16_t contigIdx= 0;
    bcf_hdr_append(this->vcf_header, "##fileformat=VCFv4.2");
    bcf_hdr_append(this->vcf_header, ("##source=tandemTwister" + std::string(this->version)).c_str());
    // cmmand line
    std::string command_line = "##command=" + this->commandline;
    bcf_hdr_append(this->vcf_header, command_line.c_str());
    while (std::getline(reference_file, line)) {
        std::istringstream iss(line);
        std::string contig, length;
        iss >> contig >> length;
        this->contig_lines[contig] = contigIdx;
        ++contigIdx;
        // Add contig information to the header
        std::string contig_line = "##contig=<ID=" + contig + ",length=" + length + ">";
        bcf_hdr_append(this->vcf_header, contig_line.c_str());
    }
    // Add columns for ALT and REF alleles
   

    // Add INFO fields to the header
    bcf_hdr_append(this->vcf_header, "##INFO=<ID=TR_type,Number=1,Type=String,Description=\"TR type STR/VNTR\">");

    bcf_hdr_append(this->vcf_header, "##INFO=<ID=MOTIFS,Number=1,Type=String,Description=\"Tandem repeat motif\">");
    bcf_hdr_append(this->vcf_header, "##INFO=<ID=UNIT_LENGTH_AVG,Number=1,Type=String,Description=\"Average Length of the repeat unit\">");
    
    if (this->analysis_type == "somatic" || this->analysis_type == "germline"){
        if (this->somatic == false ){
            bcf_hdr_append(this->vcf_header, "##INFO=<ID=CN_H1,Number=1,Type=String,Description=\"Number of repeats for each haplotpye\">");
            bcf_hdr_append(this->vcf_header, "##INFO=<ID=CN_H2,Number=1,Type=String,Description=\"Number of repeats for each haplotpye\">");
            bcf_hdr_append(this->vcf_header, "##INFO=<ID=MOTIF_IDs_H1,Number=1,Type=String,Description=\"Motif ids for each haplotype\">");
            bcf_hdr_append(this->vcf_header, "##INFO=<ID=MOTIF_IDs_H2,Number=1,Type=String,Description=\"Motif ids for each haplotype\">");
        }
        else {
            bcf_hdr_append(this->vcf_header, "##INFO=<ID=CN_alleles,Number=1,Type=String,Description=\"list of copy number of repeats for each haplotpye\">");
            bcf_hdr_append(this->vcf_header, "##INFO=<ID=MOTIF_IDs_alleles,Number=1,Type=String,Description=\"Motif ids for each haplotype\">");
        }
        bcf_hdr_append(this->vcf_header, "##INFO=<ID=CN_ref,Number=1,Type=String,Description=\"Number of repeats for each haplotpye\">");
        bcf_hdr_append(this->vcf_header, "##INFO=<ID=MOTIF_IDs_REF,Number=1,Type=String,Description=\"Motif ids for each haplotype\">");

        bcf_hdr_append(this->vcf_header, "##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Number of Reads supporting each allele\">");
        bcf_hdr_append(this->vcf_header, "##FORMAT=<ID=SP,Number=1,Type=String,Description=\"Span of the TR on each allel\">");

    }
    else{
        bcf_hdr_append(this->vcf_header, "##INFO=<ID=CN_hap,Number=1,Type=Integer,Description=\"Number of repeats for the given allele\">");
        bcf_hdr_append(this->vcf_header, "##INFO=<ID=CN_ref,Number=1,Type=Integer,Description=\"Number of repeats for the reference allele\">");
        bcf_hdr_append(this->vcf_header, "##INFO=<ID=MOTIF_IDs_REF,Number=1,Type=String,Description=\"Motif ids for each haplotype\">");
        bcf_hdr_append(this->vcf_header, "##INFO=<ID=MOTIF_IDs_H,Number=1,Type=String,Description=\"Motif ids for each haplotype\">");
        bcf_hdr_append(this->vcf_header, "##FORMAT=<ID=SP,Number=2,Type=String,Description=\"Span of the TR on each allel\">");

    }
    // Add FORMAT fields to the header
    bcf_hdr_append(this->vcf_header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    int32_t pass_filter_id = bcf_hdr_id2int(this->vcf_header, BCF_DT_ID, "PASS");
    int32_t lowqual_filter_id = bcf_hdr_id2int(this->vcf_header, BCF_DT_ID, "LOWQUAL");

    if (pass_filter_id < 0) {
        bcf_hdr_append(this->vcf_header, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
    }

    if (lowqual_filter_id < 0) {
        bcf_hdr_append(this->vcf_header, "##FILTER=<ID=LOWQUAL,Description=\"Low quality\">");
    }

    bcf_hdr_add_sample(this->vcf_header, this->sampleName.c_str());
    // Write the header to the VCF file

}



std::vector<int16_t> genotypeAllele_reads(std::string ref_seq, std::vector<std::string> alleles_seq){

    /*
    * This function is used to genotype the alleles of a tandem repeat region based on the sequence of both alleles
    * @param ref_seq: The reference sequence
    * @param alleles_seq: The sequence of the alleles
    * @return: a vector of genotypes
    */

    std::vector<int16_t> genotypes;
    unordered_map<std::string, uint16_t> not_ref_alleles = {};

    for (size_t i = 0; i < alleles_seq.size(); ++i) {
 
        if (alleles_seq[i] == ref_seq) {
            genotypes.push_back(0);
        }
        else {
            if (not_ref_alleles.find(alleles_seq[i]) != not_ref_alleles.end()) {
                genotypes.push_back(not_ref_alleles[alleles_seq[i]]);
            }
            else {
                // count the number of 0s until the current allele
                int16_t num_zeros = 0;
                for (size_t j = 0; j < i; ++j) {
                    if (genotypes[j] == 0) {
                        ++num_zeros;
                    }
                }
                genotypes.push_back(i +1 - num_zeros);
                not_ref_alleles[alleles_seq[i]] = i+1 - num_zeros;
            }
         }
    }
    return genotypes;
}


std::string formatStringIntervals(const std::vector<std::string>& intervals) {

    /*
    * This function is used to format the motif intervals as a string
    * @param intervals: The intervals to format
    * @return: The formatted string
    *   
    */

    std::ostringstream oss;
    for (size_t i = 0; i < intervals.size(); ++i) {
        oss <<  intervals[i];
        if (i < intervals.size() - 1) {
            oss << "_";
        }
    }
    return oss.str();
}

std::string get_motifs_list(const std::vector<std::string> &motifs){
    std::string motifs_str = std::accumulate(motifs.begin(), motifs.end(), std::string(),
        [](const std::string& a, const std::string& b) {
            return a.empty() ? b : a + "," + b;
        });

    return motifs_str;
}



void TandemTwister::writeRecordsToVcf(std::vector<vcfRecordInfoReads> & recordsInfo,std::string vcf_file){
    htsFile* vcf_temp = hts_open(vcf_file.c_str(), "w");
    if (vcf_temp == NULL){
        std::cerr << "Error opening VCF file" << std::endl;
    }
    // add the header to the file
    if (bcf_hdr_write(vcf_temp, this->vcf_header) != 0) {
        fprintf(stderr, "Error writing VCF header.\n");
        return;
    }

    for (auto const &recordInfo: recordsInfo){
     try{
        bcf1_t *vcf_record = bcf_init();   

        bcf_update_id(this->vcf_header, vcf_record, recordInfo.region.c_str());

        bool homozygous = false;
        std::vector<char*> alleles;
        std::vector<std::string> allele_strings;

        for (size_t i = 0; i < recordInfo.motifs_intervals_alleles.size(); ++i) {
            std::string allele_str = get_motifs_list(recordInfo.motifs_intervals_alleles[i]);
            allele_strings.push_back(allele_str);
            alleles.push_back(strdup(allele_strings.back().c_str()));
        }

        // sort the alleles based on their similarity to the reference sequence
        std::vector<size_t> indices(alleles.size());
        std::vector<size_t> move_to_beginning;
        std::vector<char*> sorted_alleles(alleles.size());
        std::vector<int> sorted_alleles_CN(recordInfo.alleles_CN.size());
        std::vector<std::vector<uint16_t>> sorted_motif_occurrences_alleles(recordInfo.motif_occurrences_alleles.size());
        std::vector<std::vector<std::string>> sorted_motifs_intervals_alleles(recordInfo.motifs_intervals_alleles.size());
        std::vector<std::string> sorted_sequences(recordInfo.alleles_seq.size());

        for (size_t i = 0; i < recordInfo.alleles_seq.size(); ++i) {
            if (recordInfo.alleles_seq[i] == recordInfo.reference_seq) {
                move_to_beginning.push_back(i);
            }
        }
        for (size_t i = 0; i < recordInfo.alleles_seq.size(); ++i) {
            indices[i] = i;
        }

        for (size_t i = 0; i < move_to_beginning.size(); ++i) {
            std::swap(indices[move_to_beginning[i]], indices[i] );
        }
        
        for (size_t i = 0; i < indices.size(); ++i) {
            sorted_alleles[i] = alleles[indices[i]];
            sorted_alleles_CN[i] = recordInfo.alleles_CN[indices[i]];
            sorted_motif_occurrences_alleles[i] = recordInfo.motif_occurrences_alleles[indices[i]];
            sorted_motifs_intervals_alleles[i] = recordInfo.motifs_intervals_alleles[indices[i]];
            sorted_sequences[i] = recordInfo.alleles_seq[indices[i]];
        }

        // 
  
        alleles.insert(alleles.begin(), strdup(recordInfo.reference_seq.c_str()));
        std::vector<const char*> const_alleles(alleles.begin(), alleles.end());
        bcf_update_alleles(this->vcf_header, vcf_record, const_alleles.data(), const_alleles.size());
        vcf_record->pos = recordInfo.start_pos;
        vcf_record->qual = std::min(std::round((static_cast<float>(recordInfo.numReadsInRegion) / 20.0f)* 10.0f)/10.0f, 1.0f); // change this to the expected number of reads based on the coverage
        vcf_record->rid = this->contig_lines[recordInfo.chr];
     
        int32_t sum = 0;
        for (auto const &motif: recordInfo.motifs){
            sum += motif.length();
        }
        int32_t unit_length = sum / recordInfo.motifs.size();

        // tansofrm the motif list to a string separated by commas
        std::string motifs = get_motifs_list(recordInfo.motifs);
        bcf_update_info_string(this->vcf_header, vcf_record, "TR_type", recordInfo.TR_type.c_str());
        bcf_update_info_string(this->vcf_header, vcf_record, "MOTIFS", motifs.c_str());
        bcf_update_info_string(this->vcf_header, vcf_record, "UNIT_LENGTH_AVG", std::to_string(unit_length).c_str());

        // update the alleles
        vcf_record->n_allele = sorted_sequences.size() + 1;

        //update the reference allele
        vcf_record->d.allele[0] = strdup(recordInfo.reference_seq.c_str());


        std::string motif_ids_ref = get_motif_ids_str(recordInfo.motif_occurrences_REF);
        if (motif_ids_ref.empty()){
            motif_ids_ref = "";
        }
        
        bcf_update_info_string(this->vcf_header,vcf_record, "MOTIF_IDs_REF", motif_ids_ref.c_str());
        std::string CN_ref = std::to_string(recordInfo.ref_CN);
        if (this->somatic){
            std::string motif_ids_allesles = "";

            for (size_t i = 0; i < sorted_motif_occurrences_alleles.size(); ++i) {
                std::string motif_ids = get_motif_ids_str(sorted_motif_occurrences_alleles[i]);
                if (motif_ids.empty()){
                    motif_ids = ".";
                }
                if (i == sorted_motif_occurrences_alleles.size() - 1){
                    motif_ids_allesles += motif_ids;
                }
                else{
                    motif_ids_allesles += motif_ids + "|";
                }
            }
            if (motif_ids_allesles.empty()){
                motif_ids_allesles = ".";
            }
            bcf_update_info_string(this->vcf_header,vcf_record, "MOTIF_IDs_alleles", motif_ids_allesles.c_str());
            bcf_update_info_string(this->vcf_header, vcf_record, "CN_ref", CN_ref.c_str());
            std::string CNs = "";
            for (size_t i = 0; i < sorted_alleles_CN.size(); ++i) {
                std::string CN = std::to_string(sorted_alleles_CN[i]);
                if (i == sorted_alleles_CN.size() - 1){
                    CNs += CN;
                }
                else{
                    CNs += CN + "|";
                }
            }
            if (CNs.empty()){
                CNs = ".";
            }
            bcf_update_info_string(this->vcf_header,vcf_record, "CN_alleles", CNs.c_str());
   
        }
        else{

            std::string CN_H1 = std::to_string(sorted_alleles_CN[0]);
      
            std::string CN_H2 = sorted_alleles_CN.size() == 2 ? std::to_string(sorted_alleles_CN[1]) : std::to_string(sorted_alleles_CN[0]);
       
            bcf_update_info_string(this->vcf_header, vcf_record, "CN_ref", CN_ref.c_str());
            bcf_update_info_string(this->vcf_header, vcf_record, "CN_H1", CN_H1.c_str());
            bcf_update_info_string(this->vcf_header, vcf_record, "CN_H2", CN_H2.c_str());  
      
            std::string motif_ids_h1 = get_motif_ids_str(sorted_motif_occurrences_alleles[0]);
            std::string motif_ids_h2 = "";
            if (sorted_motif_occurrences_alleles.size() >  1){
                motif_ids_h2 = get_motif_ids_str(sorted_motif_occurrences_alleles[1]);
            }
      
            if (motif_ids_h1.empty()){
                motif_ids_h1 = "";
            }
            if (motif_ids_h2.empty()){
                motif_ids_h2 = "";
            }
            bcf_update_info_string(this->vcf_header,vcf_record, "MOTIF_IDs_H1", motif_ids_h1.c_str()); 
            
            if (sorted_motifs_intervals_alleles[0] == sorted_motifs_intervals_alleles[1]){
                homozygous = true;
                bcf_update_info_string(this->vcf_header,vcf_record, "MOTIF_IDs_H2", "");
            }
            else{
                bcf_update_info_string(this->vcf_header,vcf_record, "MOTIF_IDs_H2", motif_ids_h2.c_str()); 
            }

        
        }
        for (size_t i = 0; i < sorted_sequences.size(); ++i) {
            if (sorted_sequences[i].empty()){
                sorted_sequences[i] = recordInfo.reference_seq[0];
            }
        }
        for (size_t i = 0; i < sorted_sequences.size(); ++i) {
            vcf_record->d.allele[i + 1] = strdup(sorted_sequences[i].c_str());
        }

        std::vector<int16_t> genotype = genotypeAllele_reads(recordInfo.reference_seq, sorted_sequences);
        std::vector<int> tmpia;
   
        // Resize tmpia to accommodate the number of alleles
        tmpia.resize(genotype.size());
        for (size_t i = 0; i < genotype.size(); ++i) {
            if (genotype[i] == -1) {
                tmpia[i] = bcf_gt_missing;
            } else {
                tmpia[i] = bcf_gt_unphased(genotype[i]);
            }
        }

        
        bcf_update_genotypes(this->vcf_header, vcf_record, tmpia.data(), tmpia.size());

        std::string cluster_sizes = "";
        if (homozygous) {
            tmpia[0] = recordInfo.cluster_sizes[0];
            bcf_update_format_int32(this->vcf_header, vcf_record, "DP", tmpia.data(), 1);
        } else {
            for (size_t i = 0; i < recordInfo.cluster_sizes.size(); ++i) {
            tmpia[i] = recordInfo.cluster_sizes[i];
            }

            // update the cluster sizes which is the number of reads supporting each allele
            bcf_update_format_int32(this->vcf_header, vcf_record, "DP", tmpia.data(), recordInfo.cluster_sizes.size());
        }
    
        // update the span of the TR on each allele
        std::string span = formatStringIntervals(recordInfo.motifs_intervals_REF);
        for (size_t i = 0; i < sorted_motifs_intervals_alleles.size(); ++i) {
            if (homozygous && i == 1) {
                break;
            }
            span += "," + formatStringIntervals(sorted_motifs_intervals_alleles[i]);
        }
        if (span.empty()){
            span = "";
        }

    
        const char* spans[1] = { span.c_str() };

        int ret = bcf_update_format_string(this->vcf_header, vcf_record, "SP", spans, bcf_hdr_nsamples(this->vcf_header));
        if (ret != 0){
            std::cerr << "Error updating SP field" << std::endl;
        }

        // toatal number of reads in the clusters sum tmpia
        auto nReadsClusters = std::accumulate(recordInfo.cluster_sizes.begin(), recordInfo.cluster_sizes.end(), 0);
        
        int32_t tmpi = 0; 
        for (size_t i = 0; i < recordInfo.cluster_sizes.size(); ++i) {
            tmpi += recordInfo.cluster_sizes[i];
        }

        if (nReadsClusters >= 8) { // consider changing this to a parameter
            tmpi = bcf_hdr_id2int(this->vcf_header, BCF_DT_ID, "PASS");
        } else {
            tmpi = bcf_hdr_id2int(this->vcf_header, BCF_DT_ID, "LOWQUAL");
        }
        bcf_update_filter(this->vcf_header, vcf_record, &tmpi, 1);

        // make sure the record has all the required information
        if (bcf_write1(vcf_temp, this->vcf_header, vcf_record) != 0) {
            fprintf(stderr, "Error writing VCF record temp.\n");
            bcf_destroy(vcf_record);
            hts_close(vcf_temp);
            continue;
        }
 
        bcf_destroy(vcf_record);
        } catch (const std::exception& e) {
            std::cout << "region: " << recordInfo.region << std::endl;
            std::cerr << "Error writing VCF record: " << e.what() << std::endl;
        }
    }
    hts_close(vcf_temp);
}


uint16_t genotype_allele_assembly(const std::vector<uint16_t> &ref_ids, const std::vector<uint16_t> &allele_ids){
    // in the case of the assembly theres is only one haplotype at a time so we can just compare the motif ids of the reference and the allel

    bool ref_is_allele = false;
    bool ref_not_allele = false;

    if (ref_ids.size() == allele_ids.size()){
        ref_is_allele = std::equal(ref_ids.begin(), ref_ids.end(), allele_ids.begin());
    }
    if (!ref_is_allele){
        ref_not_allele = true;
    }
    if (ref_is_allele){
        return 0;
    }
    else if (ref_not_allele){
        return 1;
    }
    else{
        return 2;
    }
}



void TandemTwister::writeRecordsToVcfAssembly(std::vector<vcfRecordInfoAssembly> & recordsInfo,std::string vcf_file){
    htsFile* vcf_temp = hts_open(vcf_file.c_str(), "w");
    if (vcf_temp == NULL){
        std::cerr << "Error opening VCF file" << std::endl;
    }
    if (bcf_hdr_write(vcf_temp, this->vcf_header) != 0) {
        fprintf(stderr, "Error writing VCF header.\n");
        return;
    }
    for (auto const &recordInfo: recordsInfo){

        bcf1_t *vcf_record = bcf_init();   
        bcf_update_id(this->vcf_header, vcf_record, recordInfo.region.c_str());
        const char * allele[2] = {recordInfo.ref_seq.c_str(),recordInfo.allele_seq.c_str() };
        bcf_update_alleles(this->vcf_header,vcf_record, allele, 2);
        vcf_record->pos = recordInfo.start_pos;
        vcf_record->rid = this->contig_lines[recordInfo.chr];
        // update the contig information
        bcf_update_info_string(this->vcf_header, vcf_record, "TR_type", recordInfo.TR_type.c_str());
        std::string motifs = get_motifs_list(recordInfo.motifs);

        // add the motif ids to the info field motif_occurrences_hap
        std::string motif_ids_H = get_motif_ids_str(recordInfo.motif_occurrences_allele);
        std::string motif_ids_REF = get_motif_ids_str(recordInfo.motif_occurrences_ref);
        bcf_update_info_string(this->vcf_header, vcf_record, "MOTIF_IDs_H", motif_ids_H.c_str());
        bcf_update_info_string(this->vcf_header, vcf_record, "MOTIF_IDs_REF", motif_ids_REF.c_str());

        
        bcf_update_info_string(this->vcf_header, vcf_record, "MOTIFS", motifs.c_str());
        // calculate the average length of the motifs in the list
        int32_t unit_length = std::accumulate(recordInfo.motifs.begin(), recordInfo.motifs.end(), 0, [](int32_t sum, const std::string& motif) {
            return sum + motif.length();
        }) / recordInfo.motifs.size();

        bcf_update_info_int32(this->vcf_header, vcf_record, "UNIT_LENGTH", &unit_length, 1);
        int32_t repeats_value[2] = {recordInfo.allele_CN, recordInfo.ref_CN}; 
        bcf_update_info_int32(this->vcf_header, vcf_record, "CN_hap", &repeats_value[0], 1);
        bcf_update_info_int32(this->vcf_header, vcf_record, "CN_ref", &repeats_value[1], 1);
        
        uint32_t *tmpia = (uint32_t*)malloc(bcf_hdr_nsamples(this->vcf_header)*1*sizeof(uint32_t));
        // Set genotypes for each sample
        uint16_t genotype = genotype_allele_assembly(recordInfo.motif_occurrences_ref, recordInfo.motif_occurrences_allele);
        switch (genotype){
            case 0:
                tmpia[0] = bcf_gt_unphased(0);
                break;
            case 1:
                tmpia[0] = bcf_gt_unphased(1);
                break;
            default:
                std::cerr << "Error in genotyping" << std::endl;

        }

        bcf_update_genotypes(this->vcf_header, vcf_record, tmpia, bcf_hdr_nsamples(this->vcf_header) *1);
        std::string spanREF  = formatStringIntervals(recordInfo.motifs_intervals_ref);
        std::string spanH = formatStringIntervals(recordInfo.motifs_intervals_allele);
        if (spanREF.empty()){
            spanREF = ".";
        }
        if (spanH.empty()){
            spanH = ".";
        }
        std::string span =  spanREF + "," + spanH;
        const char* spans[1] = { span.c_str() };
        bcf_update_format_string(this->vcf_header,vcf_record,"SP", spans,bcf_hdr_nsamples(this->vcf_header)*1 );
        

        if (bcf_write1(vcf_temp, this->vcf_header, vcf_record) != 0) {
            fprintf(stderr, "Error writing VCF record.\n");
            return;  // You might choose a different error code
        }
        bcf_destroy(vcf_record);
    }

    if (hts_close(vcf_temp) != 0) {
        fprintf(stderr, "Error closing VCF file.\n");
    }

}
