#include "include/TandemTwister.hpp"


std::vector<std::vector<std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>>>> TandemTwister::createIterators(unsigned int num_processes, unsigned int chunk_size) {
    /**
     * @brief Create iterators for the chunks
     * 
     * This function creates iterators for the chunks.
     * 
     * @param num_processes The number of processes
     * @param chunk_size The chunk size
     * @return The iterators
     * 
     */
    std::vector<std::vector<std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>>>> iterators;
    unsigned int total_chromosomes = this->input_chromosome_names.size();
    unsigned int push_idx = 0; // Index zum Hinzufügen der Chromosomen zu den Iteratoren

    for (unsigned int i = 0; i < num_processes; ++i) {
        std::vector<std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>>> chunk;
        for (unsigned int j = 0; j < chunk_size && push_idx < total_chromosomes; ++j) {
            std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>> chromosome_info;
            chromosome_info[this->input_chromosome_names[push_idx]] = this->TR_regions_by_chr.at(this->input_chromosome_names[push_idx]);
            chunk.push_back(chromosome_info);
            ++push_idx;
        }
        if (!chunk.empty()) {
            iterators.push_back(chunk);
        }
    }

    // Falls noch Chromosomen übrig sind, diese in den letzten Chunk hinzufügen
    if (push_idx < total_chromosomes) {
        std::vector<std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>>> chunk;
        while (push_idx < total_chromosomes) {
            std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>> chromosome_info;
            chromosome_info[this->input_chromosome_names[push_idx]] = this->TR_regions_by_chr.at(this->input_chromosome_names[push_idx]);
            chunk.push_back(chromosome_info);
            ++push_idx;
        }
        iterators.push_back(chunk);
    }

    return iterators;
}

std::unordered_map<unsigned int, unsigned int> TandemTwister::process_contig(size_t reference_start, uint32_t* cigar_data, size_t n_cigar, size_t stop_pos) {
    /**        
       * @brief Process the contig.
        * 
        * This function processes the contig and creates a map from each reference position to the corresponding read position.
        * this map is used to cut the contig in the tandem repeat regions.
        *
        * @param reference_start The reference start position.
        * @param cigar_data The cigar data.
        * @param n_cigar The number of cigar data.
        * @param stop_pos The stop position.
        * @return The reference to read mapping.
    */

    size_t current_ref_position = reference_start;
    size_t current_read_position = 0;
    std::unordered_map<unsigned int, unsigned int> ref_to_read;
    ref_to_read.reserve(n_cigar * 2);

    for (unsigned int i = 0; i < n_cigar; ++i) {
        size_t cigartype = bam_cigar_op(cigar_data[i]);
        size_t cigar_length = bam_cigar_oplen(cigar_data[i]);
        if (current_ref_position >= stop_pos) {
            break;
        }

        switch (cigartype) {
            case 0: // match
            case 7: // sequence match
            case 8: // mismatch
                for (size_t j = 0; j < cigar_length; ++j) {
                    // check if the current reference position is already in the ref_to_read dictionary if then skip the current reference position
                    if (ref_to_read.find(current_ref_position + j) != ref_to_read.end()) {
                        continue;
                    }
                    ref_to_read[current_ref_position + j] = current_read_position + j;
                }
                current_ref_position += cigar_length;
                current_read_position += cigar_length;
                break;
            case 1: // insertion
                ref_to_read[current_ref_position] = current_read_position;
                current_read_position += cigar_length;
                break;
            case 3: // skipped region
            case 4: // soft clipping
                current_read_position += cigar_length;
                break;
            case 2: // deletion
                for (size_t j = 0; j < cigar_length; ++j) {
                    if (ref_to_read.find(current_ref_position + j) != ref_to_read.end()) {
                        continue;
                    }
                    ref_to_read[current_ref_position + j] = current_read_position;
                }
                current_ref_position += cigar_length;
                break;
            case 5: // hard clipping
                // Do nothing
                break;
        }
    }

    return ref_to_read;
}

unsigned int  TandemTwister::getLastRegionPosInRead(unsigned int lastReadPos,const std::vector<std::tuple<std::string,std::vector<std::string>>>& chunk, unsigned int chunk_idx){
    /** 
    * @brief Get the last region position in the read.
    *
    * This function gets the last region position in the read based on the last read position.
    *   
    * @param lastReadPos The last read position.
    * @param chunk The chunk.
    * @param chunk_idx The chunk index.
    * @return The last region position in the read.
    */

    // get the region from the chunk which end is smaller than the lastReadPos
    // if the region is not found return 0
    std::string chr, start, end;
    unsigned int stop_pos = 0;
    for (unsigned int i = chunk_idx; i < chunk.size(); ++i){
        std::tie(chr, start, end) = parse_region(chunk,i); 
        if (static_cast<unsigned int>(std::stoi(end)) > lastReadPos){
           break;
        }
        else{
            stop_pos = std::stoi(end);
        }
    }
    return stop_pos;

}



std::vector<std::string> TandemTwister::process_chunk_assembly(const std::vector<std::tuple<std::string,std::vector<std::string>>>& chunk,std::string chromosome,  hts_idx_t* idx, bam_hdr_t* h, samFile* fp, faidx_t *fai,std::vector<vcfRecordInfoAssembly> &vcfRecords,std::string &cut_reads_fasta) {
    /**
    * @brief Process the chunk for assembly.

    * This function processes the chunk for assembly by cutting the reads in the tandem repeat regions, aligning the motifs to the contig and finding the best motif combination.
    *
    * @param chunk The chunk.
    * @param chromosome The chromosome.
    * @param idx The index.
    * @param h The header.
    * @param fp The file pointer.
    * @param fai The FASTA index.
    * @param vcfRecords The VCF records.
    * @param cut_reads_fasta The cut reads in FASTA format.
    * @return The results of the chunk.
    */

    auto reads = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>(bam_init1(), bam_destroy1);

    std::vector<std::string> chunk_results= {};
    std::string region_result;
    std::string chr, start, end;
    size_t reference_start = 0;
    uint32_t start_pos = 0;
    uint32_t end_pos = 0;
    uint32_t start_pos_in_contig = 0;
    uint32_t end_pos_in_contig = 0;
    hts_itr_t* iter = sam_itr_querys(idx, h, chromosome.c_str());
    if (!iter) {
        throw std::runtime_error("Failed to fetch iterator for region " + chromosome);
    }
    
    unsigned int regions_idx = 0;

    std::string TR_type ="";

    while (sam_itr_next (fp, iter, reads.get()) >= 0) {
        if (reads->core.flag & BAM_FUNMAP || reads->core.flag & BAM_FDUP  || reads->core.flag & BAM_FQCFAIL ) { // || reads->core.flag & BAM_FSECONDARY 
            continue;
        }
        if (regions_idx >= chunk.size()){
            break;
        }
        uint32_t* cigar_data = bam_get_cigar(reads);
        reference_start = reads->core.pos;
        uint8_t* seq_data = bam_get_seq(reads); 
        if (seq_data == nullptr) {
            std::cerr << "Error: bam_get_seq returned nullptr\n";
            exit(1);
        }
        
        std::unordered_map<unsigned int , unsigned int> ref_to_read =  process_contig(reference_start,cigar_data , reads->core.n_cigar,reads->core.pos + reads->core.l_qseq);
        unsigned int seq_length = reads->core.l_qseq; // get the sequence length
        auto qseq = std::make_unique<char[]>(seq_length+1);
        for (size_t i = 0; i < seq_length; ++i) {
            qseq.get()[i] = seq_nt16_str[bam_seqi(seq_data, i)]; // convert the sequence to string
        }
        qseq.get()[seq_length] = '\0'; // null terminate the string
        std::string ref_seq = "";
        
        std::string read_name = bam_get_qname(reads);
        while (regions_idx < chunk.size()){
            std::tie(chr, start, end) = parse_region(chunk,regions_idx);  
            start_pos = std::stoi(start) -this->padding; 
            end_pos = std::stoi(end) + this->padding; 
            if (start_pos == end_pos || end_pos < start_pos ) {
                spdlog::warn("WARNING: The start and end positions are the same or the end position is smaller than the start position in region {}:{}-{}", chr, start_pos, end_pos);
                ++regions_idx;
                continue;
            }
            // check if the start and end position are already in the ref_to_read dictionary, in this case we can cut the read from the contig(only accept if the region is completely covered by the read)
            if (ref_to_read.find(start_pos)  != ref_to_read.end() && ref_to_read.find(end_pos) != ref_to_read.end()){ 
                spdlog::info("Processing region: {}:{}-{}", chr, start_pos, end_pos);

                start_pos_in_contig = ref_to_read[start_pos];
                end_pos_in_contig = ref_to_read[end_pos];
                // check if there is an insertion at the start_pos -1 position by comparing it to start_pos
                if (ref_to_read.find(start_pos -1) != ref_to_read.end()) {
                    if ((start_pos_in_contig - ref_to_read[start_pos -1]) != 1) {
                        start_pos_in_contig = ref_to_read[start_pos -1];
                        spdlog::warn("Insertion at the start_pos -1 position found, new start position in contig: {}", start_pos_in_contig);
                    }

                }
                spdlog::debug("Start position in contig: {}", start_pos_in_contig);
                spdlog::debug("End position in contig: {}", end_pos_in_contig);
                std::string cut_read(qseq.get() + start_pos_in_contig -1, end_pos_in_contig - start_pos_in_contig +1);
                spdlog::debug("Cut read: {}", cut_read); 
                spdlog::debug("--------------------------------");
                if (this-> keep_cut_sequence){
                    // get the region, read name and the cut read sequence
                    cut_reads_fasta +=  ">" + chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos)  + "_" + read_name +"\n" +  cut_read + "\n";
                }
                std::vector<Interval> path = {}; 

                if (std::get<1>(chunk[regions_idx])[0].length() < 7)
                {
                    TR_type = "STR";
                }
                else
                {
                    TR_type = "VNTR";
                }

                std::vector<std::string> motifs_str = std::get<1>(chunk[regions_idx]);
                const std::vector<uint16_t> lcp_motifs = computeLCP(motifs_str);

             
                path = findBestPath(cut_read, motifs_str, match_score, mismatch_score, gap_score, lcp_motifs);
          
                std::string motifs_strs = "";
                for (uint16_t i = 0; i < motifs_str.size(); ++i) {
                    if (i == motifs_str.size() - 1) {
                        motifs_strs += motifs_str[i];
                    } 
                    else {
                        motifs_strs += motifs_str[i] + ",";
                    }
               
                }
                int len = 0; 
                std::string  updated_region = chr + ":" + std::to_string(start_pos) + "-" + std::to_string(end_pos);
                ref_seq = std::string(fai_fetch(fai, updated_region.c_str(), &len));
                if (ref_seq == "") {
                    std::cerr << "Failed to fetch sequence for region " << chr << std::endl;
                    ++regions_idx;
                    continue;
                }
                //std::unordered_map<char, uint16_t> nucleotides_occurences_reference;
                auto path_reference = findBestPath(ref_seq, motifs_str, match_score, mismatch_score, gap_score, lcp_motifs);
                int16_t ref_CN = path_reference.size();
                if (ref_CN == 0) {
                    spdlog::warn("Unexpected behavior: no motif was found in the reference sequence for region {}:{}-{}", chr, start_pos, end_pos);
                    ++regions_idx;
                    continue;
                }
                
                std::vector<uint16_t> motif_occurrences_allele = get_motif_ids(path);
                std::vector<uint16_t> motif_occurrences_ref = get_motif_ids(path_reference);
    
                std::vector<std::string> motifs_span_allele = {};
                std::vector<std::string> motifs_span_ref = {};

                for (const auto& interval : path) {
                    motifs_span_allele.push_back(("(" + std::to_string(interval.start) + "-" + std::to_string(interval.end) + ")"));
                }
                for (const auto& interval : path_reference) {
                    motifs_span_ref.push_back(("(" + std::to_string(interval.start) + "-" + std::to_string(interval.end) + ")"));
                }

                std::string motif_ids_str = get_motif_ids_str(motif_occurrences_allele);

        
                region_result = std::get<0>(chunk[regions_idx]) + "\t" + motifs_strs + "\t" + std::to_string(path.size()) + "\t"  + std::to_string(ref_CN) +  "\t" + motif_ids_str + "\n";;

                if (path.size() == 0) {
                    vcfRecordInfoAssembly recordInfo(std::get<0>(chunk[regions_idx]),motifs_str,motif_occurrences_allele,motif_occurrences_ref,motifs_span_allele,motifs_span_ref,TR_type,chr,start_pos,static_cast<uint16_t>(path.size()),ref_CN,0,0,cut_read,ref_seq);
                    vcfRecords.push_back(recordInfo);
                }
                else {
                    vcfRecordInfoAssembly recordInfo(std::get<0>(chunk[regions_idx]),motifs_str,motif_occurrences_allele,motif_occurrences_ref,motifs_span_allele,motifs_span_ref,TR_type,chr,start_pos,static_cast<uint16_t>(path.size()),ref_CN,path[0].start,path.back().end,cut_read,ref_seq);
                    vcfRecords.push_back(recordInfo);
                }
                
                ref_seq.clear();
                chunk_results.push_back(region_result);
                path.clear();
                region_result.clear();
                ++regions_idx;
            }
            else{

                // check if the read covers the region
                if ((reference_start > start_pos) || (reference_start > end_pos)  ){
                    ++regions_idx;
                }
                // in this case we need to move to the next contig, because the current contig does not contain anymore regions
                else{
                    ref_to_read.clear();                  
                    break;
                }
            }
        }
        ref_to_read.clear();
    }

    return chunk_results;

}

void TandemTwister::processRegionsForAssemblyInput() {
    /**
     * Process regions for assembly input
     * 
     * This function processes the regions for assembly input.
     */

    std::vector<pid_t> pids;
    unsigned int num_processes = std::min(this->num_threads, static_cast<unsigned int>(this->TR_regions_by_chr.size())); // number of processes to create is the minimum between the number of threads and the number of provided chromosomes
    unsigned int chunk_size = this->TR_regions_by_chr.size() / this->num_threads; // number of chromosomes to process in each process
    chunk_size = chunk_size == 0 ? 1 : chunk_size; // if the number of chromosomes is smaller than the number of threads then the chunk size is 1
    std::vector<std::vector<std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>>>> iterators = createIterators(num_processes, chunk_size);
    // sort the iterators by chromosome name
    std::sort(iterators.begin(), iterators.end(), [](const std::vector<std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>>> &a, const std::vector<std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>>> &b) {
        return a[0].begin()->first > b[0].begin()->first;
    });
    

    for (size_t i = 0; i < iterators.size(); ++i) {
        pid_t pid = fork();
        if (pid == 0) {
            try {
                samFile* fp = NULL;
                bam_hdr_t* h = NULL;
                hts_idx_t* idx = NULL;
                faidx_t* fai = NULL;
                open_bam_file(this->input_bam, fp, h, idx);
                open_reference_file(this->input_reference, fai);
                auto chunk = iterators[i];
                // sort the chunk by chromosome name
                std::sort(chunk.begin(), chunk.end(), [](const std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>> &a, const std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>>>> &b) {
                    return a.begin()->first > b.begin()->first;
                });
                
                for (size_t j = 0; j < iterators[i].size(); ++j) {

                    std::string cut_reads_fasta = "";
                    auto& chunk = iterators[i];
                    std::string chromosome = chunk[j].begin()->first;
                    
                    std::vector<std::tuple<std::string, std::vector<std::string>>> chunkElement_vector = std::move(chunk[j][chromosome]);
                    std::vector<std::string> results_chunk = {};
                    // Reservieren Sie Speicher für results_chunk, um häufige Speicherzuweisungen zu vermeiden
                    std::vector<vcfRecordInfoAssembly> chunk_genotype_records = {};
                    

                    results_chunk = process_chunk_assembly(chunkElement_vector, chromosome, idx, h, fp, fai, chunk_genotype_records, cut_reads_fasta);
                    std::string vcf_file_chunk = this->output_path + this->sampleName + "_" + chromosome + ".vcf";

     
                    std::string cut_reads_file_chunk = this->output_path + "cut_reads_" + this->sampleName + "_" + chromosome + ".fasta";
  
                    if (this->keep_cut_sequence){
                        std::ofstream infile_cutReads(cut_reads_file_chunk);
                        if (!infile_cutReads.is_open()) {
                            std::cerr << "Failed to open file " << cut_reads_file_chunk << std::endl;
                            return;
                        }
                        infile_cutReads << cut_reads_fasta;
                        infile_cutReads.close();
                    }
                    writeRecordsToVcfAssembly(chunk_genotype_records, vcf_file_chunk);

                    
                    chunk_genotype_records.clear(); // Stellen Sie sicher, dass chunk_genotype_records ebenfalls geleert wird
                }

                fai_destroy(fai);
                bam_hdr_destroy(h);
                hts_idx_destroy(idx);
                sam_close(fp);

                exit(0);           
            
            } catch (const std::runtime_error& e) {
                std::cerr << e.what() << std::endl;
                return;
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

    /*
    * merge the results from the different processes together
    * write the results to the output files
    * gzip the vcf file
    * index the vcf file
    * remove the intermediate files
    *      
    */

    std::string cut_reads_fasta = "";
    std::string genotype_result = "";
    // open the main vcf file
    
    // Ensure the directory exists
    std::string output_dir = this->output_file_vcf.substr(0, this->output_file_vcf.find_last_of('/'));
    if (!output_dir.empty()) {
        std::filesystem::create_directories(output_dir);
    }

    htsFile* outfile_vcf = bcf_open(this->output_file_vcf.c_str(), "wz");
    // write the this->vcf_header to the vcf file

    if (outfile_vcf == NULL) {
        std::cerr << " at the end " << std::endl;
        std::cerr << "Failed to open file " << this->output_file_vcf << std::endl;
        exit(1);
    }
    if (bcf_hdr_write(outfile_vcf, this->vcf_header) != 0) {
        fprintf(stderr, "Error writing VCF header.\n");
        return;
    }

    
    bcf_hdr_t *hdr = nullptr;
    std::vector<bcf1_t*> all_records;
    std::vector<std::tuple<std::string, int, int, bcf1_t*>> records_with_pos; // chr, pos, rid, record      
    // open the output files of the different processes and merge the results
    for (size_t i = 0; i < this->input_chromosome_names.size(); ++i) {
        std::string chromosome = this->input_chromosome_names[i];
        std::string cut_reads_file_chunk = this->output_path + "cut_reads_" + this->sampleName + "_" + chromosome + ".fasta";
        std::ifstream infile_cutReads(cut_reads_file_chunk);
        if (this->keep_cut_sequence){
            if (!infile_cutReads.is_open()) {
                std::cerr << "Failed to open file " << cut_reads_file_chunk << " results of process " << i << " will be skipped" << std::endl;
                exit(1);
            }
        }
 

        
        std::string vcf_file_chunk = this->output_path + this->sampleName + "_" + chromosome + ".vcf";

        htsFile* infile_vcf = hts_open(vcf_file_chunk.c_str(), "r");
       
        if (infile_vcf == NULL) {
            std::cerr << "Failed to open file " << vcf_file_chunk << std::endl;
            exit(1);
        }
        
        hdr = bcf_hdr_read(infile_vcf);
        
        bcf1_t *rec = bcf_init();
        while (bcf_read(infile_vcf, hdr, rec) == 0) {
            bcf_unpack(rec, BCF_UN_ALL);
            
            // Create a copy of the record
            bcf1_t *rec_copy = bcf_dup(rec);
            bcf_unpack(rec_copy, BCF_UN_ALL);
            
            std::string chrom = bcf_hdr_id2name(hdr, rec_copy->rid);
            int rid = rec_copy->rid; // Reference ID for proper chromosome ordering
            
            records_with_pos.emplace_back(chrom, rec_copy->pos, rid, rec_copy);
        }
        

        bcf_destroy(rec);
        hts_close(infile_vcf);
        std::stringstream buffer;
        buffer << infile_cutReads.rdbuf();
        cut_reads_fasta += buffer.str();
        std::stringstream buffer4;
        infile_cutReads.close();
        std::remove(cut_reads_file_chunk.c_str());
        std::remove(vcf_file_chunk.c_str());
    }
    // Sort all records by chromosome and position
    std::sort(records_with_pos.begin(), records_with_pos.end(),
    [](const auto& a, const auto& b) {
        // First compare by reference ID (rid) - this ensures correct chromosome order
        if (std::get<2>(a) != std::get<2>(b)) {
            return std::get<2>(a) < std::get<2>(b);
        }
        // Then compare by position
        return std::get<1>(a) < std::get<1>(b);
    });

    
    
    // Write sorted records to final VCF
    for (auto& record_data : records_with_pos) {
        bcf1_t* rec = std::get<3>(record_data);
        if (bcf_write(outfile_vcf, this->vcf_header, rec) < 0) {
            std::cerr << "Error writing VCF record to file." << std::endl;
        }
        bcf_destroy(rec);
    }
    bcf_hdr_destroy(hdr);
    // gzip the vcf file
    if (this->keep_cut_sequence){
        this->outfile_cutReads << cut_reads_fasta;
        outfile_cutReads.flush();
        outfile_cutReads.close();
    }
    this->outfile << genotype_result;
    this->outfile.flush();

    this->outfile.close();

    bcf_hdr_destroy(this->vcf_header);
    bcf_close(outfile_vcf);
    
    // index the vcf file
    std::string vcf_file_index = this->output_file_vcf + ".tbi";
    if (bcf_index_build(this->output_file_vcf.c_str(), 0) != 0) {
        std::cerr << "Failed to index the vcf file " << this->output_file_vcf << std::endl;  
    }
    this -> outfile.close();
}

