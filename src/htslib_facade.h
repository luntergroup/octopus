//
//  htslib_facade.h
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__htslib_facade__
#define __Octopus__htslib_facade__

#include <string>
#include <vector>
#include <cstddef> // std::size_t
#include <cstdint> // std::uint_fast32_t etc
#include <algorithm> // std::find
#include <unordered_map>
#include <stdexcept> // std::runtime_error
#include <memory> // std::unique_ptr

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "read_reader_impl.h"
#include "aligned_read.h"

class AlignedRead;

using std::uint_fast32_t;
using std::uint_fast8_t;
using std::uint8_t;
using std::uint32_t;
using std::int32_t;

auto htslib_file_deleter     = [] (htsFile* the_file) { hts_close(the_file); };
auto htslib_header_deleter   = [] (bam_hdr_t* the_header) { bam_hdr_destroy(the_header); };
auto htslib_index_deleter    = [] (hts_idx_t* the_index) { hts_idx_destroy(the_index); };
auto htslib_iterator_deleter = [] (hts_itr_t* the_iterator) { sam_itr_destroy(the_iterator); };
auto htslib_bam1_deleter     = [] (bam1_t* b) { bam_destroy1(b); };

class HtslibFacade : public IReadReaderImpl
{
public:
    HtslibFacade() = delete;
    HtslibFacade(const std::string& htslib_file_path);
    ~HtslibFacade() noexcept override = default;
    
    HtslibFacade(const HtslibFacade&)            = delete;
    HtslibFacade& operator=(const HtslibFacade&) = delete;
    HtslibFacade(HtslibFacade&&)                 = default;
    HtslibFacade& operator=(HtslibFacade&&)      = default;
    
    std::vector<std::string> get_sample_ids() override;
    std::vector<std::string> get_read_groups_in_sample(const std::string& a_sample_id);
    std::size_t get_num_reads(const GenomicRegion& a_region) override;
    SampleIdToReadsMap fetch_reads(const GenomicRegion& a_region) override;
    uint_fast32_t get_num_reference_contigs() noexcept override;
    std::vector<std::string> get_reference_contig_names() override;
    uint_fast32_t get_reference_contig_size(const std::string& contig_name) override;
    std::vector<GenomicRegion> get_regions_in_file() override;
    void close() override;
    
private:
    class HtslibIterator
    {
    public:
        HtslibIterator() = delete;
        HtslibIterator(HtslibFacade& hts_facade, const GenomicRegion& a_region);
        ~HtslibIterator() noexcept = default;
        
        HtslibIterator(const HtslibIterator&) = delete;
        HtslibIterator& operator=(const HtslibIterator&) = delete;
        
        int operator++();
        std::pair<AlignedRead, std::string> operator*() const;
        
    private:
        HtslibFacade& hts_facade_;
        std::unique_ptr<hts_itr_t, decltype(htslib_iterator_deleter)> the_iterator_;
        std::unique_ptr<bam1_t, decltype(htslib_bam1_deleter)> the_bam1_;
        
        uint_fast32_t get_read_start() const noexcept;
        uint32_t get_sequence_length() const noexcept;
        char get_base(uint8_t* a_htslib_sequence, uint32_t index) const noexcept;
        std::string get_sequence() const;
        std::vector<uint_fast8_t> get_qualities() const;
        uint32_t get_cigar_length() const noexcept;
        CigarString make_cigar_string() const;
        std::string get_read_group() const;
        std::string get_contig_name(int32_t htslib_tid) const;
        std::string get_read_name() const;
        AlignedRead::SupplementaryData get_flags() const;
        AlignedRead::MatePair::SupplementaryData get_mate_flags() const;
    };
    
    // No getting around these constants. I'll put them here so they are in plain sight.
    static constexpr const char* Read_group_tag {"RG"};
    static constexpr const char* Read_group_id_tag {"ID"};
    static constexpr const char* Sample_id_tag {"SM"};
    
    using HtsTidToContigNameMap  = std::unordered_map<int32_t, std::string>;
    using ReadGroupToSampleIdMap = std::unordered_map<std::string, std::string>;
    
    std::string the_filename_;
    std::unique_ptr<htsFile, decltype(htslib_file_deleter)> the_file_;
    std::unique_ptr<bam_hdr_t, decltype(htslib_header_deleter)> the_header_;
    std::unique_ptr<hts_idx_t, decltype(htslib_index_deleter)> the_index_;
    
    HtsTidToContigNameMap contig_name_map_;
    ReadGroupToSampleIdMap sample_id_map_;
    
    std::string get_reference_contig_name(int32_t hts_tid) const;
    ReadGroupToSampleIdMap get_read_group_to_sample_id_map() const;
    HtsTidToContigNameMap get_htslib_tid_to_contig_name_map() const;
    int32_t get_htslib_tid(const std::string& reference_contig_name) const;
    // These methods have const char* parameters so they can be declared constexpr
    bool is_type(const std::string& header_line, const char* tag) const;
    bool has_tag(const std::string& header_line, const char* tag) const;
    std::string get_tag_value(const std::string& line, const char* tag) const;
    uint64_t get_num_mapped_reads(const std::string& reference_contig_name) const;
};

#endif /* defined(__Octopus__htslib_facade__) */
