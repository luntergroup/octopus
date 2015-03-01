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
#include <cstdint>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>
#include <memory>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "read_reader_impl.h"
#include "aligned_read.h"

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

inline uint_fast32_t HtslibFacade::get_num_reference_contigs() noexcept
{
    return the_header_->n_targets;
}

inline uint_fast32_t HtslibFacade::get_reference_contig_size(const std::string& contig_name)
{
    return the_header_->target_len[get_htslib_tid(contig_name)];
}

inline int HtslibFacade::HtslibIterator::operator++()
{
    return sam_itr_next(hts_facade_.the_file_.get(), the_iterator_.get(), the_bam1_.get()) >= 0;
}

inline uint_fast32_t HtslibFacade::HtslibIterator::get_read_start() const noexcept
{
    return static_cast<uint_fast32_t>(the_bam1_->core.pos);
}

inline uint32_t HtslibFacade::HtslibIterator::get_sequence_length() const noexcept
{
    return the_bam1_->core.l_qseq;
}

inline char HtslibFacade::HtslibIterator::get_base(uint8_t* a_htslib_sequence,
                                                   uint32_t index) const noexcept
{
    static constexpr const char* symbol_table {"=ACMGRSVTWYHKDBN"};
    return symbol_table[bam_seqi(a_htslib_sequence, index)];
}

inline std::string HtslibFacade::HtslibIterator::get_sequence() const
{
    std::string result {};
    auto length = get_sequence_length();
    result.reserve(length);
    for (uint32_t i = 0; i < length; ++i) {
        result.push_back(get_base(bam_get_seq(the_bam1_), i));
    }
    return result;
}

inline std::vector<uint_fast8_t> HtslibFacade::HtslibIterator::get_qualities() const
{
    auto qualities = bam_get_qual(the_bam1_);
    auto length = get_sequence_length();
    std::vector<uint_fast8_t> result;
    result.reserve(length);
    result.insert(result.begin(), qualities, qualities + length);
    return result;
}

inline uint32_t HtslibFacade::HtslibIterator::get_cigar_length() const noexcept
{
    return the_bam1_->core.n_cigar;
}

inline CigarString HtslibFacade::HtslibIterator::make_cigar_string() const
{
    static constexpr const char* operation_table {"MIDNSHP=X"};
    auto cigar_operations = bam_get_cigar(the_bam1_);
    auto length = get_cigar_length();
    std::vector<CigarOperation> result;
    result.reserve(length);
    for (uint32_t i {0}; i < length; ++i) {
        result.emplace_back(bam_cigar_oplen(cigar_operations[i]),
                            operation_table[bam_cigar_op(cigar_operations[i])]);
    };
    return CigarString(std::move(result));
}

inline std::string HtslibFacade::HtslibIterator::get_read_group() const
{
    return std::string(bam_aux2Z(bam_aux_get(the_bam1_.get(), Read_group_tag)));
}

inline std::string HtslibFacade::HtslibIterator::get_contig_name(int32_t htslib_tid) const
{
    return hts_facade_.contig_name_map_[htslib_tid];
}

inline std::string HtslibFacade::HtslibIterator::get_read_name() const
{
    return std::string {bam_get_qname(the_bam1_)};
}

inline std::string HtslibFacade::get_reference_contig_name(int32_t hts_tid) const
{
    return std::string(the_header_->target_name[hts_tid]);
}

inline bool HtslibFacade::is_type(const std::string& header_line, const char* tag) const
{
    return header_line.compare(1, 2, tag) == 0;
}

inline bool HtslibFacade::has_tag(const std::string& header_line, const char* tag) const
{
    return header_line.find(tag) != std::string::npos;
}

inline std::string HtslibFacade::get_tag_value(const std::string& line, const char* tag) const
{
    // format TAG:VALUE\t
    auto tag_position = line.find(tag);
    if (tag_position != std::string::npos) {
        auto value_position = line.find(':', tag_position) + 1;
        auto tag_value_size = line.find('\t', value_position) - value_position;
        return line.substr(value_position, tag_value_size);
    }
    throw std::runtime_error {"no " + std::string {tag} + " tag"};
}

inline uint64_t HtslibFacade::get_num_mapped_reads(const std::string& reference_contig_name) const
{
    uint64_t num_mapped, num_unmapped;
    hts_idx_get_stat(the_index_.get(), get_htslib_tid(reference_contig_name), &num_mapped, &num_unmapped);
    return num_mapped;
}

#endif /* defined(__Octopus__htslib_facade__) */
