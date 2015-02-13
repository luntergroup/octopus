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
#include <set>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "read_reader_implementor.h"
#include "aligned_read.h"

using std::uint_fast32_t;
using std::uint_fast8_t;
using std::uint8_t;
using std::uint32_t;
using std::int32_t;

class HtslibFacade : IReadReaderImplementor
{
public:
    HtslibFacade() = delete;
    HtslibFacade(const std::string& htslib_file_path);
    ~HtslibFacade() override;
    
    HtslibFacade(const HtslibFacade&)            = delete;
    HtslibFacade& operator=(const HtslibFacade&) = delete;
    HtslibFacade(HtslibFacade&&)                 = default;
    HtslibFacade& operator=(HtslibFacade&&)      = default;
    
    std::set<AlignedRead> fetch_reads(const GenomicRegion& a_region) override;
    uint_fast32_t get_num_reference_contigs() noexcept override;
    std::unordered_set<std::string> get_reference_contig_names() override;
    uint_fast32_t get_reference_contig_size(const std::string& contig_name) override;
    
private:
    class HtslibIterator
    {
    public:
        HtslibIterator() = delete;
        HtslibIterator(HtslibFacade& hts_facade, const GenomicRegion& a_region);
        HtslibIterator(const HtslibIterator&) = delete;
        ~HtslibIterator();
        HtslibIterator& operator=(const HtslibIterator&) = delete;
        int operator++();
        AlignedRead operator*() const;
        
    private:
        HtslibFacade& hts_facade_;
        hts_itr_t* the_htslib_HtslibIterator_;
        bam1_t* b_;
        
        char get_base(uint8_t* a_htslib_sequence, uint_fast32_t index) const noexcept;
        std::string make_sequence(uint8_t* a_htslib_sequence, uint_fast32_t sequence_length) const;
        std::vector<uint_fast8_t> make_qualities(uint8_t* qualities, uint_fast32_t sequence_length) const;
        CigarString make_cigar_string(std::uint32_t* cigar_operations, std::uint32_t cigar_length) const;
        std::string get_contig_name(int32_t htslib_tid) const;
    };
    
    std::string the_filename_;
    htsFile* the_hts_file_;
    bam_hdr_t* the_header_;
    hts_idx_t* the_index_;
    std::unordered_map<int32_t, std::string> htslib_tid_to_contig_name_;
    
    std::unordered_map<int32_t, std::string> get_htslib_tid_to_contig_name_mappings() const;
    int32_t get_htslib_tid(const std::string& reference_contig_name) const;
};

inline
uint_fast32_t HtslibFacade::get_num_reference_contigs() noexcept
{
    return the_header_->n_targets;
}

inline
uint_fast32_t HtslibFacade::get_reference_contig_size(const std::string& contig_name)
{
    return the_header_->target_len[get_htslib_tid(contig_name)];
}

inline
int HtslibFacade::HtslibIterator::operator++()
{
    return sam_itr_next(hts_facade_.the_hts_file_, the_htslib_HtslibIterator_, b_)  >= 0;
}

inline
char HtslibFacade::HtslibIterator::get_base(uint8_t* a_htslib_sequence, uint_fast32_t index) const noexcept
{
    static constexpr const char* symbol_table {"=ACMGRSVTWYHKDBN"};
    return symbol_table[bam_seqi(a_htslib_sequence, index)];
}

inline
std::string HtslibFacade::HtslibIterator::make_sequence(uint8_t* a_htslib_sequence,
                                                  uint_fast32_t sequence_length) const
{
    std::string result {};
    result.reserve(sequence_length);
    for (uint_fast32_t i = 0; i < sequence_length; ++i) {
        result.push_back(get_base(a_htslib_sequence, i));
    }
    return result;
}

inline std::vector<uint_fast8_t>
HtslibFacade::HtslibIterator::make_qualities(uint8_t* qualities, uint_fast32_t sequence_length) const
{
    std::vector<uint_fast8_t> result(sequence_length);
    std::copy(qualities, qualities + sequence_length, std::begin(result));
    return result;
}

inline
CigarString HtslibFacade::HtslibIterator::make_cigar_string(uint32_t* cigar_operations,
                                          uint32_t cigar_length) const
{
    static constexpr const char* operation_table {"MIDNSHP=X"};
    std::vector<CigarString::CigarOperation> result;
    result.reserve(cigar_length);
    for (uint32_t i {0}; i < cigar_length; ++i) {
        result.emplace_back(bam_cigar_oplen(cigar_operations[i]),
                            operation_table[bam_cigar_op(cigar_operations[i])]);
    };
    return CigarString(std::move(result));
}

inline
std::string HtslibFacade::HtslibIterator::get_contig_name(int32_t htslib_tid) const
{
    return hts_facade_.htslib_tid_to_contig_name_[htslib_tid];
}

#endif /* defined(__Octopus__htslib_facade__) */
