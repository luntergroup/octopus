//
//  htslib_sam_facade.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__htslib_sam_facade__
#define __Octopus__htslib_sam_facade__

#include <string>
#include <vector>
#include <unordered_map>
#include <cstddef>   // size_t
#include <cstdint>   // std::uint_fast32_t etc
#include <algorithm> // std::find
#include <stdexcept> // std::runtime_error
#include <memory>    // std::unique_ptr
#include <tuple>     // std::pair
#include <boost/filesystem/path.hpp>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "read_reader_impl.hpp"
#include "aligned_read.hpp"

namespace fs = boost::filesystem;

using std::uint_fast32_t;
using std::uint_fast8_t;
using std::uint8_t;
using std::uint32_t;
using std::int32_t;

auto hts_file_deleter        = [] (htsFile* file)       { hts_close(file); };
auto htslib_header_deleter   = [] (bam_hdr_t* header)   { bam_hdr_destroy(header); };
auto htslib_index_deleter    = [] (hts_idx_t* index)    { hts_idx_destroy(index); };
auto htslib_iterator_deleter = [] (hts_itr_t* iterator) { sam_itr_destroy(iterator); };
auto htslib_bam1_deleter     = [] (bam1_t* b)           { bam_destroy1(b); };

class HtslibSamFacade : public IReadReaderImpl
{
public:
    using SequenceType    = AlignedRead::SequenceType;
    using SampleIdType    = IReadReaderImpl::SampleIdType;
    using Reads           = IReadReaderImpl::Reads;
    using SampleReadMap   = IReadReaderImpl::SampleReadMap;
    using SizeType        = IReadReaderImpl::SizeType;
    using ReadGroupIdType = std::string;
    
    HtslibSamFacade() = delete;
    HtslibSamFacade(const fs::path& file_path);
    ~HtslibSamFacade() noexcept override = default;
    
    HtslibSamFacade(const HtslibSamFacade&)            = delete;
    HtslibSamFacade& operator=(const HtslibSamFacade&) = delete;
    HtslibSamFacade(HtslibSamFacade&&)                 = default;
    HtslibSamFacade& operator=(HtslibSamFacade&&)      = default;
    
    std::vector<SampleIdType> get_samples() override;
    std::vector<ReadGroupIdType> get_read_groups_in_sample(const SampleIdType& sample) override;
    size_t count_reads(const GenomicRegion& region) override;
    size_t count_reads(const SampleIdType& sample, const GenomicRegion& region) override;
    GenomicRegion find_head_region(const GenomicRegion& region, size_t target_coverage) override;
    SampleReadMap fetch_reads(const GenomicRegion& region) override;
    Reads fetch_reads(const SampleIdType& sample, const GenomicRegion& region) override;
    unsigned get_num_reference_contigs() noexcept override;
    std::vector<std::string> get_reference_contig_names() override;
    SizeType get_reference_contig_size(const std::string& contig_name) override;
    std::vector<GenomicRegion> get_possible_regions_in_file() override;
    
    void open() override;
    void close() override;
    
private:
    using HtsTidType = int32_t;
    
    class HtslibIterator
    {
    public:
        HtslibIterator() = delete;
        HtslibIterator(HtslibSamFacade& hts_facade, const GenomicRegion& region);
        ~HtslibIterator() noexcept = default;
        
        HtslibIterator(const HtslibIterator&) = delete;
        HtslibIterator& operator=(const HtslibIterator&) = delete;
        
        bool operator++();
        AlignedRead operator*() const;
        HtslibSamFacade::ReadGroupIdType get_read_group() const;
        
    private:
        HtslibSamFacade& hts_facade_;
        std::unique_ptr<hts_itr_t, decltype(htslib_iterator_deleter)> hts_iterator_;
        std::unique_ptr<bam1_t, decltype(htslib_bam1_deleter)> hts_bam1_;
    };
    
    static constexpr const char* Read_group_tag    {"RG"};
    static constexpr const char* Read_group_id_tag {"ID"};
    static constexpr const char* Sample_id_tag     {"SM"};
    
    fs::path file_path_;
    std::unique_ptr<htsFile, decltype(hts_file_deleter)> hts_file_;
    std::unique_ptr<bam_hdr_t, decltype(htslib_header_deleter)> hts_header_;
    std::unique_ptr<hts_idx_t, decltype(htslib_index_deleter)> hts_index_;
    
    std::unordered_map<std::string, HtsTidType> hts_tid_map_;
    std::unordered_map<HtsTidType, std::string> contig_name_map_;
    std::unordered_map<ReadGroupIdType, SampleIdType> sample_map_;
    
    void init_maps();
    HtsTidType get_htslib_tid(const std::string& contig_name) const;
    const std::string& get_contig_name(HtsTidType hts_tid) const;
    uint64_t get_num_mapped_reads(const std::string& contig_name) const;
};

#endif /* defined(__Octopus__htslib_sam_facade__) */
