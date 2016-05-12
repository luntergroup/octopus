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
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include <boost/filesystem/path.hpp>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "read_reader_impl.hpp"
#include "aligned_read.hpp"

class GenomicRegion;
class ContigRegion;

class HtslibSamFacade : public IReadReaderImpl
{
public:
    using Path = boost::filesystem::path;
    
    using IReadReaderImpl::SampleIdType;
    using IReadReaderImpl::ReadContainer;
    using IReadReaderImpl::SampleReadMap;
    using IReadReaderImpl::SizeType;
    using IReadReaderImpl::CoveragePair;
    
    using SequenceType = AlignedRead::SequenceType;
    
    using ReadGroupIdType = std::string;
    
    HtslibSamFacade() = delete;
    HtslibSamFacade(Path file_path);
    ~HtslibSamFacade() noexcept override = default;
    
    HtslibSamFacade(const HtslibSamFacade&)            = delete;
    HtslibSamFacade& operator=(const HtslibSamFacade&) = delete;
    HtslibSamFacade(HtslibSamFacade&&)                 = default;
    HtslibSamFacade& operator=(HtslibSamFacade&&)      = default;
    
    bool is_open() const noexcept override;
    void open() override;
    void close() override;
    
    std::vector<SampleIdType> extract_samples() override;
    std::vector<ReadGroupIdType> extract_read_groups_in_sample(const SampleIdType& sample) override;
    
    bool has_contig_reads(const GenomicRegion::ContigNameType& contig) override;
    
    std::size_t count_reads(const GenomicRegion& region) override;
    std::size_t count_reads(const SampleIdType& sample, const GenomicRegion& region) override;
    
    CoveragePair find_covered_subregion(const GenomicRegion& region, std::size_t max_coverage) override;
    CoveragePair find_covered_subregion(const SampleIdType& sample, const GenomicRegion& region,
                                        std::size_t max_coverage) override;
    CoveragePair find_covered_subregion(const std::vector<SampleIdType>& samples, const GenomicRegion& region,
                                        std::size_t max_coverage) override;
    
    SampleReadMap fetch_reads(const GenomicRegion& region) override;
    ReadContainer fetch_reads(const SampleIdType& sample, const GenomicRegion& region) override;
    SampleReadMap fetch_reads(const std::vector<SampleIdType>& samples, const GenomicRegion& region) override;
    
    unsigned count_reference_contigs() override;
    std::vector<std::string> extract_reference_contig_names() override;
    SizeType get_reference_contig_size(const std::string& contig_name) override;
    std::vector<GenomicRegion> extract_possible_regions_in_file() override;
    
private:
    using HtsTidType = std::int32_t;
    
    static constexpr std::size_t default_reserve_ {10'000'000};
    
    class HtslibIterator
    {
    public:
        HtslibIterator() = delete;
        HtslibIterator(HtslibSamFacade& hts_facade, const GenomicRegion& region);
        HtslibIterator(HtslibSamFacade& hts_facade, const GenomicRegion::ContigNameType& contig);
        ~HtslibIterator() noexcept = default;
        
        HtslibIterator(const HtslibIterator&) = delete;
        HtslibIterator& operator=(const HtslibIterator&) = delete;
        
        bool operator++();
        AlignedRead operator*() const;
        
        HtslibSamFacade::ReadGroupIdType get_read_group() const;
        
        bool is_good() const noexcept;
        std::size_t get_begin() const noexcept;
        
    private:
        struct HtsIteratorDeleter
        {
            void operator()(hts_itr_t* iterator) const { sam_itr_destroy(iterator); }
        };
        struct HtsBam1Deleter
        {
            void operator()(bam1_t* b) const { bam_destroy1(b); }
        };
        
        HtslibSamFacade& hts_facade_;
        std::unique_ptr<hts_itr_t, HtsIteratorDeleter> hts_iterator_;
        std::unique_ptr<bam1_t, HtsBam1Deleter> hts_bam1_;
    };
    
    struct HtsFileDeleter
    {
        void operator()(htsFile* file) const { hts_close(file); }
    };
    struct HtsHeaderDeleter
    {
        void operator()(bam_hdr_t* header) const { bam_hdr_destroy(header); }
    };
    struct HtsIndexDeleter
    {
        void operator()(hts_idx_t* index) const { hts_idx_destroy(index); }
    };
    
    Path file_path_;
    
    std::unique_ptr<htsFile, HtsFileDeleter> hts_file_;
    std::unique_ptr<bam_hdr_t, HtsHeaderDeleter> hts_header_;
    std::unique_ptr<hts_idx_t, HtsIndexDeleter> hts_index_;
    
    std::unordered_map<std::string, HtsTidType> hts_tids_;
    std::unordered_map<HtsTidType, std::string> contig_names_;
    std::unordered_map<ReadGroupIdType, SampleIdType> sample_names_;
    
    std::vector<SampleIdType> samples_;
    
    void init_maps();
    HtsTidType get_htslib_tid(const std::string& contig_name) const;
    const std::string& get_contig_name(HtsTidType hts_tid) const;
    uint64_t get_num_mapped_reads(const std::string& contig_name) const;
};

#endif /* defined(__Octopus__htslib_sam_facade__) */
