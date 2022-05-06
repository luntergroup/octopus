// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef htslib_sam_facade_hpp
#define htslib_sam_facade_hpp

#include <string>
#include <vector>
#include <unordered_map>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "basics/aligned_read.hpp"
#include "read_reader_impl.hpp"

namespace octopus {

class GenomicRegion;
class ContigRegion;

namespace io {

class HtslibSamFacade : public IReadReaderImpl
{
public:
    using Path = boost::filesystem::path;
    
    using IReadReaderImpl::SampleName;
    using IReadReaderImpl::ReadContainer;
    using IReadReaderImpl::SampleReadMap;
    using IReadReaderImpl::PositionList;
    using IReadReaderImpl::AlignedReadReadVisitor;
    using IReadReaderImpl::ContigRegionVisitor;
    
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    
    using ReadGroupIdType = std::string;
    
    HtslibSamFacade() = delete;
    
    HtslibSamFacade(Path file_path);
    HtslibSamFacade(Path sam_out, Path sam_template);
    
    HtslibSamFacade(const HtslibSamFacade&)            = delete;
    HtslibSamFacade& operator=(const HtslibSamFacade&) = delete;
    HtslibSamFacade(HtslibSamFacade&&)                 = default;
    HtslibSamFacade& operator=(HtslibSamFacade&&)      = default;
    
    ~HtslibSamFacade() override;
    
    bool is_open() const noexcept override;
    void open() override;
    void close() override;
    
    std::vector<SampleName> extract_samples() const override;
    std::vector<ReadGroupIdType> extract_read_groups(const SampleName& sample) const override;
    
    bool iterate(const GenomicRegion& region,
                 AlignedReadReadVisitor visitor) const override;
    bool iterate(const SampleName& sample,
                 const GenomicRegion& region,
                 AlignedReadReadVisitor visitor) const override;
    bool iterate(const std::vector<SampleName>& samples,
                 const GenomicRegion& region,
                 AlignedReadReadVisitor visitor) const override;
    
    bool iterate(const GenomicRegion& region,
                 ContigRegionVisitor visitor) const override;
    bool iterate(const SampleName& sample,
                 const GenomicRegion& region,
                 ContigRegionVisitor visitor) const override;
    bool iterate(const std::vector<SampleName>& samples,
                 const GenomicRegion& region,
                 ContigRegionVisitor visitor) const override;
    
    bool has_reads(const GenomicRegion& region) const override;
    bool has_reads(const SampleName& sample,
                   const GenomicRegion& region) const override;
    bool has_reads(const std::vector<SampleName>& samples,
                   const GenomicRegion& region) const override;
    
    std::size_t count_reads(const GenomicRegion& region) const override;
    std::size_t count_reads(const SampleName& sample,
                            const GenomicRegion& region) const override;
    std::size_t count_reads(const std::vector<SampleName>& samples,
                            const GenomicRegion& region) const override;
    
    PositionList extract_read_positions(const GenomicRegion& region,
                                        std::size_t max_reads) const override;
    PositionList extract_read_positions(const SampleName& sample,
                                        const GenomicRegion& region,
                                        std::size_t max_reads) const override;
    PositionList extract_read_positions(const std::vector<SampleName>& samples,
                                        const GenomicRegion& region,
                                        std::size_t max_reads) const override;
    
    SampleReadMap fetch_reads(const GenomicRegion& region) const override;
    ReadContainer fetch_reads(const SampleName& sample,
                              const GenomicRegion& region) const override;
    SampleReadMap fetch_reads(const std::vector<SampleName>& samples,
                              const GenomicRegion& region) const override;
    
    GenomicRegion::Size reference_size(const GenomicRegion::ContigName& contig) const override;
    std::vector<GenomicRegion::ContigName> reference_contigs() const override;
    boost::optional<std::vector<GenomicRegion::ContigName>> mapped_contigs() const override;
    
    void write(const AlignedRead& read);
    
private:
    using HtsTid = std::int32_t;
    
    static constexpr std::size_t defaultReserve_ {1'000'000};
    
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
    struct HtsBam1Deleter
    {
        void operator()(bam1_t* b) const { bam_destroy1(b); }
    };
    
    class HtslibIterator
    {
    public:
        HtslibIterator() = delete;
        
        HtslibIterator(const HtslibSamFacade& hts_facade, const GenomicRegion& region);
        HtslibIterator(const HtslibSamFacade& hts_facade, const GenomicRegion::ContigName& contig);
        
        HtslibIterator(const HtslibIterator&) = delete;
        HtslibIterator& operator=(const HtslibIterator&) = delete;
        
        ~HtslibIterator() noexcept = default;
        
        bool operator++();
        AlignedRead operator*() const;
        
        HtslibSamFacade::ReadGroupIdType read_group() const;
        
        bool is_good() const noexcept;
        ContigRegion region() const;
        std::size_t begin() const noexcept;
    
    private:
        struct HtsIteratorDeleter
        {
            void operator()(hts_itr_t* iterator) const { sam_itr_destroy(iterator); }
        };
        struct HtsBam1Deleter
        {
            void operator()(bam1_t* b) const { bam_destroy1(b); }
        };
        
        const HtslibSamFacade& hts_facade_;
        
        std::unique_ptr<hts_itr_t, HtsIteratorDeleter> hts_iterator_;
        std::unique_ptr<bam1_t, HtsBam1Deleter> hts_bam1_;
    };
    
    Path file_path_;
    
    std::unique_ptr<htsFile, HtsFileDeleter> hts_file_;
    std::unique_ptr<bam_hdr_t, HtsHeaderDeleter> hts_header_;
    std::unique_ptr<hts_idx_t, HtsIndexDeleter> hts_index_;
    
    std::unordered_map<GenomicRegion::ContigName, HtsTid> hts_targets_;
    std::unordered_map<HtsTid, GenomicRegion::ContigName> contig_names_;
    std::unordered_map<ReadGroupIdType, SampleName> sample_names_;
    
    std::vector<SampleName> samples_;
    
    void init_maps();
    HtsTid get_htslib_target(const GenomicRegion::ContigName& contig) const;
    const GenomicRegion::ContigName& get_contig_name(HtsTid target) const;
    std::uint64_t get_num_mapped_reads(const GenomicRegion::ContigName& contig) const;
    ReadContainer fetch_all_reads(const GenomicRegion& region) const;
    void set_fixed_length_data(const AlignedRead& read, bam1_t* result) const;
    void write(const AlignedRead& read, bam1_t* result) const;
};

} // namespace io
} // namespace octopus

#endif
