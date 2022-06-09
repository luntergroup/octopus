// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_reader_impl_hpp
#define read_reader_impl_hpp

#include <string>
#include <vector>
#include <cstddef>
#include <unordered_map>
#include <utility>
#include <functional>

#include <boost/optional.hpp>

#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"

namespace octopus { namespace io {

class IReadReaderImpl
{
public:
    using SampleName    = std::string;
    using ReadContainer = std::vector<AlignedRead>;
    using SampleReadMap = std::unordered_map<SampleName, ReadContainer>;
    using PositionList  = std::vector<GenomicRegion::Position>;
    using AlignedReadReadVisitor = std::function<bool(const SampleName&, AlignedRead)>;
    using ContigRegionVisitor = std::function<bool(const SampleName&, ContigRegion)>;
    
    virtual ~IReadReaderImpl() noexcept = default;
    
    virtual bool is_open() const noexcept = 0;
    virtual void open() = 0;
    virtual void close() = 0;
    
    virtual std::vector<SampleName> extract_samples() const = 0;
    virtual std::vector<std::string> extract_read_groups(const SampleName& sample) const = 0;
    
    virtual bool iterate(const GenomicRegion& region,
                         AlignedReadReadVisitor visitor) const = 0;
    virtual bool iterate(const SampleName& sample,
                         const GenomicRegion& region,
                         AlignedReadReadVisitor visitor) const = 0;
    virtual bool iterate(const std::vector<SampleName>& samples,
                         const GenomicRegion& region,
                         AlignedReadReadVisitor visitor) const = 0;
    
    virtual bool iterate(const GenomicRegion& region,
                         ContigRegionVisitor visitor) const = 0;
    virtual bool iterate(const SampleName& sample,
                         const GenomicRegion& region,
                         ContigRegionVisitor visitor) const = 0;
    virtual bool iterate(const std::vector<SampleName>& samples,
                         const GenomicRegion& region,
                         ContigRegionVisitor visitor) const = 0;
    
    virtual bool has_reads(const GenomicRegion& region) const = 0;
    virtual bool has_reads(const SampleName& sample,
                           const GenomicRegion& region) const = 0;
    virtual bool has_reads(const std::vector<SampleName>& samples,
                           const GenomicRegion& region) const = 0;
    
    virtual std::size_t count_reads(const GenomicRegion& region) const = 0;
    virtual std::size_t count_reads(const SampleName& sample,
                                    const GenomicRegion& region) const = 0;
    virtual std::size_t count_reads(const std::vector<SampleName>& sample,
                                    const GenomicRegion& region) const = 0;
    
    virtual PositionList extract_read_positions(const GenomicRegion& region,
                                                std::size_t max_reads) const = 0;
    virtual PositionList extract_read_positions(const SampleName& sample,
                                                const GenomicRegion& region,
                                                std::size_t max_reads) const = 0;
    virtual PositionList extract_read_positions(const std::vector<SampleName>& samples,
                                                const GenomicRegion& region,
                                                std::size_t max_reads) const = 0;
    
    virtual SampleReadMap fetch_reads(const GenomicRegion& region) const = 0;
    virtual ReadContainer fetch_reads(const SampleName& sample,
                                      const GenomicRegion& region) const = 0;
    virtual SampleReadMap fetch_reads(const std::vector<SampleName>& samples,
                                      const GenomicRegion& region) const = 0;
    virtual SampleReadMap fetch_reads(const std::vector<GenomicRegion>& regions) const = 0;
    virtual ReadContainer fetch_reads(const SampleName& sample,
                                      const std::vector<GenomicRegion>& regions) const = 0;
    virtual SampleReadMap fetch_reads(const std::vector<SampleName>& samples,
                                      const std::vector<GenomicRegion>& regions) const = 0;
    
    virtual std::vector<GenomicRegion::ContigName> reference_contigs() const = 0;
    virtual GenomicRegion::Size reference_size(const GenomicRegion::ContigName& contig) const = 0;
    
    virtual boost::optional<std::vector<GenomicRegion::ContigName>> mapped_contigs() const { return boost::none; };
    virtual boost::optional<std::vector<GenomicRegion>> mapped_regions() const { return boost::none; };
};

} // namespace io
} // namespace octopus

#endif
