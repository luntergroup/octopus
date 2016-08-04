// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef Octopus_read_reader_impl_hpp
#define Octopus_read_reader_impl_hpp

#include <string>
#include <vector>
#include <cstddef>
#include <unordered_map>
#include <utility>

#include <basics/genomic_region.hpp>
#include <basics/aligned_read.hpp>

namespace octopus { namespace io {

class IReadReaderImpl
{
public:
    using SampleName      = std::string;
    using ReadContainer   = std::vector<AlignedRead>;
    using SampleReadMap   = std::unordered_map<SampleName, ReadContainer>;
    using CoveragePair    = std::pair<GenomicRegion, std::vector<unsigned>>;
    
    virtual ~IReadReaderImpl() noexcept = default;
    
    virtual bool is_open() const noexcept = 0;
    virtual void open() = 0;
    virtual void close() = 0;
    
    virtual std::vector<SampleName> extract_samples() const = 0;
    virtual std::vector<std::string> extract_read_groups_in_sample(const SampleName& sample) const = 0;
    
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
    
    virtual CoveragePair find_covered_subregion(const GenomicRegion& region,
                                                std::size_t max_coverage) const = 0;
    virtual CoveragePair find_covered_subregion(const SampleName& sample,
                                                const GenomicRegion& region,
                                                std::size_t max_coverage) const = 0;
    virtual CoveragePair find_covered_subregion(const std::vector<SampleName>& samples,
                                                const GenomicRegion& region,
                                                std::size_t max_coverage) const = 0;
    
    virtual SampleReadMap fetch_reads(const GenomicRegion& region) const = 0;
    virtual ReadContainer fetch_reads(const SampleName& sample,
                                      const GenomicRegion& region) const = 0;
    virtual SampleReadMap fetch_reads(const std::vector<SampleName>& samples,
                                      const GenomicRegion& region) const = 0;
    
    virtual unsigned count_reference_contigs() const = 0;
    virtual std::vector<std::string> extract_reference_contig_names() const = 0;
    virtual ContigRegion::Size get_reference_contig_size(const std::string& contig_name) const = 0;
    virtual std::vector<GenomicRegion> extract_possible_regions_in_file() const = 0;
};

} // namespace io
} // namespace octopus

#endif
