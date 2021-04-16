// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef downsampler_hpp
#define downsampler_hpp

#include <cstddef>
#include <unordered_map>

#include "config/common.hpp"
#include "concepts/mappable.hpp"
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "containers/mappable_flat_set.hpp"
#include "containers/mappable_flat_multi_set.hpp"
#include "containers/mappable_map.hpp"

namespace octopus { namespace readpipe {

/**
 Downsampler removes reads from a collection in order to reduce the coverage at positions
 which exceed a maximum bound.
 */
class Downsampler
{
public:
    struct Report
    {
        class DownsampleRegion : public Mappable<DownsampleRegion>
        {
            GenomicRegion region_;
            std::size_t num_reads_;
        public:
            const GenomicRegion& mapped_region() const noexcept { return region_; }
            std::size_t num_reads() const noexcept { return num_reads_; }
            DownsampleRegion(GenomicRegion region, std::size_t num_reads);
        };
        using DownsampleRegionSet = MappableFlatSet<DownsampleRegion>;
        DownsampleRegionSet downsampled_regions;
    };
    
    Downsampler() = default;
    
    Downsampler(unsigned trigger_coverage, unsigned target_coverage);
    
    Downsampler(const Downsampler&)            = default;
    Downsampler& operator=(const Downsampler&) = default;
    Downsampler(Downsampler&&)                 = default;
    Downsampler& operator=(Downsampler&&)      = default;
    
    ~Downsampler() = default;
    
    // Returns the number of reads removed
    Report downsample(ReadContainer& reads) const;
    
private:
    unsigned trigger_coverage_ = 10'000;
    unsigned target_coverage_  = 10'000;
};

using DownsamplerReportMap = std::unordered_map<SampleName, Downsampler::Report>;

template <typename Map>
auto downsample(Map& reads, const Downsampler& downsampler)
{
    DownsamplerReportMap result {};
    result.reserve(reads.size());
    for (auto& p : reads) {
        result.emplace(p.first, downsampler.downsample(p.second));
    }
    return result;
}

std::size_t count_downsampled_reads(const DownsamplerReportMap& reports);
std::size_t count_downsampled_reads(const DownsamplerReportMap& reports, const GenomicRegion& region);

} // namespace readpipe
} // namespace octopus

#endif
