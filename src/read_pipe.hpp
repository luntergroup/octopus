//
//  read_pipe.hpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef read_pipe_hpp
#define read_pipe_hpp

#include <vector>
#include <unordered_map>
#include <cstddef>

#include <boost/optional.hpp>

#include "common.hpp"
#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "read_filter.hpp"
#include "read_transform.hpp"
#include "downsampler.hpp"

/*
 ReadPipe provides a wrapper for the basic manipulation of AlignedRead's. It is responsable for
 actually fetching the reads from file, and applying any filters, transforms etc. The result is 
 then a set of reads that can be used for calling.
 
 Rather than just fetching all reads for all samples in one go and processing, we can be a bit more 
 clever and fetch reads in batches (e.g. each sample, or regions in samples). This could potentially 
 decrease avergae memory consumption (and also increase runtime performance) by only minimising the
 number of 'bad' reads in memory. If we are really short on memory we could even compress filtered 
 read batches while we process other batches.
 */

namespace Octopus {

class ReadPipe
{
public:
    ReadPipe() = delete;
    explicit ReadPipe(ReadManager& read_manager,
                      ReadFilterer read_filter,
                      boost::optional<Downsampler> downsampler,
                      ReadTransform read_transform);
    ~ReadPipe() = default;
    
    ReadPipe(const ReadPipe&)            = delete;
    ReadPipe& operator=(const ReadPipe&) = delete;
    ReadPipe(ReadPipe&&)                 = default;
    ReadPipe& operator=(ReadPipe&&)      = default;
    
    ReadMap fetch_reads(std::vector<SampleIdType> samples, const GenomicRegion& region);
    
private:
    ReadManager& read_manager_;
    ReadFilterer read_filter_;
    boost::optional<Downsampler> downsampler_;
    ReadTransform read_transform_;
    
    using ContigIdType = GenomicRegion::ContigNameType;
    
    struct ContigStats
    {
        size_t min_coverage_;
        size_t max_coverage_;
        size_t average_coverage_;
        size_t max_filtered_reads_;
        size_t min_filtered_reads_;
        size_t average_filtered_reads_;
    };
    
    using ContigStatsMap = std::unordered_map<ContigIdType, ContigStats>;
    using SampleStatsMap = std::unordered_map<SampleIdType, ContigStatsMap>;
    
    SampleStatsMap sample_stats_;
    
    double average_coverage(const SampleIdType& sample) const noexcept;
    double average_coverage(const SampleIdType& sample, const ContigIdType& contig) const noexcept;
};
    
} // namespace Octopus

#endif /* read_pipe_hpp */
