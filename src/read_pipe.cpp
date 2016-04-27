//
//  read_pipe.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "read_pipe.hpp"

#include <utility>

#include "read_utils.hpp"
#include "mappable_algorithms.hpp"
#include "logging.hpp"

#include <iostream>
#include "timers.hpp"

namespace Octopus
{

// public members

ReadPipe::ReadPipe(ReadManager& read_manager, ReadFilterer read_filter,
                   boost::optional<Downsampler> downsampler, ReadTransform read_transform,
                   std::vector<SampleIdType> samples)
:
read_manager_ {read_manager},
read_filter_ {std::move(read_filter)},
downsampler_ {std::move(downsampler)},
read_transform_ {std::move(read_transform)},
samples_ {std::move(samples)}
{}

std::vector<std::vector<SampleIdType>> batch_samples(std::vector<SampleIdType> samples)
{
    std::vector<std::vector<SampleIdType>> result {};
    
    result.emplace_back(std::move(samples)); // TODO: find a better strategy for this
    
    return result;
}

void ReadPipe::set_read_manager(ReadManager& read_manager) noexcept
{
    read_manager_ = read_manager;
}

unsigned ReadPipe::num_samples() const noexcept
{
    return static_cast<unsigned>(samples_.size());
}

const std::vector<SampleIdType>& ReadPipe::get_samples() const noexcept
{
    return samples_;
}

ReadMap ReadPipe::fetch_reads(const GenomicRegion& region)
{
    ReadMap result {samples_.size()};
    
    const auto batches = batch_samples(samples_);
    
    for (const auto& batch : batches) {
        auto batch_reads = read_manager_.get().fetch_reads(batch, region);
        
        if (DEBUG_MODE) {
            Logging::DebugLogger log {};
            stream(log) << "Fetched " << count_reads(batch_reads) << " unfiltered reads from " << region;
        }
        
        erase_filtered_reads(batch_reads, filter(batch_reads, read_filter_));
        
        if (DEBUG_MODE) {
            Logging::DebugLogger log {};
            stream(log) << "There are " << count_reads(batch_reads) << " reads in " << region
                << " after filtering";
        }
        
        if (downsampler_) {
            const auto n = downsampler_->downsample(batch_reads);
            
            if (DEBUG_MODE) {
                Logging::DebugLogger log {};
                stream(log) << "Downsampling removed " << n << " reads from " << region;
            }
        }
        
        transform_reads(batch_reads, read_transform_);
        
        for (auto&& sample_batch : batch_reads) {
            result.emplace(sample_batch.first, std::move(sample_batch.second));
        }
    }
    
    return result;
}

ReadMap ReadPipe::fetch_reads(const std::vector<GenomicRegion>& regions)
{
    // TODO: improve this
    return fetch_reads(encompassing_region(regions));
}

} // namespace Octopus
