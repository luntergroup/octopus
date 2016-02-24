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

#include <iostream> // DEBUG

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

ReadMap ReadPipe::fetch_reads(const GenomicRegion& region)
{
    ReadMap result {samples_.size()};
    
    const auto batches = batch_samples(samples_);
    
    for (const auto& batch : batches) {
        auto batch_reads = read_manager_.get().fetch_reads(batch, region);
        
        auto filtered_batch = filter_reads(std::move(batch_reads), read_filter_).first;
        
        if (downsampler_) {
            filtered_batch = downsampler_->sample(std::move(filtered_batch));
        }
        
        transform_reads(filtered_batch, read_transform_);
        
        for (auto&& sample_batch : filtered_batch) {
            result.emplace(std::move(sample_batch.first), std::move(sample_batch.second));
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
