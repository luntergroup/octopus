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

#include <iostream> // DEBUG

namespace Octopus {

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
    result.reserve(samples.size());
    
    for (auto&& sample : samples) {
        result.emplace_back(std::vector<SampleIdType> {std::move(sample)});
    }
    
    return result;
}

void ReadPipe::set_read_manager(ReadManager& read_manager) noexcept
{
    read_manager_ = read_manager;
}

ReadMap ReadPipe::fetch_reads(const GenomicRegion& region)
{
    ReadMap result {};
    result.reserve(samples_.size());
    
    //auto batches = batch_samples(std::move(samples));
    auto batches = std::vector<std::vector<SampleIdType>> {samples_};
    
    for (const auto& batch : batches) {
        auto batch_reads = read_manager_.get().fetch_reads(batch, region);
        
        //std::cout << "fetched " << count_reads(batch_reads) << " batch reads" << std::endl;
        
//        for (auto r : filter_reads(batch_reads, read_filter_).second.cbegin()->second) {
//            std::cout << r << std::endl;
//        }
//        exit(0);
        
        auto filtered_batch = filter_reads(std::move(batch_reads), read_filter_).first;
        
        //std::cout << "found " << count_reads(filtered_batch) << " good batch reads" << std::endl;
        
        if (downsampler_) {
            filtered_batch = downsampler_->sample(std::move(filtered_batch));
        }
        
        transform_reads(filtered_batch, read_transform_);
        
        for (auto&& sample_batch : filtered_batch) {
            result.emplace(std::move(sample_batch.first), std::move(sample_batch.second));
        }
    }
    
    std::cout << "fetched " << count_reads(result) << " total reads" << std::endl;
    
    return result;
}

} // namespace Octopus
