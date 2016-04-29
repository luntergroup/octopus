//
//  read_pipe.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "read_pipe.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <cassert>

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
    
std::vector<GenomicRegion> join_close_regions(const std::vector<GenomicRegion>& regions,
                                              const GenomicRegion::SizeType max_distance)
{
    std::vector<GenomicRegion> result {};
    
    if (regions.empty()) return result;
    
    result.reserve(regions.size());
    
    auto tmp = regions.front();
    
    std::for_each(std::next(std::cbegin(regions)), std::cend(regions),
                  [&tmp, &result, max_distance] (const auto& region) {
                      if (inner_distance(tmp, region) <= max_distance) {
                          tmp = encompassing_region(tmp, region);
                      } else {
                          result.emplace_back(tmp);
                          tmp = region;
                      }
                  });
    
    if (result.empty() || ends_equal(tmp, result.back())) {
        result.emplace_back(std::move(tmp));
    }
    
    return result;
}

ReadMap ReadPipe::fetch_reads(const std::vector<GenomicRegion>& regions)
{
    assert(std::is_sorted(std::cbegin(regions), std::cend(regions)));
    
    const auto fetch_regions = join_close_regions(extract_covered_regions(regions), 500);
    
    ReadMap result {samples_.size()};
    
    for (const auto& sample : samples_) {
        result[sample].reserve(100000);
    }
    
    for (const auto& region : fetch_regions) {
        auto tmp = fetch_reads(region);
        
        for (const auto& sample : samples_) {
            result.at(sample).insert(std::make_move_iterator(std::begin(tmp.at(sample))),
                                     std::make_move_iterator(std::end(tmp.at(sample))));
            tmp.at(sample).clear();
            tmp.at(sample).shrink_to_fit();
        }
    }
    
    for (const auto& sample : samples_) {
        result[sample].shrink_to_fit();
    }
    
    return result;
}

} // namespace Octopus
