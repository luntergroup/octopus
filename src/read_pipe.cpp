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

namespace Octopus
{

// public members

ReadPipe::ReadPipe(ReadManager& read_manager, ReadTransform read_transform, ReadFilterer read_filter,
                   boost::optional<Downsampler> downsampler, std::vector<SampleIdType> samples)
:
read_manager_ {read_manager},
read_filter_ {std::move(read_filter)},
downsampler_ {std::move(downsampler)},
read_transform_ {std::move(read_transform)},
samples_ {std::move(samples)},
debug_log_ {}
{
    if (DEBUG_MODE) debug_log_ = Logging::DebugLogger {};
}

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

const std::vector<SampleIdType>& ReadPipe::samples() const noexcept
{
    return samples_;
}

namespace
{
    template <typename Container>
    void move_construct(Container&& src, ReadMap::mapped_type& dst)
    {
        dst.insert(std::make_move_iterator(std::begin(src)),
                   std::make_move_iterator(std::end(src)));
    }
    
    template <>
    void move_construct<ReadMap::mapped_type>(ReadMap::mapped_type&& src, ReadMap::mapped_type& dst)
    {
        dst = std::move(src);
    }
    
    template <typename Container>
    void insert_each(Container&& src, ReadMap& dst)
    {
        for (auto& p : src) {
            auto& sample_dst = dst.at(p.first);
            
            if (sample_dst.empty()) {
                move_construct(std::move(p.second), sample_dst);
            } else {
                sample_dst.insert(std::make_move_iterator(std::begin(p.second)),
                                  std::make_move_iterator(std::end(p.second)));
            }
            
            p.second.clear();
            p.second.shrink_to_fit();
        }
        
        src.clear();
    }
    
    void shrink_to_fit(ReadMap& reads)
    {
        for (auto& p : reads) {
            p.second.shrink_to_fit();
        }
    }
} // namespace

ReadMap ReadPipe::fetch_reads(const GenomicRegion& region)
{
    ReadMap result {samples_.size()};
    
    for (const auto& sample : samples_) {
        result.emplace(std::piecewise_construct, std::forward_as_tuple(sample), std::forward_as_tuple());
    }
    
    const auto batches = batch_samples(samples_);
    
    for (const auto& batch : batches) {
        auto batch_reads = read_manager_.get().fetch_reads(batch, region);
        
        if (debug_log_) {
            stream(*debug_log_) << "Fetched " << count_reads(batch_reads) << " unfiltered reads from " << region;
        }
        
        // transforms should be done first as they may affect which reads are filtered
        transform_reads(batch_reads, read_transform_);
        
        if (debug_log_) {
            SampleFilterCountMap<SampleIdType, decltype(read_filter_)> filter_counts {};
            filter_counts.reserve(samples_.size());
            
            for (const auto& sample : samples_) {
                filter_counts[sample].reserve(read_filter_.num_filters());
            }
            
            erase_filtered_reads(batch_reads, filter(batch_reads, read_filter_, filter_counts));
            
            for (const auto& p : filter_counts) {
                stream(*debug_log_) << "In sample " << p.first;
                for (const auto& c : p.second) {
                    stream(*debug_log_) << c.second << " reads were removed by the " << c.first << " filter";
                }
            }
        } else {
            erase_filtered_reads(batch_reads, filter(batch_reads, read_filter_));
        }
        
        if (debug_log_) {
            stream(*debug_log_) << "There are " << count_reads(batch_reads) << " reads in " << region
                            << " after filtering";
        }
        
        if (downsampler_) {
            auto reads = make_mappable_map(std::move(batch_reads));
            
            const auto n = downsampler_->downsample(reads);
            
            if (debug_log_) stream(*debug_log_) << "Downsampling removed " << n << " reads from " << region;
            
            insert_each(std::move(reads), result);
        } else {
            insert_each(std::move(batch_reads), result);
        }
    }
    
    shrink_to_fit(result); // TODO: should we make this conditional on extra capacity?
    
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
        result.emplace(std::piecewise_construct, std::forward_as_tuple(sample), std::forward_as_tuple());
    }
    
    for (const auto& region : fetch_regions) {
        insert_each(fetch_reads(region), result);
    }
    
    shrink_to_fit(result);
    
    return result;
}

} // namespace Octopus
