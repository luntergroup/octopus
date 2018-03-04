// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_pipe.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <cassert>

#include "utils/read_stats.hpp"
#include "utils/mappable_algorithms.hpp"

namespace octopus {

// public members

ReadPipe::ReadPipe(const ReadManager& source, std::vector<SampleName> samples)
: ReadPipe {source, {}, {}, boost::none, std::move(samples)}
{}

ReadPipe::ReadPipe(const ReadManager& source, ReadTransformer transformer, ReadFilterer filterer,
                   boost::optional<Downsampler> downsampler, std::vector<SampleName> samples)
: source_ {source}
, prefilter_transformer_ {std::move(transformer)}
, filterer_ {std::move(filterer)}
, postfilter_transformer_ {}
, downsampler_ {std::move(downsampler)}
, samples_ {std::move(samples)}
, debug_log_ {}
{
    if (DEBUG_MODE) debug_log_ = logging::DebugLogger {};
}

ReadPipe::ReadPipe(const ReadManager& source, ReadTransformer prefilter_transformer,
                   ReadFilterer filterer, ReadTransformer postfilter_transformer,
                   boost::optional<Downsampler> downsampler, std::vector<SampleName> samples)
: source_ {source}
, prefilter_transformer_ {std::move(prefilter_transformer)}
, filterer_ {std::move(filterer)}
, postfilter_transformer_ {std::move(postfilter_transformer)}
, downsampler_ {std::move(downsampler)}
, samples_ {std::move(samples)}
, debug_log_ {}
{
    if (DEBUG_MODE) debug_log_ = logging::DebugLogger {};
}

std::vector<std::vector<SampleName>> batch_samples(std::vector<SampleName> samples)
{
    std::vector<std::vector<SampleName>> result {};
    result.emplace_back(std::move(samples)); // TODO: find a better strategy for this
    return result;
}

const ReadManager& ReadPipe::read_manager() const noexcept
{
    return source_;
}

void ReadPipe::set_read_manager(const ReadManager& source) noexcept
{
    source_ = source;
}

unsigned ReadPipe::num_samples() const noexcept
{
    return static_cast<unsigned>(samples_.size());
}

const std::vector<SampleName>& ReadPipe::samples() const noexcept
{
    return samples_;
}

namespace {

template <typename Map>
void sort_each(Map& reads)
{
    for (auto& p : reads) {
        std::sort(std::begin(p.second), std::end(p.second));
    }
}

auto fetch_batch(const ReadManager& rm, const std::vector<SampleName>& samples, const GenomicRegion& region)
{
    auto result = rm.fetch_reads(samples, region);
    sort_each(result);
    return result;
}

template <typename Container>
void move_construct(Container&& src, ReadMap::mapped_type& dst)
{
    dst.insert(std::make_move_iterator(std::begin(src)), std::make_move_iterator(std::end(src)));
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

ReadMap ReadPipe::fetch_reads(const GenomicRegion& region) const
{
    using namespace readpipe;
    ReadMap result {samples_.size()};
    for (const auto& sample : samples_) {
        result.emplace(std::piecewise_construct, std::forward_as_tuple(sample), std::forward_as_tuple());
    }
    for (const auto& batch : batch_samples(samples_)) {
        auto batch_reads = fetch_batch(source_, batch, region);
        if (debug_log_) {
            stream(*debug_log_) << "Fetched " << count_reads(batch_reads) << " unfiltered reads from " << region;
        }
        transform_reads(batch_reads, prefilter_transformer_);
        if (debug_log_) {
            SampleFilterCountMap<SampleName, decltype(filterer_)> filter_counts {};
            filter_counts.reserve(samples_.size());
            for (const auto& sample : samples_) {
                filter_counts[sample].reserve(filterer_.num_filters());
            }
            erase_filtered_reads(batch_reads, filter(batch_reads, filterer_, filter_counts));
            if (filterer_.num_filters() > 0) {
                for (const auto& p : filter_counts) {
                    stream(*debug_log_) << "In sample " << p.first;
                    if (!p.second.empty()) {
                        for (const auto& c : p.second) {
                            stream(*debug_log_) << c.second << " failed the " << c.first << " filter";
                        }
                    } else {
                        *debug_log_ << "No reads were filtered";
                    }
                }
            }
        } else {
            erase_filtered_reads(batch_reads, filter(batch_reads, filterer_));
        }
        if (postfilter_transformer_) {
            transform_reads(batch_reads, *postfilter_transformer_);
        }
        if (debug_log_) {
            stream(*debug_log_) << "There are " << count_reads(batch_reads) << " reads in " << region
                            << " after filtering";
        }
        if (downsampler_) {
            auto reads = make_mappable_map(std::move(batch_reads));
            const auto n = downsample(reads, *downsampler_);
            if (debug_log_) stream(*debug_log_) << "Downsampling removed " << n << " reads from " << region;
            insert_each(std::move(reads), result);
        } else {
            insert_each(std::move(batch_reads), result);
        }
    }
    shrink_to_fit(result); // TODO: should we make this conditional on extra capacity?
    return result;
}

ReadMap ReadPipe::fetch_reads(const std::vector<GenomicRegion>& regions) const
{
    assert(std::is_sorted(std::cbegin(regions), std::cend(regions)));
    
    const auto covered_regions = extract_covered_regions(regions);
    const auto fetch_regions = join(covered_regions, 10000);
    
    ReadMap result {samples_.size()};
    const auto total_fetch_bp = sum_region_sizes(fetch_regions);
    for (const auto& sample : samples_) {
        const auto p = result.emplace(std::piecewise_construct, std::forward_as_tuple(sample),
                                      std::forward_as_tuple());
        p.first->second.reserve(20 * total_fetch_bp); // TODO: use estimated coverage
    }
    
    for (const auto& region : fetch_regions) {
        auto reads = fetch_reads(region);
        
        const auto request_regions = contained_range(covered_regions, region);
        const auto removal_regions = extract_intervening_regions(request_regions, region);
        
        std::for_each(std::crbegin(removal_regions), std::crend(removal_regions),
                      [&reads] (const auto& region) {
                          for (auto& p : reads) {
                              p.second.erase_contained(region);
                          }
                      });
        insert_each(std::move(reads), result);
    }
    
    shrink_to_fit(result);
    
    return result;
}

} // namespace octopus
