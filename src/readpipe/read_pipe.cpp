// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_pipe.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <cassert>

#include "utils/read_stats.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/append.hpp"

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
, fragment_size_ {}
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
, fragment_size_ {}
, debug_log_ {}
{
    if (DEBUG_MODE) debug_log_ = logging::DebugLogger {};
}

ReadPipe::ReadPipe(const ReadManager& source, GenomicRegion::Size fragment_size, ReadTransformer prefilter_transformer,
                   ReadFilterer filterer, ReadTransformer postfilter_transformer,
                   boost::optional<Downsampler> downsampler, std::vector<SampleName> samples)
: source_ {source}
, prefilter_transformer_ {std::move(prefilter_transformer)}
, filterer_ {std::move(filterer)}
, postfilter_transformer_ {std::move(postfilter_transformer)}
, downsampler_ {std::move(downsampler)}
, samples_ {std::move(samples)}
, fragment_size_ {fragment_size}
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

template <typename Regions>
auto fetch_batch(const ReadManager& rm, const std::vector<SampleName>& samples, const Regions& region)
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

void fragment(std::vector<AlignedRead>& reads, const unsigned fragment_length)
{
    if (reads.empty()) return;
    std::vector<AlignedRead> fragments {};
    fragments.reserve((sequence_size(reads.front()) / fragment_length + 1) * reads.size());
    for (const auto& read : reads) {
        auto chunks = split(read, fragment_length);
        if (read.barcode().empty()) {
            for (auto& chunk : chunks) chunk.set_barcode(read.name());
        }
        utils::append(std::move(chunks), fragments);
    }
    reads = std::move(fragments);
}

void erase_non_overlapped(std::vector<AlignedRead>& reads, const GenomicRegion& region)
{
    const auto not_overlapped = [&] (const auto& read) { return !overlaps(read, region); };
    reads.erase(std::remove_if(std::begin(reads), std::end(reads), not_overlapped), std::end(reads));
}

void erase_non_overlapped(std::vector<AlignedRead>& reads, const std::vector<GenomicRegion>& regions)
{
    const auto not_overlapped = [&] (const auto& read) { 
        return std::none_of(std::cbegin(regions), std::cend(regions), [&] (const auto& region) { 
            return overlaps(read, region); }); };
    reads.erase(std::remove_if(std::begin(reads), std::end(reads), not_overlapped), std::end(reads));
}

template <typename Region>
void fragment(std::vector<AlignedRead>& reads, const unsigned fragment_length, const Region& region)
{
    fragment(reads, fragment_length);
    std::sort(std::begin(reads), std::end(reads));
    erase_non_overlapped(reads, region);
}

template <typename Region>
void fragment(ReadManager::SampleReadMap& reads, const unsigned fragment_length, const Region& region)
{
    for (auto& p : reads) fragment(p.second, fragment_length, region);
}

struct IsMappingQualityZero
{
    bool operator()(const AlignedRead& read) const noexcept { return read.mapping_quality() == 0; }
};

} // namespace

ReadMap ReadPipe::fetch_reads(const GenomicRegion& region, boost::optional<Report&> report) const
{
    using namespace readpipe;
    ReadMap result {samples_.size()};
    for (const auto& sample : samples_) {
        result.emplace(std::piecewise_construct, std::forward_as_tuple(sample), std::forward_as_tuple());
    }
    if (report) report->raw_depths.reserve(samples_.size());
    for (const auto& batch : batch_samples(samples_)) {
        auto batch_reads = fetch_batch(source_, batch, region);
        if (debug_log_) {
            stream(*debug_log_) << "Fetched " << count_reads(batch_reads) << " unfiltered reads from " << region;
        }
        if (report) {
            for (const auto& p : batch_reads) {
                report->raw_depths.emplace(p.first, make_coverage_tracker(p.second));
                report->mapping_quality_zero_depths.emplace(p.first, make_coverage_tracker(p.second, IsMappingQualityZero {}));
            }
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
                            stream(*debug_log_) << c.second << " reads failed the " << c.first << " filter";
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
        if (fragment_size_) {
            fragment(batch_reads, *fragment_size_, region);
            if (debug_log_) {
                stream(*debug_log_) << "Fragmented reads from " << region << " into " << count_reads(batch_reads) << " reads";
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
                                stream(*debug_log_) << c.second << " read fragments failed the " << c.first << " filter";
                            }
                        } else {
                            *debug_log_ << "No read fragments were filtered";
                        }
                    }
                }
            } else {
                erase_filtered_reads(batch_reads, filter(batch_reads, filterer_));
            }
        }
        if (downsampler_) {
            auto reads = make_mappable_map(std::move(batch_reads));
            auto downsample_reports = downsample(reads, *downsampler_);
            if (debug_log_) stream(*debug_log_) << "Downsampling removed " << count_downsampled_reads(downsample_reports) << " reads from " << region;
            if (report) {
                report->downsample_report = std::move(downsample_reports);
            }
            insert_each(std::move(reads), result);
        } else {
            insert_each(std::move(batch_reads), result);
        }
    }
    shrink_to_fit(result); // TODO: should we make this conditional on extra capacity?
    return result;
}

auto sort_and_merge(std::vector<GenomicRegion> regions)
{
    std::sort(std::begin(regions), std::end(regions));
    return extract_covered_regions(regions);
}

ReadMap ReadPipe::fetch_reads(const std::vector<GenomicRegion>& regions, boost::optional<Report&> report) const
{
    if (regions.size() == 1) { return fetch_reads(regions.front(), report); }
    const auto covered_regions = sort_and_merge(regions);
    if (covered_regions.size() == 1) { return fetch_reads(covered_regions.front(), report); }
    using namespace readpipe;
    ReadMap result {samples_.size()};
    for (const auto& sample : samples_) {
        result.emplace(std::piecewise_construct, std::forward_as_tuple(sample), std::forward_as_tuple());
    }
    if (regions.empty()) return result;
    if (report) report->raw_depths.reserve(samples_.size());
    for (const auto& batch : batch_samples(samples_)) {
        auto batch_reads = fetch_batch(source_, batch, covered_regions);
        // if (debug_log_) {
        //     stream(*debug_log_) << "Fetched " << count_reads(batch_reads) << " unfiltered reads from " << region;
        // }
        if (report) {
            for (const auto& p : batch_reads) {
                report->raw_depths.emplace(p.first, make_coverage_tracker(p.second));
                report->mapping_quality_zero_depths.emplace(p.first, make_coverage_tracker(p.second, IsMappingQualityZero {}));
            }
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
                            stream(*debug_log_) << c.second << " reads failed the " << c.first << " filter";
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
        // if (debug_log_) {
        //     stream(*debug_log_) << "There are " << count_reads(batch_reads) << " reads in " << region
        //                     << " after filtering";
        // }
        if (fragment_size_) {
            fragment(batch_reads, *fragment_size_, covered_regions);
            erase_filtered_reads(batch_reads, filter(batch_reads, filterer_));
        }
        if (downsampler_) {
            auto reads = make_mappable_map(std::move(batch_reads));
            auto downsample_reports = downsample(reads, *downsampler_);
            // if (debug_log_) stream(*debug_log_) << "Downsampling removed " << count_downsampled_reads(downsample_reports) << " reads from " << region;
            if (report) {
                report->downsample_report = std::move(downsample_reports);
            }
            insert_each(std::move(reads), result);
        } else {
            insert_each(std::move(batch_reads), result);
        }
    }
    shrink_to_fit(result); // TODO: should we make this conditional on extra capacity?
    return result;
}

} // namespace octopus
