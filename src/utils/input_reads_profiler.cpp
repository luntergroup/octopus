// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "input_reads_profiler.hpp"

#include <random>
#include <deque>
#include <iterator>
#include <algorithm>
#include <utility>
#include <cassert>
#include <iostream>

#include "mappable_algorithms.hpp"
#include "maths.hpp"
#include "append.hpp"
#include "random_select.hpp"
#include "read_stats.hpp"

namespace octopus {

namespace {

auto get_covered_sample_regions(const std::vector<SampleName>& samples, const InputRegionMap& input_regions,
                                const ReadManager& read_manager)
{
    InputRegionMap result {};
    result.reserve(input_regions.size());
    for (const auto& p : input_regions) {
        InputRegionMap::mapped_type contig_regions {};
        std::copy_if(std::cbegin(p.second), std::cend(p.second),
                     std::inserter(contig_regions, std::begin(contig_regions)),
                     [&] (const auto& region) { return read_manager.has_reads(samples, region); });
        if (!contig_regions.empty()) {
            result.emplace(p.first, std::move(contig_regions));
        }
    }
    return result;
}

auto choose_sample_region(const GenomicRegion& from, GenomicRegion::Size max_size)
{
    if (size(from) <= max_size) return from;
    const auto max_begin = from.end() - max_size;
    static std::mt19937 generator {42};
    std::uniform_int_distribution<GenomicRegion::Position> dist {from.begin(), max_begin};
    return GenomicRegion {from.contig_name(), dist(generator), from.end()};
}

auto draw_sample(const SampleName& sample, const InputRegionMap& regions,
                 const ReadManager& source, const ReadSetProfileConfig& config)
{
    const auto contig_itr = random_select(std::cbegin(regions), std::cend(regions));
    assert(!contig_itr->second.empty());
    const auto region_itr = random_select(std::cbegin(contig_itr->second), std::cend(contig_itr->second));
    const auto sample_region = choose_sample_region(*region_itr, 10 * config.max_sample_size);
    auto test_region = source.find_covered_subregion(sample, sample_region, config.max_sample_size);
    if (is_empty(test_region)) {
        test_region = expand_rhs(test_region, 1);
    }
    return source.fetch_reads(sample, test_region);
}

auto draw_sample_from_begin(const SampleName& sample, const InputRegionMap& regions,
                            const ReadManager& source, const ReadSetProfileConfig& config)
{
    const auto contig_itr = random_select(std::cbegin(regions), std::cend(regions));
    assert(!contig_itr->second.empty());
    const auto region_itr = random_select(std::cbegin(contig_itr->second), std::cend(contig_itr->second));
    auto test_region = source.find_covered_subregion(sample, *region_itr, config.max_sample_size);
    if (is_empty(test_region)) {
        test_region = expand_rhs(test_region, 1);
    }
    return source.fetch_reads(sample, test_region);
}

using ReadSetSamples = std::vector<ReadManager::ReadContainer>;

bool all_empty(const ReadSetSamples& samples)
{
    return std::all_of(std::cbegin(samples), std::cend(samples), [] (const auto& reads) { return reads.empty(); });
}

auto draw_samples(const SampleName& sample, const InputRegionMap& regions,
                  const ReadManager& source, const ReadSetProfileConfig& config)
{
    ReadSetSamples result {};
    result.reserve(config.max_samples_per_sample);
    std::generate_n(std::back_inserter(result), config.max_samples_per_sample,
                    [&] () { return draw_sample(sample, regions, source, config); });
    if (all_empty(result)) {
        result.back() = draw_sample_from_begin(sample, regions, source, config);
    }
    return result;
}

auto draw_samples(const std::vector<SampleName>& samples, const InputRegionMap& regions,
                  const ReadManager& source, const ReadSetProfileConfig& config)
{
    std::vector<ReadSetSamples> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.push_back(draw_samples(sample, regions, source, config));
    }
    return result;
}

auto get_read_bytes(const std::vector<ReadSetSamples>& read_sets)
{
    std::deque<std::size_t> result {};
    for (const auto& set : read_sets) {
        for (const auto& reads : set) {
            std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(result),
                           [] (const auto& read) noexcept { return footprint(read).bytes(); });
        }
    }
    return result;
}

auto fragmented_footprint(const AlignedRead& read, const AlignedRead::NucleotideSequence::size_type fragment_size)
{
    const auto fragments = split(read, fragment_size);
    const static auto add_footprint = [] (auto total, const auto& read) { return total + footprint(read); };
    return std::accumulate(std::cbegin(fragments), std::cend(fragments), MemoryFootprint {0}, add_footprint);
}

auto compute_fragmented_template_bytes(const std::vector<ReadSetSamples>& read_sets, const AlignedRead::NucleotideSequence::size_type fragment_size)
{
    std::deque<std::size_t> result {};
    for (const auto& set : read_sets) {
        for (const auto& reads : set) {
            std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(result),
                           [fragment_size] (const auto& read) noexcept { return fragmented_footprint(read, fragment_size).bytes(); });
        }
    }
    return result;
}

auto calculate_median_fragmented_template_bytes(const std::vector<ReadSetSamples>& read_sets, const AlignedRead::NucleotideSequence::size_type fragment_size)
{
    return maths::median(compute_fragmented_template_bytes(read_sets, fragment_size));
}

template <typename T>
auto copy_positive(const std::deque<T>& values)
{
    std::vector<T> result {};
    result.reserve(values.size());
    std::copy_if(std::cbegin(values), std::cend(values), std::back_inserter(result), [] (T value) { return value > 0; });
    return result;
}

template <typename Range>
auto make_depth_distribution(const Range& depths)
{
    std::vector<unsigned> counts {};
    if (!depths.empty()) {
        auto max_depth = *std::max_element(std::cbegin(depths), std::cend(depths));
        counts.resize(max_depth + 1);
    }
    for (const auto& depth : depths) ++counts[depth];
    std::vector<double> result(counts.size());
    const auto calculate_frequency = [&] (auto count) noexcept { return static_cast<double>(count) / depths.size(); };
    std::transform(std::cbegin(counts), std::cend(counts), std::begin(result), calculate_frequency);
    return result;
}

} // namespace

boost::optional<ReadSetProfile>
profile_reads(const std::vector<SampleName>& samples,
              const InputRegionMap& input_regions,
              const ReadManager& source,
              ReadSetProfileConfig config)
{
    if (input_regions.empty()) return boost::none;
    const auto sampling_regions = get_covered_sample_regions(samples, input_regions, source);
    if (sampling_regions.empty()) return boost::none;
    const auto read_sets = draw_samples(samples, sampling_regions, source, config);
    if (read_sets.empty()) return boost::none;
    const auto bytes = get_read_bytes(read_sets);
    if (bytes.empty()) return boost::none;
    ReadSetProfile result {};
    result.mean_read_bytes = maths::mean(bytes);
    result.read_bytes_stdev = maths::stdev(bytes);
    if (config.fragment_size) {
        result.fragmented_template_median_bytes = calculate_median_fragmented_template_bytes(read_sets, *config.fragment_size);
    } else {
        result.fragmented_template_median_bytes = boost::none;
    }
    result.samples = samples;
    result.sample_depth_distribution.resize(samples.size());
    result.sample_mean_depth.resize(samples.size());
    result.sample_median_depth.resize(samples.size());
    result.sample_depth_stdev.resize(samples.size());
    result.sample_positive_depth_stdev.resize(samples.size());
    result.sample_median_positive_depth.resize(samples.size());
    result.sample_mean_positive_depth.resize(samples.size());
    std::deque<unsigned> depths {};
    std::vector<unsigned> read_lengths {};
    std::vector<AlignedRead::MappingQuality> mapping_qualities {};
    for (std::size_t s {0}; s < samples.size(); ++s) {
        std::deque<unsigned> sample_depths {};
        for (const auto& reads : read_sets[s]) {
            if (!reads.empty()) {
                utils::append(calculate_positional_coverage(reads), sample_depths);
                read_lengths.reserve(read_lengths.size() + reads.size());
                std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(read_lengths),
                               [] (const auto& read) { return sequence_size(read); });
                std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(mapping_qualities),
                               [] (const auto& read) { return read.mapping_quality(); });
            }
        }
        if (!sample_depths.empty()) {
            result.sample_depth_distribution[s] = make_depth_distribution(sample_depths);
            result.sample_mean_depth[s] = maths::mean(sample_depths);
            result.sample_median_depth[s] = maths::median(sample_depths);
            result.sample_depth_stdev[s] = maths::stdev(sample_depths);
            const auto sample_positive_depths = copy_positive(sample_depths);
            if (!sample_positive_depths.empty()) {
                result.sample_median_positive_depth[s] = maths::median(sample_positive_depths);
                result.sample_mean_positive_depth[s] = maths::mean(sample_positive_depths);
                result.sample_positive_depth_stdev[s] = maths::stdev(sample_positive_depths);
            } else {
                result.sample_median_positive_depth[s] = 0;
                result.sample_mean_positive_depth[s] = 0;
                result.sample_positive_depth_stdev[s] = 0;
            }
        } else {
            result.sample_mean_depth[s] = 0;
            result.sample_median_depth[s] = 0;
            result.sample_depth_stdev[s] = 0;
            result.sample_median_positive_depth[s] = 0;
        }
        utils::append(std::move(sample_depths), depths);
    }
    assert(!depths.empty());
    result.depth_distribution = make_depth_distribution(depths);
    result.mean_depth = maths::mean(depths);
    result.median_depth = maths::median(depths);
    result.depth_stdev = maths::stdev(depths);
    const auto positive_depths = copy_positive(depths);
    if (!positive_depths.empty()) {
        result.median_positive_depth = maths::median(positive_depths);
        result.mean_positive_depth = maths::mean(positive_depths);
        result.positive_depth_stdev = maths::stdev(positive_depths);
    } else {
        result.median_positive_depth = 0;
        result.mean_positive_depth = 0;
        result.positive_depth_stdev = 0;
    }
    result.max_read_length = *std::max_element(std::cbegin(read_lengths), std::cend(read_lengths));
    result.median_read_length = maths::median(read_lengths);
    if (!mapping_qualities.empty()) {
        result.max_mapping_quality = *std::max_element(std::cbegin(mapping_qualities), std::cend(mapping_qualities));
        result.median_mapping_quality = maths::median(mapping_qualities);
        result.rmq_mapping_quality = maths::rmq(mapping_qualities);
    } else {
        result.max_mapping_quality = 0;
        result.median_mapping_quality = 0;
        result.rmq_mapping_quality = 0;
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const ReadSetProfile& profile)
{
    os << "Bytes:\n";
    os << "\tmean = " << profile.mean_read_bytes << '\n';
    os << "\tstdev = " << profile.read_bytes_stdev << '\n';
    
    for (std::size_t s {0}; s < profile.sample_mean_depth.size(); ++s) {
        os << "Sample " << s << " depth:\n";
        os << "\tdistribution = "; for (auto p : profile.sample_depth_distribution[s]) os << p << ' '; os << '\n';
        os << "\tmean = " << profile.sample_mean_depth[s] << '\n';
        os << "\tmean positive = " << profile.sample_mean_positive_depth[s] << '\n';
        os << "\tmedian = " << profile.sample_median_depth[s] << '\n';
        os << "\tmedian positive = " << profile.sample_median_positive_depth[s] << '\n';
        os << "\tstdev = " << profile.sample_depth_stdev[s] << '\n';
        os << "\tpositive stdev = " << profile.sample_positive_depth_stdev[s] << '\n';
    }
    
    os << "Total depth:\n";
    os << "\tdistribution = "; for (auto p : profile.depth_distribution) os << p << ' '; os << '\n';
    os << "\tmean = " << profile.mean_depth << '\n';
    os << "\tmean positive = " << profile.mean_positive_depth << '\n';
    os << "\tmedian = " << profile.median_depth << '\n';
    os << "\tmedian positive = " << profile.median_positive_depth << '\n';
    os << "\tstdev = " << profile.depth_stdev << '\n';
    os << "\tpositive stdev = " << profile.positive_depth_stdev << '\n';
    
    os << "Read length:\n";
    os << "\tmax = " << profile.max_read_length << '\n';
    os << "\tmedian = " << profile.median_read_length << '\n';
    
    os << "Mapping quality:\n";
    os << "\tmax = " << static_cast<unsigned>(profile.max_mapping_quality + 33) << '\n';
    os << "\tmedian = " << static_cast<unsigned>(profile.median_mapping_quality + 33) << '\n';
    os << "\trmq = " << static_cast<unsigned>(profile.rmq_mapping_quality + 33) << '\n';
    
    return os;
}

} // namespace octopus
