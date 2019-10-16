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

auto draw_sample(const InputRegionMap& regions, std::discrete_distribution<>& contig_sampling_distribution)
{
    static std::mt19937 generator {42};
    return std::next(std::cbegin(regions), contig_sampling_distribution(generator));
}

auto choose_sample_window(const GenomicRegion& target, const GenomicRegion::Size window_size)
{
    if (size(target) <= window_size) return target;
    const auto max_begin = target.end() - window_size;
    static std::mt19937 generator {42};
    std::uniform_int_distribution<GenomicRegion::Position> dist {target.begin(), max_begin};
    return GenomicRegion {target.contig_name(), dist(generator), target.end()};
}

auto choose_sample_region(const SampleName& sample, const InputRegionMap::mapped_type& regions,
                          const ReadManager& source, const ReadSetProfileConfig& config)
{
    const auto region_itr = random_select(std::cbegin(regions), std::cend(regions));
    const auto window = choose_sample_window(*region_itr, 100 * config.max_reads_per_draw);
    return source.find_covered_subregion(sample, window, config.max_reads_per_draw);
}

auto choose_sample_region(const SampleName& sample, const InputRegionMap& regions,
                          const ReadManager& source, const ReadSetProfileConfig& config,
                          std::discrete_distribution<>& contig_sampling_distribution)
{
    return choose_sample_region(sample, draw_sample(regions, contig_sampling_distribution)->second, source, config);
}

struct SamplingSummary
{
    InputRegionMap sampled_regions;
};

ReadManager::ReadContainer
draw_next_sample(const SampleName& sample, const InputRegionMap& regions,
                 const ReadManager& source, const ReadSetProfileConfig& config,
                 SamplingSummary& sampling_summary,
                 std::discrete_distribution<>& contig_sampling_distribution)
{
    for (const auto& p : regions) {
        if (sampling_summary.sampled_regions[p.first].size() < config.min_draws_per_contig) {
            auto sample_region = choose_sample_region(sample, p.second, source, config);
            sampling_summary.sampled_regions[sample_region.contig_name()].insert(sample_region);
            return source.fetch_reads(sample, sample_region);
        }
    }
    auto sample_region = choose_sample_region(sample, regions, source, config, contig_sampling_distribution);
    sampling_summary.sampled_regions[sample_region.contig_name()].insert(sample_region);
    return source.fetch_reads(sample, sample_region);
}

auto make_contig_sampling_distribution(const InputRegionMap& regions)
{
    std::vector<unsigned> contig_weights(regions.size());
    std::transform(std::cbegin(regions), std::cend(regions), std::begin(contig_weights),
                   [] (const auto& p) { return sum_region_sizes(p.second); });
    return std::discrete_distribution<> {std::cbegin(contig_weights), std::cend(contig_weights)};
}

auto fragmented_footprint(const AlignedRead& read, const AlignedRead::NucleotideSequence::size_type fragment_size)
{
    const auto fragments = split(read, fragment_size);
    const static auto add_footprint = [] (auto total, const auto& read) { return total + footprint(read); };
    return std::accumulate(std::cbegin(fragments), std::cend(fragments), MemoryFootprint {0}, add_footprint);
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

template <typename Range, typename T>
void fill_summary_stats(const Range& values, ReadSetProfile::SummaryStats<T>& result)
{
    assert(!values.empty());
    result.min = *std::min_element(std::cbegin(values), std::cend(values));
    result.max = *std::max_element(std::cbegin(values), std::cend(values));
    result.mean = maths::mean(values);
    result.median = maths::median(values);
    result.stdev = maths::stdev(values);
}

void fill_summary_stats(const std::deque<MemoryFootprint>& footprints, ReadSetProfile::ReadMemoryStats& result)
{
    std::vector<std::size_t> bytes(footprints.size());
    std::transform(std::cbegin(footprints), std::cend(footprints), std::begin(bytes), [] (auto footprint) { return footprint.bytes(); });
    fill_summary_stats(bytes, result);
}

template <typename Range>
void fill_depth_stats(const Range& depths, ReadSetProfile::DepthStats& result)
{
    result.distribution = make_depth_distribution(depths);
    fill_summary_stats(depths, result.all);
    fill_summary_stats(copy_positive(depths), result.positive);
}

template <typename Range>
auto make_depth_stats(const Range& depths)
{
    ReadSetProfile::DepthStats result {};
    fill_depth_stats(depths, result);
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
    
    ReadSetProfile result {};
    std::deque<MemoryFootprint> memory_footprints {}, fragmented_memory_footprints {};
    std::deque<unsigned> depths {};
    std::unordered_map<GenomicRegion::ContigName, std::deque<unsigned>> contig_depths {};
    std::deque<unsigned> read_lengths {};
    std::deque<AlignedRead::MappingQuality> mapping_qualities {};
    auto contig_sampling_distribution = make_contig_sampling_distribution(sampling_regions);
    
    for (const auto& sample : samples) {
        std::deque<unsigned> sample_depths {};
        std::unordered_map<GenomicRegion::ContigName, std::deque<unsigned>> sample_contig_depths {};
        SamplingSummary sampling_summary {};
        for (unsigned n {0}; n < config.max_draws_per_sample; ++n) {
            const auto reads = draw_next_sample(sample, sampling_regions, source, config, sampling_summary, contig_sampling_distribution);
            if (!reads.empty()) {
                auto read_depths = calculate_positional_coverage(reads);
                utils::append(read_depths, sample_contig_depths[contig_name(reads.front())]);
                utils::append(std::move(read_depths), sample_depths);
                for (const auto& read : reads) {
                    read_lengths.push_back(sequence_size(read));
                    mapping_qualities.push_back(read.mapping_quality());
                    memory_footprints.push_back(footprint(read));
                    if (config.fragment_size) {
                        fragmented_memory_footprints.push_back(fragmented_footprint(read, *config.fragment_size));
                    }
                }
            }
        }
        ReadSetProfile::GenomeContigDepthStatsPair sample_depth_stats {};
        if (!sample_depths.empty()) {
            fill_depth_stats(sample_depths, sample_depth_stats.genome);
            for (const auto& p : sample_contig_depths) {
                sample_depth_stats.contig.emplace(p.first, make_depth_stats(p.second));
            }
        }
        result.depth_stats.sample.emplace(sample, std::move(sample_depth_stats));
        utils::append(std::move(sample_depths), depths);
        for (auto& p : sample_contig_depths) {
            utils::append(std::move(p.second), contig_depths[p.first]);
        }
    }
    if (memory_footprints.empty()) return boost::none;
    fill_summary_stats(memory_footprints, result.memory_stats);
    if (config.fragment_size) {
        result.fragmented_memory_stats = ReadSetProfile::ReadMemoryStats {};
        fill_summary_stats(fragmented_memory_footprints, *result.fragmented_memory_stats);
    }
    fill_summary_stats(read_lengths, result.length_stats);
    fill_summary_stats(mapping_qualities, result.mapping_quality_stats);
    if (!depths.empty()) {
        fill_depth_stats(depths, result.depth_stats.combined.genome);
        for (const auto& p : contig_depths) {
            result.depth_stats.combined.contig.emplace(p.first, make_depth_stats(p.second));
        }
    }
    return result;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const ReadSetProfile::SummaryStats<T>& stats)
{
    return os << "min: " << stats.min << "mmax: " << stats.max << " mean: " << stats.mean << " median: " << stats.median << " stdev: " << stats.stdev;
}

std::ostream& operator<<(std::ostream& os, const ReadSetProfile::DepthStats::DiscreteDistribution& dist)
{
    for (auto depth : dist) os << depth << ' ';
    return os;
}

std::ostream& operator<<(std::ostream& os, const ReadSetProfile::DepthStats& stats)
{
    return os << " summary: all [" << stats.all << "] positive [" << stats.positive << ']'
              << " distribution: " << stats.distribution;
}

std::ostream& operator<<(std::ostream& os, const ReadSetProfile::GenomeContigDepthStatsPair& stats)
{
    os << "Genome: " << stats.genome << '\n';
    os << "Contigs:" << '\n';
    for (const auto& p : stats.contig) {
        os << '\t' << p.first << ": " << p.second << '\n';
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const ReadSetProfile::SampleCombinedDepthStatsPair& stats)
{
    os << "Combined: " << stats.combined << '\n';
    os << "Sample:" << '\n';
    for (const auto& p : stats.sample) {
        os << '\t' << p.first << ": " << p.second << '\n';
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const ReadSetProfile& profile)
{
    os << "Read memory stats: " << profile.memory_stats << '\n';
    if (profile.fragmented_memory_stats) {
        os << "Fragmented read memory stats: " << *profile.fragmented_memory_stats << '\n';
    }
    os << "Depth stats: " << profile.depth_stats << '\n';
    os << "Mapping quality stats: " << profile.mapping_quality_stats << '\n';
    os << "Read length stats: " << profile.length_stats << std::endl;
    return os;
}

} // namespace octopus
