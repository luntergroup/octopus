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

auto choose_sample_region(const GenomicRegion& from, GenomicRegion::Size max_size)
{
    if (size(from) <= max_size) return from;
    const auto max_begin = from.end() - max_size;
    static std::mt19937 generator {42};
    std::uniform_int_distribution<GenomicRegion::Position> dist {from.begin(), max_begin};
    return GenomicRegion {from.contig_name(), dist(generator), from.end()};
}

template <typename Range>
auto draw_sample(const SampleName& sample, const Range& regions,
                 const ReadManager& source, const ReadSetProfileConfig& config)
{
    const auto region_itr = random_select(std::cbegin(regions), std::cend(regions));
    const auto sample_region = choose_sample_region(*region_itr, 10 * config.max_reads_per_draw);
    auto test_region = source.find_covered_subregion(sample, sample_region, config.max_reads_per_draw);
    if (is_empty(test_region)) {
        test_region = expand_rhs(test_region, 1);
    }
    return source.fetch_reads(sample, test_region);
}

auto draw_sample(const SampleName& sample, const InputRegionMap& regions,
                 const ReadManager& source, const ReadSetProfileConfig& config,
                 std::discrete_distribution<>& contig_sampling_distribution)
{
    return draw_sample(sample, draw_sample(regions, contig_sampling_distribution)->second, source, config);
}

template <typename Range>
auto draw_sample_from_begin(const SampleName& sample, const Range& regions,
                            const ReadManager& source, const ReadSetProfileConfig& config)
{
    const auto region_itr = random_select(std::cbegin(regions), std::cend(regions));
    auto test_region = source.find_covered_subregion(sample, *region_itr, config.max_reads_per_draw);
    if (is_empty(test_region)) {
        test_region = expand_rhs(test_region, 1);
    }
    return source.fetch_reads(sample, test_region);
}

auto draw_sample_from_begin(const SampleName& sample, const InputRegionMap& regions,
                            const ReadManager& source, const ReadSetProfileConfig& config,
                            std::discrete_distribution<>& contig_sampling_distribution)
{
    return draw_sample_from_begin(sample, draw_sample(regions, contig_sampling_distribution)->second, source, config);
}

using ReadSetSamples = std::vector<ReadManager::ReadContainer>;

bool all_empty(const ReadSetSamples& samples)
{
    return std::all_of(std::cbegin(samples), std::cend(samples), [] (const auto& reads) { return reads.empty(); });
}

auto draw_samples(const SampleName& sample, const InputRegionMap& regions,
                  const ReadManager& source, const ReadSetProfileConfig& config,
                  std::discrete_distribution<>& contig_sampling_distribution)
{
    ReadSetSamples result {};
    result.reserve(config.max_draws_per_sample);
    auto remaining_draws = config.max_draws_per_sample;
    // Draw from each contig first to ensure all contigs get sampled
    for (const auto& p : regions) {
	    std::generate_n(std::back_inserter(result), config.min_draws_per_contig,
	                    [&] () { return draw_sample(sample, p.second, source, config); });
        if (remaining_draws > 0) --remaining_draws;
    }
    // Then sample contigs randomly
    std::generate_n(std::back_inserter(result), remaining_draws,
                    [&] () { return draw_sample(sample, regions, source, config, contig_sampling_distribution); });
    if (all_empty(result)) {
        result.back() = draw_sample_from_begin(sample, regions, source, config, contig_sampling_distribution);
    }
    return result;
}

auto draw_samples(const std::vector<SampleName>& samples, const InputRegionMap& regions,
                  const ReadManager& source, const ReadSetProfileConfig& config)
{
    std::vector<unsigned> contig_weights(regions.size());
    std::transform(std::cbegin(regions), std::cend(regions), std::begin(contig_weights),
                   [] (const auto& p) { return sum_region_sizes(p.second); });
    std::discrete_distribution<> contig_sampling_distribution(std::cbegin(contig_weights), std::cend(contig_weights));
    std::vector<ReadSetSamples> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.push_back(draw_samples(sample, regions, source, config, contig_sampling_distribution));
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
    const auto read_sets = draw_samples(samples, sampling_regions, source, config);
    if (read_sets.empty()) return boost::none;
    const auto bytes = get_read_bytes(read_sets);
    if (bytes.empty()) return boost::none;
    ReadSetProfile result {};
    fill_summary_stats(bytes, result.memory_stats);
    if (config.fragment_size) {
        const auto fragmented_bytes = compute_fragmented_template_bytes(read_sets, *config.fragment_size);
        result.fragmented_memory_stats = ReadSetProfile::ReadMemoryStats {};
        fill_summary_stats(fragmented_bytes, *result.fragmented_memory_stats);
    }
    std::deque<unsigned> depths {};
    std::unordered_map<GenomicRegion::ContigName, std::deque<unsigned>> contig_depths {};
    std::vector<unsigned> read_lengths {};
    std::vector<AlignedRead::MappingQuality> mapping_qualities {};
    for (std::size_t s {0}; s < samples.size(); ++s) {
        std::deque<unsigned> sample_depths {};
        std::unordered_map<GenomicRegion::ContigName, std::deque<unsigned>> sample_contig_depths {};
        for (const auto& reads : read_sets[s]) {
            if (!reads.empty()) {
                auto read_depths = calculate_positional_coverage(reads);
                utils::append(read_depths, sample_contig_depths[contig_name(reads.front())]);
                utils::append(std::move(read_depths), sample_depths);
                read_lengths.reserve(read_lengths.size() + reads.size());
                std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(read_lengths),
                               [] (const auto& read) { return sequence_size(read); });
                std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(mapping_qualities),
                               [] (const auto& read) { return read.mapping_quality(); });
            }
        }
        ReadSetProfile::GenomeContigDepthStatsPair sample_depth_stats {};
        if (!sample_depths.empty()) {
            fill_depth_stats(sample_depths, sample_depth_stats.genome);
            for (const auto& p : sample_contig_depths) {
                sample_depth_stats.contig.emplace(p.first, make_depth_stats(p.second));
            }
        }
        result.depth_stats.sample.emplace(samples[s], std::move(sample_depth_stats));
        utils::append(std::move(sample_depths), depths);
        for (auto& p : sample_contig_depths) {
            utils::append(std::move(p.second), contig_depths[p.first]);
        }
    }
    if (!depths.empty()) {
        fill_depth_stats(depths, result.depth_stats.combined.genome);
        for (const auto& p : contig_depths) {
            result.depth_stats.combined.contig.emplace(p.first, make_depth_stats(p.second));
        }
    }
    if (!read_lengths.empty()) {
        fill_summary_stats(read_lengths, result.length_stats);
    }
    if (!mapping_qualities.empty()) {
        fill_summary_stats(mapping_qualities, result.mapping_quality_stats);
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
