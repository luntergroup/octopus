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
#include "coverage_tracker.hpp"

namespace octopus {

namespace {

auto draw_sample(const InputRegionMap& regions, std::discrete_distribution<>& contig_sampling_distribution)
{
    static std::mt19937 generator {42};
    return std::next(std::cbegin(regions), contig_sampling_distribution(generator));
}

auto choose_sample_window(const GenomicRegion& target)
{
    static std::mt19937 generator {42};
    std::uniform_int_distribution<GenomicRegion::Position> dist {target.begin(), target.end()};
    return GenomicRegion {target.contig_name(), dist(generator), target.end()};
}

auto choose_sample_region(const SampleName& sample, const InputRegionMap::mapped_type& regions)
{
    assert(!regions.empty());
    return choose_sample_window(*random_select(std::cbegin(regions), std::cend(regions)));
}

auto choose_sample_region(const SampleName& sample, const InputRegionMap& regions,
                          std::discrete_distribution<>& contig_sampling_distribution)
{
    return choose_sample_region(sample, draw_sample(regions, contig_sampling_distribution)->second);
}

struct SamplingSummary
{
    InputRegionMap sampled_regions;
    std::size_t num_samples;
};

boost::optional<GenomicRegion>
choose_next_sample_region(const SampleName& sample,
                          const InputRegionMap& regions,
                          const ReadSetProfileConfig& config,
                          std::discrete_distribution<>& contig_sampling_distribution,
                          SamplingSummary& sampling_summary)
{
    if (!regions.empty() && sampling_summary.num_samples < std::max(config.max_draws_per_sample, regions.size() * config.min_draws_per_contig)) {
        for (const auto& p : regions) {
            if (!p.second.empty() && sampling_summary.sampled_regions[p.first].size() < config.min_draws_per_contig) {
                auto sample_region = choose_sample_region(sample, p.second);
                sampling_summary.sampled_regions[sample_region.contig_name()].insert(sample_region);
                return sample_region;
            }
        }
        return choose_sample_region(sample, regions, contig_sampling_distribution);
    } else {
        return boost::none;
    }
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

void cut(const GenomicRegion& region, MappableFlatSet<GenomicRegion>& regions)
{
    const auto overlapped = bases(overlap_range(regions, region));
    std::vector<GenomicRegion> overlapped_copy(std::cbegin(overlapped), std::cend(overlapped));
    regions.erase(std::cbegin(overlapped), std::cend(overlapped));
    for (const auto& overlap : overlapped_copy) {
        if (!contains(region, overlap)) {
            if (begins_before(overlap, region)) {
                regions.insert(left_overhang_region(overlap, region));
            }
            if (ends_before(region, overlap)) {
                regions.insert(right_overhang_region(overlap, region));
            }
        }
    }
}

template <typename Range, typename T>
void fill_summary_stats(const Range& values, ReadSetProfile::SummaryStats<T>& result)
{
    if (!values.empty()) {
        result.min = *std::min_element(std::cbegin(values), std::cend(values));
        result.max = *std::max_element(std::cbegin(values), std::cend(values));
        result.mean = maths::mean(values);
        result.median = maths::median(values);
        result.stdev = maths::stdev(values);
    }
}

void fill_summary_stats(const std::deque<MemoryFootprint>& footprints, ReadSetProfile::ReadMemoryStats& result)
{
    std::vector<std::size_t> bytes(footprints.size());
    std::transform(std::cbegin(footprints), std::cend(footprints), std::begin(bytes), [] (auto footprint) { return footprint.bytes(); });
    fill_summary_stats(bytes, result);
}

template <typename Range>
auto copy_positive(const Range& values)
{
    using T = typename Range::value_type;
    std::vector<T> result {};
    result.reserve(values.size());
    std::copy_if(std::cbegin(values), std::cend(values), std::back_inserter(result), [] (T value) { return value > 0; });
    return result;
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

template <typename DepthType>
boost::optional<ReadSetProfile>
profile_reads_helper(const std::vector<SampleName>& samples,
                     const InputRegionMap& regions,
                     const ReadManager& source,
                     ReadSetProfileConfig config)
{
    ReadSetProfile result {};
    std::deque<MemoryFootprint> memory_footprints {}, fragmented_memory_footprints {};
    std::unordered_map<GenomicRegion::ContigName, std::vector<DepthType>> contig_depths {};
    std::deque<unsigned> read_lengths {};
    std::deque<AlignedRead::MappingQuality> mapping_qualities {};
    const auto contig_sampling_distribution = make_contig_sampling_distribution(regions);
    for (const auto& sample : samples) {
        std::vector<DepthType> sample_depths {};
        std::unordered_map<GenomicRegion::ContigName, std::vector<DepthType>> sample_contig_depths {};
        SamplingSummary sampling_summary {};
        auto remaining_sampling_regions = regions;
        auto sample_contig_sampling_distribution = contig_sampling_distribution;
        while (true) {
            const auto target_sampling_region = choose_next_sample_region(sample, remaining_sampling_regions, config, sample_contig_sampling_distribution, sampling_summary);
            if (!target_sampling_region) break;
            CoverageTracker<GenomicRegion, DepthType> depth_tracker {true};
            auto remaining_reads = static_cast<int>(config.max_reads_per_draw);
            const auto read_visitor = [&] (const SampleName& sample, AlignedRead read) {
                read_lengths.push_back(sequence_size(read));
                mapping_qualities.push_back(read.mapping_quality());
                memory_footprints.push_back(footprint(read));
                if (config.fragment_size) {
                    fragmented_memory_footprints.push_back(fragmented_footprint(read, *config.fragment_size));
                }
                depth_tracker.add(read);
                return --remaining_reads > 0;
            };
            source.iterate(sample, *target_sampling_region, read_visitor);
            auto sampled_region = *target_sampling_region;
            if (depth_tracker.any()) {
                auto sampled_reads_region = *depth_tracker.encompassing_region();
                assert(!read_lengths.empty());
                if (remaining_reads > 0) {
                    sampled_region = *target_sampling_region;
                } else {
                    assert(!is_before(sampled_reads_region, *target_sampling_region));
                    sampled_region = closed_region(*target_sampling_region, sampled_reads_region);
                    if (size(sampled_region) > read_lengths.back()) {
                        // Ignore the last read length bases to avoid adding positions undersampled because
                        // the sampled read limit was hit.
                        sampled_region = expand_rhs(sampled_region, -static_cast<GenomicRegion::Distance>(read_lengths.back()));
                    } else {
                        sampled_region = expand_rhs(head_region(*target_sampling_region), read_lengths.back() / 2);
                    }
                }
            }
            auto read_depths = depth_tracker.get(sampled_region);
            utils::append(read_depths, sample_contig_depths[sampled_region.contig_name()]);
            utils::append(std::move(read_depths), sample_depths);
            ++sampling_summary.num_samples;
            auto removal_region = sampled_region;
            if (depth_tracker.any()) {
                removal_region = encompassing_region(removal_region, *depth_tracker.encompassing_region());
            }
            cut(removal_region, remaining_sampling_regions.at(sampled_region.contig_name()));
            if (remaining_sampling_regions.at(sampled_region.contig_name()).empty()) {
                remaining_sampling_regions.erase(sampled_region.contig_name());
                sample_contig_sampling_distribution = make_contig_sampling_distribution(remaining_sampling_regions);
            }
        }
        ReadSetProfile::GenomeContigDepthStatsPair sample_depth_stats {};
        if (!sample_depths.empty()) {
            std::sort(std::begin(sample_depths), std::end(sample_depths)); // sorting means no copying from stats calculations
            fill_depth_stats(sample_depths, sample_depth_stats.genome);
            sample_depths.clear();
            sample_depths.shrink_to_fit();
            for (auto& p : sample_contig_depths) {
                std::sort(std::begin(p.second), std::end(p.second)); // sorting means no copying from stats calculations
                sample_depth_stats.contig.emplace(p.first, make_depth_stats(p.second));
                utils::append(std::move(p.second), contig_depths[p.first]);
                p.second.clear();
                p.second.shrink_to_fit();
            }
        }
        result.depth_stats.sample.emplace(sample, std::move(sample_depth_stats));
    }
    if (memory_footprints.empty()) return boost::none;
    fill_summary_stats(memory_footprints, result.memory_stats);
    if (config.fragment_size) {
        result.fragmented_memory_stats = ReadSetProfile::ReadMemoryStats {};
        fill_summary_stats(fragmented_memory_footprints, *result.fragmented_memory_stats);
    }
    fill_summary_stats(read_lengths, result.length_stats);
    fill_summary_stats(mapping_qualities, result.mapping_quality_stats);
    std::size_t total_depths {0};
    for (auto& p : contig_depths) {
        std::sort(std::begin(p.second), std::end(p.second)); // sorting means no copying from stats calculations
        result.depth_stats.combined.contig.emplace(p.first, make_depth_stats(p.second));
        total_depths += p.second.size();
    }
    if (total_depths > 0) {
        std::vector<DepthType> depths {};
        depths.reserve(total_depths / contig_depths.size());
        for (auto& p : contig_depths) {
            utils::append(std::move(p.second), depths);
            p.second.clear();
            p.second.shrink_to_fit();
        }
        contig_depths.clear();
        std::sort(std::begin(depths), std::end(depths)); // sorting means no copying from stats calculations
        fill_depth_stats(depths, result.depth_stats.combined.genome);
    }
    return result;
}

} // namespace

boost::optional<ReadSetProfile>
profile_reads(const std::vector<SampleName>& samples,
              const InputRegionMap& regions,
              const ReadManager& source,
              ReadSetProfileConfig config)
{
    boost::optional<ReadSetProfile> result {};
    try {
        result = profile_reads_helper<unsigned short>(samples, regions, source, config);
    } catch (const std::runtime_error& e) {
        try {
            result = profile_reads_helper<unsigned>(samples, regions, source, config);
        } catch (const std::runtime_error& e) {
            result = profile_reads_helper<unsigned long>(samples, regions, source, config);
        }
    }
    return result;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const ReadSetProfile::SummaryStats<T>& stats)
{
    return os << "min: " << stats.min << " max: " << stats.max << " mean: " << stats.mean << " median: " << stats.median << " stdev: " << stats.stdev;
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
