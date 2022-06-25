// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "input_reads_profiler.hpp"

#include <random>
#include <deque>
#include <iterator>
#include <algorithm>
#include <utility>
#include <cassert>
#include <iostream>
#include <future>
#include <atomic>

#include <boost/random/uniform_int_distribution.hpp>

#include "mappable_algorithms.hpp"
#include "maths.hpp"
#include "append.hpp"
#include "random_select.hpp"
#include "read_stats.hpp"
#include "coverage_tracker.hpp"
#include "sequence_utils.hpp"

namespace octopus {

namespace {

auto draw_sample(const InputRegionMap& regions, std::discrete_distribution<>& contig_sampling_distribution)
{
    thread_local std::mt19937 generator {42};
    return std::next(std::cbegin(regions), contig_sampling_distribution(generator));
}

auto choose_sample_window(const GenomicRegion& target)
{
    thread_local std::mt19937 generator {42};
    boost::random::uniform_int_distribution<GenomicRegion::Position> dist {target.begin(), target.end()};
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
void erase_non_dna_or_rna_positions(std::vector<DepthType>& depths,
                                    const GenomicRegion& region,
                                    const ReferenceGenome& reference)
{
    const auto reference_sequence = reference.fetch_sequence(region);
    assert(depths.size() <= reference_sequence.size());
    auto reference_itr = std::cbegin(reference_sequence);
    const auto not_dna_or_rna = [&reference_itr] (DepthType) { return !utils::is_dna_or_rna_nucleotide(*reference_itr++); };
    depths.erase(std::remove_if(std::begin(depths), std::end(depths), not_dna_or_rna), std::end(depths));
}

template <typename DepthType>
struct SampleReadSetProfileHelper
{
    std::deque<MemoryFootprint> memory_footprints, fragmented_memory_footprints;
    std::unordered_map<GenomicRegion::ContigName, std::vector<DepthType>> contig_depths;
    std::deque<unsigned> read_lengths;
    std::deque<AlignedRead::MappingQuality> mapping_qualities;
    ReadSetProfile::GenomeContigDepthStatsPair depth_stats;
};

template <typename DepthType>
SampleReadSetProfileHelper<DepthType>
profile_reads(const SampleName& sample,
              const ReferenceGenome& reference,
              InputRegionMap regions,
              const ReadManager& source,
              const ReadSetProfileConfig& config,
              std::discrete_distribution<> contig_sampling_distribution,
              boost::optional<std::atomic_bool&> kill = boost::none)
{
    SampleReadSetProfileHelper<DepthType> result {};
    std::vector<DepthType> sample_depths {};
    SamplingSummary sampling_summary {};
    while (true) {
        if (kill && *kill) return result;
        const auto target_sampling_region = choose_next_sample_region(sample, regions, config, contig_sampling_distribution, sampling_summary);
        if (!target_sampling_region) break;
        CoverageTracker<GenomicRegion, DepthType> depth_tracker {true};
        auto remaining_reads = static_cast<int>(config.target_reads_per_draw);
        boost::optional<GenomicRegion> critical_region {};
        const auto read_visitor = [&] (const SampleName& sample, AlignedRead read) {
            result.read_lengths.push_back(sequence_size(read));
            result.mapping_qualities.push_back(read.mapping_quality());
            result.memory_footprints.push_back(footprint(read));
            if (config.fragment_size) {
                result.fragmented_memory_footprints.push_back(fragmented_footprint(read, *config.fragment_size));
            }
            depth_tracker.add(read);
            if (!critical_region) {
                critical_region = mapped_region(read);
                if (config.min_read_lengths > 1) {
                    critical_region = expand_rhs(*critical_region, (config.min_read_lengths - 1) * size(*critical_region));
                }
            }
            if (remaining_reads > 0) --remaining_reads;
            return remaining_reads > 0 || overlaps(read, *critical_region);
        };
        source.iterate(sample, *target_sampling_region, read_visitor);
        auto sampled_region = *target_sampling_region;
        if (depth_tracker.any()) {
            auto sampled_reads_region = *depth_tracker.encompassing_region();
            assert(!result.read_lengths.empty());
            if (remaining_reads > 0) {
                sampled_region = *target_sampling_region;
            } else {
                assert(!is_before(sampled_reads_region, *target_sampling_region));
                sampled_region = closed_region(*target_sampling_region, sampled_reads_region);
                if (size(sampled_region) > result.read_lengths.back()) {
                    // Ignore the last half read length bases to avoid adding positions undersampled because
                    // the sampled read limit was hit.
                    const auto read_length = static_cast<GenomicRegion::Distance>(result.read_lengths.back());
                    sampled_region = expand_rhs(sampled_region, -read_length / 2);
                } else {
                    sampled_region = expand_rhs(head_region(*target_sampling_region), result.read_lengths.back() / 2);
                }
            }
        }
        auto read_depths = depth_tracker.get(sampled_region);
        erase_non_dna_or_rna_positions(read_depths, sampled_region, reference);
        utils::append(read_depths, result.contig_depths[sampled_region.contig_name()]);
        utils::append(std::move(read_depths), sample_depths);
        ++sampling_summary.num_samples;
        auto removal_region = sampled_region;
        if (depth_tracker.any()) {
            removal_region = encompassing_region(removal_region, *depth_tracker.encompassing_region());
        }
        cut(removal_region, regions.at(sampled_region.contig_name()));
        if (regions.at(sampled_region.contig_name()).empty()) {
            regions.erase(sampled_region.contig_name());
            contig_sampling_distribution = make_contig_sampling_distribution(regions);
        }
    }
    if (!sample_depths.empty()) {
        std::sort(std::begin(sample_depths), std::end(sample_depths)); // sorting means no copying from stats calculations
        fill_depth_stats(sample_depths, result.depth_stats.genome);
        sample_depths.clear();
        sample_depths.shrink_to_fit();
        for (auto& p : result.contig_depths) {
            std::sort(std::begin(p.second), std::end(p.second)); // sorting means no copying from stats calculations
            result.depth_stats.contig.emplace(p.first, make_depth_stats(p.second));
        }
    }
    return result;
}

template <typename DepthType>
boost::optional<ReadSetProfile>
profile_reads_helper(const std::vector<SampleName>& samples,
                     const ReferenceGenome& reference,
                     const InputRegionMap& regions,
                     const ReadManager& source,
                     ReadSetProfileConfig config,
                     boost::optional<ThreadPool&> workers)
{
    const auto contig_sampling_distribution = make_contig_sampling_distribution(regions);
    std::vector<SampleReadSetProfileHelper<DepthType>> sample_read_set_profiles {};
    sample_read_set_profiles.reserve(samples.size());
    if (samples.size() > 1 && workers) {
        std::vector<std::future<SampleReadSetProfileHelper<DepthType>>> sample_read_set_profiles_futures {};
        sample_read_set_profiles_futures.reserve(samples.size());
        std::atomic_bool kill {false};
        for (const auto& sample : samples) {
            sample_read_set_profiles_futures.push_back(workers->push([&] () { 
                return profile_reads<DepthType>(sample, reference, regions, source, config, contig_sampling_distribution, kill); }));
        }
        std::size_t future_idx {};
        try {
            for (; future_idx < sample_read_set_profiles_futures.size(); ++future_idx) {
                sample_read_set_profiles.push_back(sample_read_set_profiles_futures[future_idx].get());
            }
        } catch (...) {
            kill = true;
            for (; future_idx < sample_read_set_profiles_futures.size(); ++future_idx) {
                try {
                    sample_read_set_profiles.push_back(sample_read_set_profiles_futures[future_idx].get());
                } catch (...) {}
            }
            throw;
        }
    } else {
        for (const auto& sample : samples) {
            sample_read_set_profiles.push_back(profile_reads<DepthType>(sample, reference, regions, source, config, contig_sampling_distribution));
        }
    }
    ReadSetProfile result {};
    std::deque<MemoryFootprint> memory_footprints {}, fragmented_memory_footprints {};
    std::unordered_map<GenomicRegion::ContigName, std::vector<DepthType>> contig_depths {};
    std::deque<unsigned> read_lengths {};
    std::deque<AlignedRead::MappingQuality> mapping_qualities {};
    for (std::size_t s {0}; s < samples.size(); ++s) {
        auto& sample_profiles = sample_read_set_profiles[s];
        utils::append(std::move(sample_profiles.memory_footprints), memory_footprints);
        utils::append(std::move(sample_profiles.fragmented_memory_footprints), fragmented_memory_footprints);
        for (auto& p : sample_profiles.contig_depths) {
            utils::append(std::move(p.second), contig_depths[p.first]);
        }
        utils::append(std::move(sample_profiles.read_lengths), read_lengths);
        utils::append(std::move(sample_profiles.mapping_qualities), mapping_qualities);
        result.depth_stats.sample.emplace(samples[s], std::move(sample_profiles.depth_stats));
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
              const ReferenceGenome& reference,
              const InputRegionMap& regions,
              const ReadManager& source,
              ReadSetProfileConfig config,
              boost::optional<ThreadPool&> workers)
{
    boost::optional<ReadSetProfile> result {};
    try {
        result = profile_reads_helper<unsigned short>(samples, reference, regions, source, config, workers);
    } catch (const std::runtime_error& e) {
        try {
            result = profile_reads_helper<unsigned>(samples, reference, regions, source, config, workers);
        } catch (const std::runtime_error& e) {
            result = profile_reads_helper<unsigned long>(samples, reference, regions, source, config, workers);
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
