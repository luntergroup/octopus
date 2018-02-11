// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "input_reads_profiler.hpp"

#include <random>
#include <deque>
#include <iterator>
#include <algorithm>
#include <utility>
#include <cassert>

#include "mappable_algorithms.hpp"
#include "maths.hpp"
#include "append.hpp"

namespace octopus {

namespace {

template <typename ForwardIt, typename RandomGenerator>
ForwardIt random_select(ForwardIt first, ForwardIt last, RandomGenerator& g)
{
    if (first == last) return first;
    const auto max = static_cast<std::size_t>(std::distance(first, last));
    std::uniform_int_distribution<std::size_t> dist {0, max - 1};
    std::advance(first, dist(g));
    return first;
}

template <typename ForwardIt>
ForwardIt random_select(ForwardIt first, ForwardIt last)
{
    static std::default_random_engine gen {};
    return random_select(first, last, gen);
}

auto estimate_dynamic_size(const AlignedRead& read) noexcept
{
    return read.name().size() * sizeof(char)
    + read.read_group().size() * sizeof(char)
    + sequence_size(read) * sizeof(char)
    + sequence_size(read) * sizeof(AlignedRead::BaseQuality)
    + read.cigar().size() * sizeof(CigarOperation)
    + contig_name(read).size() * sizeof(char)
    + (read.has_other_segment() ? sizeof(AlignedRead::Segment) : 0);
}

auto estimate_read_size(const AlignedRead& read) noexcept
{
    return sizeof(AlignedRead) + estimate_dynamic_size(read);
}

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
    static std::default_random_engine gen {};
    std::uniform_int_distribution<GenomicRegion::Position> dist {from.begin(), max_begin};
    return GenomicRegion {from.contig_name(), dist(gen), from.end()};
}

auto draw_sample(const SampleName& sample, const InputRegionMap& regions,
                 const ReadManager& source, const ReadSetProfileConfig& config)
{
    const auto contig_itr = random_select(std::cbegin(regions), std::cend(regions));
    assert(!contig_itr->second.empty());
    const auto region_itr = random_select(std::cbegin(contig_itr->second), std::cend(contig_itr->second));
    const auto sample_region = choose_sample_region(*region_itr, config.max_sample_size);
    auto test_region = source.find_covered_subregion(sample, sample_region, config.max_sample_size);
    if (is_empty(test_region)) {
        test_region = expand_rhs(test_region, 1);
    }
    return source.fetch_reads(sample, test_region);
}

using ReadSetSamples = std::vector<ReadManager::ReadContainer>;

auto draw_samples(const SampleName& sample, const InputRegionMap& regions,
                  const ReadManager& source, const ReadSetProfileConfig& config)
{
    ReadSetSamples result {};
    result.reserve(config.max_samples_per_sample);
    std::generate_n(std::back_inserter(result), config.max_samples_per_sample,
                    [&] () { return draw_sample(sample, regions, source, config); });
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
            std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(result), estimate_read_size);
        }
    }
    return result;
}

} // namespace

boost::optional<ReadSetProfile> profile_reads(const std::vector<SampleName>& samples,
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
    result.sample_mean_depth.resize(samples.size());
    result.sample_depth_stdev.resize(samples.size());
    std::deque<unsigned> depths {};
    for (std::size_t s {0}; s < samples.size(); ++s) {
        std::deque<unsigned> sample_depths {};
        for (const auto& reads : read_sets[s]) {
            utils::append(calculate_positional_coverage(reads), sample_depths);
        }
        result.sample_mean_depth[s] = maths::mean(sample_depths);
        result.sample_depth_stdev[s] = maths::stdev(sample_depths);
        utils::append(std::move(sample_depths), depths);
    }
    result.mean_depth = maths::mean(depths);
    result.depth_stdev = maths::stdev(depths);
    return result;
}

boost::optional<std::size_t> estimate_mean_read_size(const std::vector<SampleName>& samples,
                                                     const InputRegionMap& input_regions,
                                                     ReadManager& read_manager,
                                                     const unsigned max_sample_size)
{
    if (input_regions.empty()) return boost::none;
    const auto sample_regions = get_covered_sample_regions(samples, input_regions, read_manager);
    if (sample_regions.empty()) return boost::none;
    const auto num_samples_per_sample = max_sample_size / samples.size();
    std::deque<std::size_t> read_size_samples {};
    // take read samples from each sample seperatly to ensure we cover each
    for (const auto& sample : samples) {
        const auto it = random_select(std::cbegin(sample_regions), std::cend(sample_regions));
        assert(!it->second.empty());
        const auto it2 = random_select(std::cbegin(it->second), std::cend(it->second));
        auto test_region = read_manager.find_covered_subregion(sample, *it2, num_samples_per_sample);
        if (is_empty(test_region)) {
            test_region = expand_rhs(test_region, 1);
        }
        const auto reads = read_manager.fetch_reads(sample, test_region);
        std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(read_size_samples),
                       estimate_read_size);
    }
    if (read_size_samples.empty()) return boost::none;
    return static_cast<std::size_t>(maths::mean(read_size_samples) + maths::stdev(read_size_samples));
}

std::size_t default_read_size_estimate() noexcept
{
    return sizeof(AlignedRead) + 300;
}

} // namespace octopus
