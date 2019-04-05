// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "bad_region_detector.hpp"

#include <deque>
#include <array>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <limits>
#include <cmath>

#include "concepts/mappable.hpp"
#include "concepts/comparable.hpp"
#include "basics/contig_region.hpp"
#include "basics/aligned_read.hpp"
#include "containers/mappable_flat_set.hpp"
#include "containers/mappable_map.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/maths.hpp"
#include "utils/append.hpp"

#include <iostream>

namespace octopus { namespace coretools {

BadRegionDetector::BadRegionDetector(Parameters params, boost::optional<ReadSetProfile> reads_profile)
: params_ {params}
, reads_profile_ {std::move(reads_profile)}
{}

std::vector<BadRegionDetector::BadRegion>
BadRegionDetector::detect(const MappableFlatSet<Variant>& candidate_variants, const ReadMap& reads,
                               boost::optional<const ReadPipe::Report&> reads_report) const
{
    auto candidate_bad_regions = get_candidate_bad_regions(candidate_variants, reads, reads_report);
    if (candidate_bad_regions.empty()) return {};
    const auto candidate_bad_regions_states = compute_states(candidate_bad_regions, candidate_variants, reads, reads_report);
    const auto is_bad_state = [this] (const auto& state) { return is_bad(state); };
    auto bad_state_itr = std::find_if(std::cbegin(candidate_bad_regions_states), std::cend(candidate_bad_regions_states), is_bad_state);
    std::vector<BadRegion> result {};
    result.reserve(result.size());
    while (bad_state_itr != std::cend(candidate_bad_regions_states)) {
        const auto state_idx = static_cast<std::size_t>(std::distance(std::cbegin(candidate_bad_regions_states), bad_state_itr));
        auto next_bad_state_itr = std::find_if(std::next(bad_state_itr), std::cend(candidate_bad_regions_states), is_bad_state);
        if (next_bad_state_itr != std::cend(candidate_bad_regions_states) && next_bad_state_itr == std::next(bad_state_itr)) {
            auto connecting_region = *intervening_region(candidate_bad_regions[state_idx], candidate_bad_regions[state_idx + 1]);
            auto connecting_state = compute_state(connecting_region, candidate_variants, reads, reads_report);
            if (is_bad(connecting_state)) {
                result.push_back({std::move(connecting_region), BadRegion::RecommendedAction::restrict_lagging});
            }
        }
        result.push_back({std::move(candidate_bad_regions[state_idx]), BadRegion::RecommendedAction::skip});
        bad_state_itr = next_bad_state_itr;
    }
    return result;
}

namespace {

template <typename Range, typename T>
std::size_t find_index(const Range& values, const T& value)
{
    return std::distance(std::cbegin(values), std::find(std::cbegin(values), std::cend(values), value));
}

auto find_high_depth_regions(const GenomicRegion& target_region, const ReadPipe::Report::DepthMap& read_depths, const ReadSetProfile& profile)
{
    std::vector<GenomicRegion> sample_high_depth_regions {};
    const auto num_samples = read_depths.size();
    assert(read_depths.size() <= profile.samples.size());
    sample_high_depth_regions.reserve(num_samples);
    for (const auto& p : read_depths) {
        const auto sample_depths = p.second.get(target_region);
        const auto sample_idx = find_index(profile.samples, p.first);
        assert(sample_idx < profile.sample_median_positive_depth.size());
        const auto depth_threshold = static_cast<unsigned>(5 * profile.sample_median_positive_depth[sample_idx]);
        utils::append(find_high_coverage_regions(sample_depths, target_region, depth_threshold), sample_high_depth_regions);
    }
    if (num_samples == 1) {
        return sample_high_depth_regions;
    } else {
        std::sort(std::begin(sample_high_depth_regions), std::end(sample_high_depth_regions));
        return extract_covered_regions(sample_high_depth_regions);
    }
}

auto find_high_depth_regions(const GenomicRegion& target_region, const ReadMap& reads, const ReadSetProfile& profile)
{
    ReadPipe::Report::DepthMap depths {};
    depths.reserve(reads.size());
    for (const auto& p : reads) depths.emplace(p.first, make_coverage_tracker(p.second));
    return find_high_depth_regions(target_region, depths, profile);
}

auto find_high_depth_regions(const ReadMap& reads, const ReadSetProfile& profile,
                             boost::optional<const ReadPipe::Report&> reads_reports)
{
    const auto target_region = encompassing_region(reads);
    if (reads_reports) {
        return find_high_depth_regions(target_region, reads_reports->raw_depths, profile);
    } else {
        return find_high_depth_regions(target_region, reads, profile);
    }
}

void merge(std::vector<GenomicRegion> src, std::vector<GenomicRegion>& dst)
{
    auto itr = utils::append(std::move(src), dst);
    std::inplace_merge(std::begin(dst), itr, std::end(dst));
}


} // namespace

std::vector<GenomicRegion>
BadRegionDetector::get_candidate_bad_regions(const MappableFlatSet<Variant>& candidate_variants, const ReadMap& reads,
                                             boost::optional<const ReadPipe::Report&> reads_report) const
{
    auto result = get_candidate_dense_regions(candidate_variants, reads, reads_report);
    if (reads_profile_) {
        merge(find_high_depth_regions(reads, *reads_profile_, reads_report), result);
        result = extract_covered_regions(result);
    }
    return result;
}

namespace {

struct AlleleBlock : public Mappable<AlleleBlock>, public Comparable<AlleleBlock>
{
    ContigRegion region;
    double log_count;
    const ContigRegion& mapped_region() const noexcept { return region; }
    AlleleBlock() = default;
    AlleleBlock(ContigRegion region, double count) noexcept : region {region}, log_count {count} {}
    AlleleBlock(const GenomicRegion& region, double count) noexcept : region {region.contig_region()}, log_count {count} {}
};

bool operator==(const AlleleBlock& lhs, const AlleleBlock& rhs) noexcept { return lhs.region < rhs.region; }
bool operator<(const AlleleBlock& lhs, const AlleleBlock& rhs) noexcept { return lhs.region < rhs.region; }

template <typename Range>
auto sum_block_counts(const Range& block_range) noexcept
{
    return std::accumulate(std::cbegin(block_range), std::cend(block_range), 0.0,
                           [] (const auto curr, const AlleleBlock& block) noexcept { return curr + block.log_count; });
}

auto calculate_positional_coverage(const MappableFlatSet<Variant>& variants, const GenomicRegion& region)
{
    std::vector<ContigRegion> variant_regions(variants.size());
    std::transform(std::cbegin(variants), std::cend(variants), std::begin(variant_regions),
                   [] (const auto& v) { return is_simple_insertion(v) ? expand(contig_region(v), 1) : contig_region(v); });
    return calculate_positional_coverage(variant_regions, region.contig_region());
}

struct BaseState
{
    std::uint8_t variant_depth;
    std::uint16_t read_depth;
    AlignedRead::MappingQuality median_mq;
};

auto compute_median(std::vector<AlignedRead::MappingQuality>& mapping_qualities)
{
    if (mapping_qualities.empty()) {
        static constexpr auto max_mapping_quality = std::numeric_limits<AlignedRead::MappingQuality>::max();
        return max_mapping_quality;
    } else {
        return maths::median<AlignedRead::MappingQuality>(mapping_qualities);
    }
}

auto compute_positional_median_mapping_qualities(const MappableFlatMultiSet<AlignedRead>& reads, const GenomicRegion& region)
{
    const auto num_positions = size(region);
    std::vector<AlignedRead::MappingQuality> result(num_positions);
    std::deque<std::vector<AlignedRead::MappingQuality>> pileup_buffer {};
    auto pileup_region = head_region(region.contig_region());
    auto result_itr = std::begin(result);
    auto pileup_buffer_begin_itr = std::begin(pileup_buffer);
    for (const AlignedRead& read : overlap_range(reads, region)) {
        const auto& read_region = contig_region(read);
        if (ends_before(pileup_region, read_region)) {
            const auto extension_size = right_overhang_size(read_region, pileup_region);
            pileup_buffer.resize(pileup_buffer.size() + extension_size);
            pileup_buffer_begin_itr = std::begin(pileup_buffer);
            pileup_region = expand_rhs(pileup_region, extension_size);
        }
        auto pileup_read_overlap_begin_itr = pileup_buffer_begin_itr;
        if (begins_before(pileup_region, read_region)) {
            const auto offset = left_overhang_size(pileup_region, read_region);
            pileup_read_overlap_begin_itr = std::next(pileup_buffer_begin_itr, offset);
        }
        const auto num_overlapped_positions = overlap_size(pileup_region, read_region);
        const auto pileup_read_overlap_end_itr = std::next(pileup_read_overlap_begin_itr, num_overlapped_positions);
        for (auto itr = pileup_read_overlap_begin_itr; itr != pileup_read_overlap_end_itr; ++itr) {
            itr->push_back(read.mapping_quality());
        }
        if (pileup_read_overlap_begin_itr != pileup_buffer_begin_itr) {
            result_itr = std::transform(pileup_buffer_begin_itr, pileup_read_overlap_begin_itr, result_itr, compute_median);
            pileup_buffer_begin_itr = pileup_buffer.erase(pileup_buffer_begin_itr, pileup_read_overlap_begin_itr);
            pileup_region = expand_lhs(pileup_region, begin_distance(read_region, pileup_region));
        }
    }
    return result;
}

auto compute_positional_mean_median_mapping_qualities(const ReadMap& reads, const GenomicRegion& region)
{
    if (reads.size() == 1) {
        return compute_positional_median_mapping_qualities(std::cbegin(reads)->second, region);
    }
    const auto num_positions = size(region);
    std::vector<float> mean_mqs(num_positions);
    for (const auto& p : reads) {
        const auto median_mqs = compute_positional_median_mapping_qualities(p.second, region);
        for (std::size_t i {0}; i < num_positions; ++i) {
            mean_mqs[i] += (static_cast<float>(median_mqs[i]) - mean_mqs[i]) / num_positions;
        }
    }
    return std::vector<AlignedRead::MappingQuality> {std::cbegin(mean_mqs), std::cend(mean_mqs)};
}

auto compute_base_states(const MappableFlatSet<Variant>& variants, const ReadMap& reads)
{
    const auto region = encompassing_region(variants);
    const auto num_bases = size(region);
    std::vector<BaseState> result(num_bases);
    {
        auto variant_depths = calculate_positional_coverage(variants, region);
        for (std::size_t i {0}; i < num_bases; ++i) {
            static constexpr unsigned max_variant_depth {std::numeric_limits<std::uint8_t>::max()};
            result[i].variant_depth = std::min(variant_depths[i], max_variant_depth);
        }
    }
    {
        auto read_depths = calculate_positional_coverage(reads, region);
        for (std::size_t i {0}; i < num_bases; ++i) {
            static constexpr unsigned max_read_depth {std::numeric_limits<std::uint16_t>::max()};
            result[i].read_depth = std::min(read_depths[i], max_read_depth);
        }
    }
    {
        auto median_mqs = compute_positional_mean_median_mapping_qualities(reads, region);
        for (std::size_t i {0}; i < num_bases; ++i) {
            result[i].median_mq = median_mqs[i];
        }
    }
    return result;
}

auto find_dense_regions(const MappableFlatSet<Variant>& variants, const ReadMap& reads,
                        const double dense_zone_log_count_threshold,
                        const double max_shared_dense_zones)
{
    const auto initial_blocks = extract_covered_regions(variants);
    MappableFlatSet<AlleleBlock> blocks {};
    for (const auto& region : initial_blocks) {
        blocks.emplace(region, std::log2(2 * count_overlapped(variants, region)));
    }
    for (const auto& p : reads) {
        for (const auto& read : p.second) {
            const auto interacting_blocks = bases(contained_range(blocks, contig_region(read)));
            if (size(interacting_blocks) > 1) {
                auto joined_block_region = closed_region(interacting_blocks.front(), interacting_blocks.back());
                auto joined_block_count = sum_block_counts(interacting_blocks);
                auto hint = blocks.erase(std::cbegin(interacting_blocks), std::cend(interacting_blocks));
                blocks.insert(hint, AlleleBlock {joined_block_region, joined_block_count});
            }
        }
    }
    std::deque<GenomicRegion> seeds {};
    const auto& contig = contig_name(variants.front());
    for (const auto& block : blocks) {
        if (block.log_count > dense_zone_log_count_threshold) {
            seeds.push_back(GenomicRegion {contig, block.region});
        }
    }
    auto joined_seeds = join_if(seeds, [&] (const auto& lhs, const auto& rhs) { return has_shared(reads, lhs, rhs); });
    std::vector<AlleleBlock> dense_blocks(joined_seeds.size());
    std::transform(std::cbegin(joined_seeds), std::cend(joined_seeds), std::begin(dense_blocks),
                   [&](const auto& region) {
                       auto contained = contained_range(blocks, region.contig_region());
                       assert(!empty(contained));
                       auto joined_block_region = closed_region(contained.front(), contained.back());
                       auto joined_block_count = sum_block_counts(contained);
                       return AlleleBlock {joined_block_region, joined_block_count};
                   });
    std::vector<GenomicRegion> result {};
    const auto join_threshold = max_shared_dense_zones * dense_zone_log_count_threshold;
    for (const auto& block : dense_blocks) {
        if (block.log_count > join_threshold) {
            result.push_back(GenomicRegion {contig, block.region});
        }
    }
    return result;
}

} // namespace

std::vector<GenomicRegion>
BadRegionDetector::get_candidate_dense_regions(const MappableFlatSet<Variant>& candidates, const ReadMap& reads,
                                               boost::optional<const ReadPipe::Report&> reads_report) const
{
    const auto average_read_length = mean_read_length(reads);
    auto expected_log_count = get_max_expected_log_allele_count_per_base();
    const auto dense_zone_log_count_threshold = expected_log_count * average_read_length;
    return find_dense_regions(candidates, reads, dense_zone_log_count_threshold, 1);
}

double BadRegionDetector::get_max_expected_log_allele_count_per_base() const noexcept
{
    using Tolerance = Parameters::Tolerance;
    double stdev_multiplier {1};
    switch (params_.density_tolerance) {
        case Tolerance::low: stdev_multiplier = 3; break;
        case Tolerance::high: stdev_multiplier = 7; break;
        case Tolerance::normal:
        default: stdev_multiplier = 6;
    }
    return params_.heterozygosity + stdev_multiplier * params_.heterozygosity_stdev;
}

BadRegionDetector::RegionState
BadRegionDetector::compute_state(const GenomicRegion& region, const MappableFlatSet<Variant>& variants,
                                 const ReadMap& reads, boost::optional<const ReadPipe::Report&> reads_report) const
{
    RegionState result {};
    result.region = region;
    result.median_mapping_quality = median_mapping_quality(reads, region);
    result.sample_mean_read_depths.reserve(reads.size());
    result.sample_mean_read_depths.reserve(reads.size());
    if (reads_report) {
        for (const auto& p : reads_report->raw_depths) {
            result.sample_mean_read_depths.emplace(p.first, p.second.mean(region));
        }
    } else {
        for (const auto& p : reads) {
            result.sample_mean_read_depths.emplace(p.first, mean_coverage(p.second, region));
        }
    }
    result.variant_count = count_contained(variants, region);
    result.variant_density = static_cast<double>(result.variant_count) / size(region);
    result.max_read_length = max_read_length(reads);
    return result;
}

std::vector<BadRegionDetector::RegionState>
BadRegionDetector::compute_states(const std::vector<GenomicRegion>& regions, const MappableFlatSet<Variant>& variants,
                                  const ReadMap& reads, boost::optional<const ReadPipe::Report&> reads_report) const
{
    std::vector<RegionState> result {};
    result.reserve(regions.size());
    for (const auto& region : regions) {
        result.push_back(compute_state(region, variants, reads, reads_report));
    }
    return result;
}

double BadRegionDetector::calculate_probability_good(const RegionState& state) const
{
    // Basic idea
    // larger regions -> lower probability
    // lower mapping quality -> higher probability
    // higher variant density -> higher probability
    // higher depth -> higher probability
    double result {1};
    if (reads_profile_) {
        for (std::size_t s {0}; s < reads_profile_->samples.size(); ++s) {
            auto mu = 5 * std::max(reads_profile_->sample_median_positive_depth[s], reads_profile_->sample_mean_positive_depth[s]);
            auto sigma = reads_profile_->sample_depth_stdev[s];
            result *= maths::normal_sf<double>(state.sample_mean_read_depths.at(reads_profile_->samples[s]), mu, sigma);
        }
    }
    const auto density_mean = size(state.region) * (params_.heterozygosity + 2 * params_.heterozygosity_stdev);
    result *= maths::poisson_sf(state.variant_count, density_mean);
    if (state.median_mapping_quality < 40) {
        result /= 2;
    }
    return result;
}

bool BadRegionDetector::is_bad(const RegionState& state) const
{
    const auto lower_bound = 0.01 * std::pow(0.9, state.sample_mean_read_depths.size() - 1);
    return state.variant_count > 0 && calculate_probability_good(state) < lower_bound;
}

} // namespace coretools
} // namespace octopus
