// Copyright (c) 2015-2021 Daniel Cooke
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

BadRegionDetector::BadRegionDetector(Parameters params, boost::optional<const ReadSetProfile&> reads_profile)
: params_ {params}
, reads_profile_ {reads_profile}
{}

std::vector<BadRegionDetector::BadRegion>
BadRegionDetector::detect(const MappableFlatSet<Variant>& variants,
                          const ReadMap& reads,
                          OptionalReadsReport reads_report) const
{
    return detect(InputData {reads, variants, reads_report}, reads_report);
}

std::vector<BadRegionDetector::BadRegion>
BadRegionDetector::detect(const ReadMap& reads, OptionalReadsReport reads_report) const
{
    return detect(InputData {reads, boost::none, reads_report}, reads_report);
}

// private methods

std::vector<BadRegionDetector::BadRegion>
BadRegionDetector::detect(const InputData& data, OptionalReadsReport reads_report) const
{
    auto candidate_bad_regions = get_candidate_bad_regions(data);
    if (candidate_bad_regions.empty()) return {};
    const auto candidate_bad_regions_states = compute_states(candidate_bad_regions, data);
    const auto is_bad_state = [&] (const auto& state) { return is_bad(state, reads_report); };
    auto bad_state_itr = std::find_if(std::cbegin(candidate_bad_regions_states), std::cend(candidate_bad_regions_states), is_bad_state);
    std::vector<BadRegion> result {};
    result.reserve(result.size());
    while (bad_state_itr != std::cend(candidate_bad_regions_states)) {
        const auto state_idx = static_cast<std::size_t>(std::distance(std::cbegin(candidate_bad_regions_states), bad_state_itr));
        auto next_bad_state_itr = std::find_if(std::next(bad_state_itr), std::cend(candidate_bad_regions_states), is_bad_state);
        if (next_bad_state_itr != std::cend(candidate_bad_regions_states) && next_bad_state_itr == std::next(bad_state_itr)) {
            auto connecting_region = *intervening_region(candidate_bad_regions[state_idx], candidate_bad_regions[state_idx + 1]);
            auto connecting_state = compute_state(connecting_region, data);
            if (is_bad(connecting_state, reads_report)) {
                result.push_back({std::move(connecting_region), BadRegion::Severity::low});
            }
        }
        result.push_back({std::move(candidate_bad_regions[state_idx]), BadRegion::Severity::high});
        bad_state_itr = next_bad_state_itr;
    }
    return result;
}

namespace {

auto find_high_depth_regions_helper(const GenomicRegion& target_region, const ReadPipe::Report::DepthMap& read_depths, const ReadSetProfile& profile)
{
    std::vector<GenomicRegion> sample_high_depth_regions {};
    const auto num_samples = read_depths.size();
    sample_high_depth_regions.reserve(num_samples);
    for (const auto& p : read_depths) {
        const auto sample_depths = p.second.get(target_region);
        auto depth_threshold = profile.depth_stats.combined.genome.positive.median;
        if (profile.depth_stats.combined.contig.count(target_region.contig_name()) == 1) {
            depth_threshold = std::max(depth_threshold, profile.depth_stats.combined.contig.at(target_region.contig_name()).positive.median);
        }
        if (profile.depth_stats.sample.count(p.first) == 1) {
            depth_threshold = std::max(depth_threshold, profile.depth_stats.sample.at(p.first).genome.positive.median);
            if (profile.depth_stats.sample.at(p.first).contig.count(target_region.contig_name()) == 1) {
                depth_threshold = std::max(depth_threshold, profile.depth_stats.sample.at(p.first).contig.at(target_region.contig_name()).positive.median);
            }
        } else if (profile.depth_stats.combined.contig.count(target_region.contig_name()) == 1) {
            depth_threshold = std::max(depth_threshold, profile.depth_stats.combined.contig.at(target_region.contig_name()).positive.median);
        }
		depth_threshold *= 4;
        utils::append(find_high_coverage_regions(sample_depths, target_region, static_cast<unsigned>(depth_threshold)), sample_high_depth_regions);
    }
    if (num_samples == 1) {
        return sample_high_depth_regions;
    } else {
        std::sort(std::begin(sample_high_depth_regions), std::end(sample_high_depth_regions));
        return extract_covered_regions(sample_high_depth_regions);
    }
}

auto find_high_depth_regions_helper(const GenomicRegion& target_region, const ReadMap& reads, const ReadSetProfile& profile)
{
    ReadPipe::Report::DepthMap depths {};
    depths.reserve(reads.size());
    for (const auto& p : reads) depths.emplace(p.first, make_coverage_tracker(p.second));
    return find_high_depth_regions_helper(target_region, depths, profile);
}

} // namespace

std::vector<GenomicRegion>
BadRegionDetector::find_high_depth_regions(const ReadMap& reads, OptionalReadsReport reads_reports) const
{
    if (reads_profile_) {
        const auto target_region = encompassing_region(reads);
        if (reads_reports) {
            return find_high_depth_regions_helper(target_region, reads_reports->raw_depths, *reads_profile_);
        } else {
            return find_high_depth_regions_helper(target_region, reads, *reads_profile_);
        }
    } else {
        return {};
    }
}

namespace {

void merge(std::vector<GenomicRegion> src, std::vector<GenomicRegion>& dst)
{
    auto itr = utils::append(std::move(src), dst);
    std::inplace_merge(std::begin(dst), itr, std::end(dst));
}

} // namespace

std::vector<GenomicRegion>
BadRegionDetector::get_candidate_bad_regions(const InputData& data) const
{
    std::vector<GenomicRegion> result {};
    if (data.variants) {
        result = get_candidate_dense_regions(*data.variants, data.reads, data.reads_report);
    }
    auto high_depth_regions = find_high_depth_regions(data.reads, data.reads_report);
    if (!high_depth_regions.empty()) {
        merge(std::move(high_depth_regions), result);
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

std::vector<GenomicRegion>
find_dense_regions(const MappableFlatSet<Variant>& variants, const ReadMap& reads,
                   const double dense_zone_log_count_threshold,
                   const double max_shared_dense_zones)
{
    if (variants.empty()) return {};
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
                                               OptionalReadsReport reads_report) const
{
    if (candidates.empty()) return {};
    const auto average_read_length = mean_read_length(reads);
    auto expected_log_count = get_max_expected_log_allele_count_per_base();
    const auto dense_zone_log_count_threshold = expected_log_count * average_read_length;
    return find_dense_regions(candidates, reads, dense_zone_log_count_threshold, 1);
}

double BadRegionDetector::get_max_expected_log_allele_count_per_base() const noexcept
{
    using Tolerance = Parameters::Tolerance;
    double stdev_multiplier {1};
    switch (params_.tolerance) {
        case Tolerance::low: stdev_multiplier = 3; break;
        case Tolerance::high: stdev_multiplier = 7; break;
        case Tolerance::normal:
        default: stdev_multiplier = 6;
    }
    return params_.heterozygosity + stdev_multiplier * params_.heterozygosity_stdev;
}

void BadRegionDetector::fill(RegionState::ReadSummaryStats& stats,
                             const ReadMap& reads,
                             const GenomicRegion& region,
                             OptionalReadsReport reads_report) const
{
    if (has_coverage(reads, region)) {
        stats.max_length = max_read_length(reads);
        stats.median_mapping_quality = median_mapping_quality(reads, region);
    } else {
        if (reads_profile_) {
            stats.median_mapping_quality = reads_profile_->mapping_quality_stats.max;
        } else {
            stats.median_mapping_quality = 60;
        }
    }
    stats.average_depths.reserve(reads.size());
    if (reads_report) {
        for (const auto& p : reads_report->raw_depths) {
            stats.average_depths.emplace(p.first, p.second.mean(region));
        }
    } else {
        for (const auto& p : reads) {
            stats.average_depths.emplace(p.first, mean_coverage(p.second, region));
        }
    }
}

void BadRegionDetector::fill(RegionState::VariantSummaryStats& stats,
                             const MappableFlatSet<Variant>& variants,
                             const GenomicRegion& region) const
{
    stats.count = count_contained(variants, region);
    stats.density = static_cast<double>(stats.count) / size(region);
}

BadRegionDetector::RegionState
BadRegionDetector::compute_state(const GenomicRegion& region, const InputData& data) const
{
    RegionState result {};
    result.region = region;
    fill(result.read_stats, data.reads, region, data.reads_report);
    if (data.variants) {
        result.variant_stats = RegionState::VariantSummaryStats {};
        fill(*result.variant_stats, *data.variants, region);
    }
    return result;
}

std::vector<BadRegionDetector::RegionState>
BadRegionDetector::compute_states(const std::vector<GenomicRegion>& regions, const InputData& data) const
{
    std::vector<RegionState> result {};
    result.reserve(regions.size());
    for (const auto& region : regions) {
        result.push_back(compute_state(region, data));
    }
    return result;
}

namespace {

template <typename PairRange>
auto sum_second(const PairRange& pairs)
{
    using T = typename PairRange::value_type::second_type;
    const static auto add_second = [] (T total, const auto& p) { return total + p.second; };
    return std::accumulate(std::cbegin(pairs), std::cend(pairs), T {}, add_second);
}

bool are_similar(const ReadSetProfile::DepthSummaryStats& lhs, const ReadSetProfile::DepthSummaryStats& rhs)
{
	const auto diff = lhs.median < rhs.median ? (rhs.median - rhs.median) : (lhs.median < rhs.median);
	return diff <= std::min(lhs.stdev, rhs.stdev);
}

template <typename Range>
auto discrete_cdf(const Range& probabilities, const std::size_t x)
{
    using ProbabilityType = typename Range::value_type;
    if (x >= probabilities.size()) return ProbabilityType {1};
    return std::accumulate(std::cbegin(probabilities), std::next(std::cbegin(probabilities), x), ProbabilityType {0});
}

template <typename Range>
auto discrete_sf(const Range& probabilities, std::size_t x)
{
    return 1 - discrete_cdf(probabilities, x);
}

auto calculate_conditional_depth_probability(const unsigned target_depth, const ReadSetProfile::DepthStats& stats)
{
    const auto low_depth = std::max(stats.positive.mean, stats.positive.median) + (stats.positive.stdev / (2 - stats.distribution[0]));
    double result {1};
    if (low_depth < target_depth) {
        result *= discrete_sf(stats.distribution, target_depth);
        result /= discrete_sf(stats.distribution, low_depth);
    }
    return result;
}

} // namespace

double BadRegionDetector::calculate_probability_good(const RegionState& state, OptionalReadsReport reads_report) const
{
    // Basic idea
    // lower mapping quality -> lower probability
    // higher variant density -> lower probability
    // higher depth -> lower probability
    double result {1};
    if (reads_profile_) {
        const auto& depth_stats = reads_profile_->depth_stats;
        const auto average_depth = sum_second(state.read_stats.average_depths) / state.read_stats.average_depths.size();
        if (depth_stats.combined.contig.count(state.region.contig_name()) == 1) {
            result *= calculate_conditional_depth_probability(average_depth, depth_stats.combined.contig.at(state.region.contig_name()));
        } else {
            result *= calculate_conditional_depth_probability(average_depth, depth_stats.combined.genome);
        }
        if (state.read_stats.median_mapping_quality < reads_profile_->mapping_quality_stats.median) {
            result /= std::max(std::min((reads_profile_->mapping_quality_stats.median - state.read_stats.median_mapping_quality) / 10, 4), 1);
        }
    } else if (state.read_stats.median_mapping_quality < 40) {
        result /= 2;
    }
    if (reads_report) {
        unsigned mapping_quality_zero_depth {0}, total_depth {0};
        for (const auto& p : reads_report->mapping_quality_zero_depths) {
            mapping_quality_zero_depth += p.second.sum(state.region);
            total_depth += reads_report->raw_depths.at(p.first).sum(state.region);
        }
        const auto average_mapping_quality_zero_fraction = static_cast<double>(mapping_quality_zero_depth) / total_depth;
        result *= std::max(1 - average_mapping_quality_zero_fraction, 0.25);
    }
    if (state.variant_stats) {
        double tolerance_factor;
        switch (params_.tolerance) {
            case Parameters::Tolerance::high: tolerance_factor = 50; break;
            case Parameters::Tolerance::normal: tolerance_factor = 40; break;
            case Parameters::Tolerance::low: tolerance_factor = 30; break;
        }
        const auto density_mean = size(state.region) * (params_.heterozygosity + tolerance_factor * params_.heterozygosity_stdev);
        result *= maths::poisson_sf(state.variant_stats->count, density_mean);
    }
    if (size(state.region) > 1'000) {
        result = std::pow(result, static_cast<double>(size(state.region)) / 1'000); // penalise large regions
    }
    return result;
}

bool BadRegionDetector::is_bad(const RegionState& state, OptionalReadsReport reads_report) const
{
    GenomicRegion::Size min_region_size {100};
    unsigned min_alleles {0};
    double min_probability_good {0.5};
    switch (params_.tolerance) {
        case Parameters::Tolerance::high: {
            min_alleles = 20;
            min_probability_good = 0.0001;
            min_region_size = 300;
            break;
        }
        case Parameters::Tolerance::normal: {
            min_alleles = 10;
            min_probability_good = 0.005;
			min_region_size = 200;
            break;
        }
        case Parameters::Tolerance::low: {
            min_alleles = 5;
            min_probability_good = 0.01;
            min_region_size = 100;
            break;
        }
    }
    return (!state.variant_stats || state.variant_stats->count >= min_alleles)
        && size(state.region) > min_region_size
        && calculate_probability_good(state, reads_report) < min_probability_good;
}

} // namespace coretools
} // namespace octopus
