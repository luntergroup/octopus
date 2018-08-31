// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "dense_variation_detector.hpp"

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

#include <iostream>

namespace octopus { namespace coretools {

DenseVariationDetector::DenseVariationDetector(double heterozygosity, double heterozygosity_stdev,
                                               boost::optional<ReadSetProfile> reads_profile)
: expected_heterozygosity_ {heterozygosity}
, heterozygosity_stdev_ {heterozygosity_stdev}
, reads_profile_ {std::move(reads_profile)}
{}

namespace {

auto get_max_expected_log_allele_count_per_base(double heterozygosity, double heterozygosity_stdev) noexcept
{
    return heterozygosity + 6 * heterozygosity_stdev;
}

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

bool all_empty(const ReadMap& reads) noexcept
{
    return std::all_of(std::cbegin(reads), std::cend(reads), [](const auto& p) noexcept { return p.second.empty(); });
}

auto mean_mapped_region_size(const ReadMap& reads) noexcept
{
    double total_read_size {0};
    std::size_t num_reads {0};
    for (const auto& p : reads) {
        total_read_size += sum_region_sizes(p.second);
        num_reads += p.second.size();
    }
    return total_read_size / num_reads;
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
    MappableFlatSet<GenomicRegion> result {};
    const auto join_threshold = max_shared_dense_zones * dense_zone_log_count_threshold;
    for (const auto& block : dense_blocks) {
        if (block.log_count > join_threshold) {
            result.insert(GenomicRegion {contig, block.region});
        }
    }
    return result;
}

struct RegionState
{
    GenomicRegion region;
    unsigned variant_count;
    double variant_density;
    unsigned mean_read_depth;
    double rmq_mapping_quality;
};

auto compute_state(const GenomicRegion& region, const MappableFlatSet<Variant>& variants, const ReadMap& reads)
{
    RegionState result {};
    result.region = region;
    result.rmq_mapping_quality = rmq_mapping_quality(reads, region);
    result.mean_read_depth = mean_coverage(reads, region) / reads.size();
    result.variant_count = count_overlapped(variants, region);
    result.variant_density = static_cast<double>(result.variant_count) / size(region);
    return result;
}

auto compute_states(const MappableFlatSet<GenomicRegion>& regions, const MappableFlatSet<Variant>& variants, const ReadMap& reads)
{
    std::vector<RegionState> result {};
    result.reserve(regions.size());
    for (const auto& region : regions) {
        result.push_back(compute_state(region, variants, reads));
    }
    return result;
}

bool should_join(const RegionState& lhs_state, const RegionState& connecting_state, const RegionState& rhs_state)
{
    if (connecting_state.variant_density > std::min(lhs_state.variant_density, rhs_state.variant_density) / 2) {
        return true;
    }
    if (size(connecting_state.region) > std::max(size(lhs_state.region), size(rhs_state.region))) {
        return false;
    }
    if (connecting_state.rmq_mapping_quality < std::min(lhs_state.rmq_mapping_quality, rhs_state.rmq_mapping_quality)) {
        return true;
    }
    if (connecting_state.mean_read_depth > std::min(lhs_state.mean_read_depth, rhs_state.mean_read_depth)) {
        return true;
    }
    return false;
}

auto join_dense_regions(const MappableFlatSet<GenomicRegion>& dense_regions,
                        const MappableFlatSet<Variant>& variants, const ReadMap& reads)
{
    if (dense_regions.size() > 1) {
        std::vector<GenomicRegion> final_regions {};
        final_regions.reserve(dense_regions.size());
        const auto dense_states = compute_states(dense_regions, variants, reads);
        final_regions.push_back(dense_regions.front());
        for (std::size_t i {1}; i < dense_regions.size(); ++i) {
            const auto connecting_region = *intervening_region(dense_regions[i - 1], dense_regions[i]);
            const auto connecting_state = compute_state(connecting_region, variants, reads);
            if (should_join(dense_states[i - 1], connecting_state, dense_states[i])) {
                final_regions.push_back(connecting_region);
            }
            final_regions.push_back(dense_regions[i]);
        }
        return extract_covered_regions(final_regions);
    } else {
        return std::vector<GenomicRegion> {dense_regions.front()};
    }
}

} // namespace

std::vector<DenseVariationDetector::DenseRegion>
DenseVariationDetector::detect(const MappableFlatSet<Variant>& variants, const ReadMap& reads,
                               boost::optional<const ReadPipe::Report&> reads_report) const
{
    const auto mean_read_size = mean_mapped_region_size(reads);
    auto expected_log_count = get_max_expected_log_allele_count_per_base(expected_heterozygosity_, heterozygosity_stdev_);
    const auto dense_zone_log_count_threshold = expected_log_count * mean_read_size;
    auto dense_regions = find_dense_regions(variants, reads, dense_zone_log_count_threshold, 1);
    if (dense_regions.empty()) return {};
    auto joined_dense_regions = join_dense_regions(dense_regions, variants, reads);
    std::vector<DenseRegion> result {};
    result.reserve(joined_dense_regions.size());
    double max_expected_coverage {};
    if (reads_profile_) {
        max_expected_coverage = 2 * std::max(reads_profile_->mean_depth, reads_profile_->median_depth) + 2 * reads_profile_->depth_stdev;
    } else {
        max_expected_coverage = 2 * mean_coverage(reads);
    }
    for (const auto& region : joined_dense_regions) {
        const auto state = compute_state(region, variants, reads);
        auto total_mean_depth = state.mean_read_depth;
        if (reads_report) {
            const auto num_downsampled_reads = count_downsampled_reads(reads_report->downsample_report, region);
            const auto mean_downsample_depth = static_cast<double>(num_downsampled_reads) / size(state.region);
            total_mean_depth += mean_downsample_depth;
        }
        if (state.variant_count > 100 && size(state.region) > mean_read_size && total_mean_depth > max_expected_coverage) {
            result.push_back({region, DenseRegion::RecommendedAction::skip});
        }
    }
    return result;
}

} // namespace coretools
} // namespace octopus
