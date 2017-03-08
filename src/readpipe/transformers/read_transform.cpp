// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_transform.hpp"

#include <algorithm>
#include <iterator>
#include <tuple>
#include <cassert>

#include "utils/maths.hpp"

namespace octopus { namespace readpipe {

void CapitaliseBases::operator()(AlignedRead& read) const noexcept
{
    capitalise_bases(read);
}

CapBaseQualities::CapBaseQualities(BaseQuality max) : max_ {max} {}

void CapBaseQualities::operator()(AlignedRead& read) const noexcept
{
    cap_qualities(read, max_);
}

void MaskOverlappedSegment::operator()(AlignedRead& read) const noexcept
{
    // Only reads in the forward direction are masked to prevent double masking
    if (read.has_other_segment() && contig_name(read) == read.next_segment().contig_name()
        && !read.next_segment().is_marked_unmapped() && !read.is_marked_reverse_mapped()) {
        const auto next_segment_begin = read.next_segment().begin();
        if (next_segment_begin < mapped_end(read)) {
            const auto overlapped_size = mapped_end(read) - next_segment_begin;
            zero_back_qualities(read, overlapped_size);
        }
    }
}

void MaskAdapters::operator()(AlignedRead& read) const noexcept
{
    if (read.has_other_segment() && read.is_marked_all_segments_in_read_aligned()
        && contig_name(read) == read.next_segment().contig_name()) {
        const auto insert_size = read.next_segment().inferred_template_length();
        const auto read_size   = sequence_size(read);
        if (insert_size < read_size) {
            const auto num_adapter_bases = read_size - insert_size;
            if (read.is_marked_reverse_mapped()) {
                zero_front_qualities(read, num_adapter_bases);
            } else {
                zero_back_qualities(read, num_adapter_bases);
            }
        }
    }
}

MaskTail::MaskTail(Length num_bases) : num_bases_ {num_bases} {}

void MaskTail::operator()(AlignedRead& read) const noexcept
{
    if (read.is_marked_reverse_mapped()) {
        zero_front_qualities(read, num_bases_);
    } else {
        zero_back_qualities(read, num_bases_);
    }
}

MaskLowQualityTails::MaskLowQualityTails(BaseQuality threshold) : threshold_ {threshold} {}

void MaskLowQualityTails::operator()(AlignedRead& read) const noexcept
{
    auto& qualities = read.base_qualities();
    const auto is_low_quality = [this] (BaseQuality q) noexcept { return q < threshold_; };
    if (read.is_marked_reverse_mapped()) {
        const auto first_high_quality = std::find_if_not(std::begin(qualities), std::end(qualities), is_low_quality);
        std::fill(std::begin(qualities), first_high_quality, 0);
    } else {
        const auto first_high_quality = std::find_if_not(std::rbegin(qualities), std::rend(qualities), is_low_quality);
        std::fill(std::rbegin(qualities), first_high_quality, 0);
    }
}

void MaskSoftClipped::operator()(AlignedRead& read) const noexcept
{
    if (is_soft_clipped(read)) {
        const auto p = get_soft_clipped_sizes(read);
        zero_front_qualities(read, p.first);
        zero_back_qualities(read, p.second);
    }
}

MaskSoftClippedBoundraryBases::MaskSoftClippedBoundraryBases(Length num_bases) : num_bases_ {num_bases} {}

void MaskSoftClippedBoundraryBases::operator()(AlignedRead& read) const noexcept
{
    if (is_soft_clipped(read)) {
        Length num_front_bases, num_back_bases;
        std::tie(num_front_bases, num_back_bases) = get_soft_clipped_sizes(read);
        if (num_front_bases > 0) {
            zero_front_qualities(read, num_front_bases + num_bases_);
        }
        if (num_back_bases > 0) {
            zero_back_qualities(read, num_back_bases + num_bases_);
        }
    }
}

namespace {

template<typename InputIterator>
void zero_if_less_than(InputIterator first, InputIterator last,
                       typename std::iterator_traits<InputIterator>::value_type value) noexcept {
    std::transform(first, last, first, [value] (auto v) noexcept { return v < value ? 0 : v; });
}

void mask_low_quality_front_bases(AlignedRead& read, std::size_t num_bases, AlignedRead::BaseQuality min_quality) noexcept
{
    auto& qualities = read.base_qualities();
    zero_if_less_than(std::begin(qualities), std::next(std::begin(qualities), std::min(num_bases, sequence_size(read))), min_quality);
}

void mask_low_quality_back_bases(AlignedRead& read, std::size_t num_bases, AlignedRead::BaseQuality min_quality) noexcept
{
    auto& qualities = read.base_qualities();
    zero_if_less_than(std::rbegin(qualities), std::next(std::rbegin(qualities), std::min(num_bases, sequence_size(read))), min_quality);
}

} // namespace

MaskLowQualitySoftClippedBases::MaskLowQualitySoftClippedBases(BaseQuality max) : max_ {max} {}

void MaskLowQualitySoftClippedBases::operator()(AlignedRead& read) const noexcept
{
    if (is_soft_clipped(read)) {
        const auto p = get_soft_clipped_sizes(read);
        mask_low_quality_front_bases(read, p.first, max_);
        mask_low_quality_back_bases(read, p.second, max_);
    }
}

MaskLowQualitySoftClippedBoundaryBases::MaskLowQualitySoftClippedBoundaryBases(Length num_bases, BaseQuality max)
: num_bases_ {num_bases}
, max_ {max}
{}

void MaskLowQualitySoftClippedBoundaryBases::operator()(AlignedRead& read) const noexcept
{
    if (is_soft_clipped(read)) {
        Length num_front_bases, num_back_bases;
        std::tie(num_front_bases, num_back_bases) = get_soft_clipped_sizes(read);
        if (num_front_bases > 0) {
            mask_low_quality_front_bases(read, num_front_bases + num_bases_, max_);
        }
        if (num_back_bases > 0) {
            mask_low_quality_back_bases(read, num_back_bases + num_bases_, max_);
        }
    }
}

MaskLowAverageQualitySoftClippedTails::MaskLowAverageQualitySoftClippedTails(BaseQuality threshold, Length min_tail_length)
: threshold_ {threshold}
, min_tail_length_ {min_tail_length}
{}

namespace {

auto get_soft_clip_tail_size(const AlignedRead& read) noexcept
{
    CigarOperation::Size front_size, back_size;
    std::tie(front_size, back_size) = get_soft_clipped_sizes(read);
    if (read.is_marked_reverse_mapped()) {
        return front_size;
    } else {
        return back_size;
    }
}

template <typename ForwardIt>
auto mean_quality(const ForwardIt first, const std::size_t num_bases) noexcept
{
    const auto last = std::next(first, num_bases);
    return maths::mean(first, last);
}

auto mean_tail_quality(const AlignedRead& read, const std::size_t num_bases) noexcept
{
    assert(num_bases > 0);
    if (read.is_marked_reverse_mapped()) {
        return mean_quality(std::cbegin(read.base_qualities()), num_bases);
    } else {
        return mean_quality(std::crbegin(read.base_qualities()), num_bases);
    }
}

void zero_tail_base_qualities(AlignedRead& read, const std::size_t num_bases) noexcept
{
    if (read.is_marked_reverse_mapped()) {
        zero_front_qualities(read, num_bases);
    } else {
        zero_back_qualities(read, num_bases);
    }
}

} // namespace

void MaskLowAverageQualitySoftClippedTails::operator()(AlignedRead& read) const noexcept
{
    const auto tail_clip_size = get_soft_clip_tail_size(read);
    if (tail_clip_size >= min_tail_length_) {
        const auto mean_quality = mean_tail_quality(read, tail_clip_size);
        if (mean_quality < threshold_) {
            zero_tail_base_qualities(read, tail_clip_size);
        }
    }
}

// template transforms

void mask_adapter_contamination(AlignedRead& forward, AlignedRead& reverse) noexcept
{
    if (begins_before(reverse, forward)) {
        const auto adapter_region = left_overhang_region(reverse, forward);
        zero_front_qualities(reverse, sequence_size(reverse, adapter_region));
    }
    if (ends_before(reverse, forward)) {
        const auto adapter_region = right_overhang_region(forward, reverse);
        zero_back_qualities(forward, sequence_size(forward, adapter_region));
    }
}

void MaskTemplateAdapters::operator()(ReadReferenceVector& read_template) const
{
    const auto template_size = read_template.size();
    if (template_size < 2) {
        return;
    } else if (template_size == 2) {
        if (read_template.front().get().is_marked_reverse_mapped()) {
            mask_adapter_contamination(read_template.back(), read_template.front());
        } else {
            mask_adapter_contamination(read_template.front(), read_template.back());
        }
    } else {
        // TODO
    }
}

namespace {

void
mask_strand_of_duplicated_bases(AlignedRead& forward, AlignedRead& reverse, const GenomicRegion& duplicated_region) {
    auto forward_qual_itr = std::rbegin(forward.base_qualities());
    if (ends_before(reverse, forward)) {
        const auto adapter_region = right_overhang_region(forward, reverse);
        const auto num_adapter_bps = sequence_size(forward, adapter_region);
        forward_qual_itr += num_adapter_bps;
    }
    auto reverse_qual_itr = std::begin(reverse.base_qualities());
    if (begins_before(reverse, forward)) {
        const auto adapter_region = left_overhang_region(reverse, forward);
        const auto num_adapter_bps = sequence_size(reverse, adapter_region);
        reverse_qual_itr += num_adapter_bps;
    }
    const auto num_forward_duplicate_bps = sequence_size(forward, duplicated_region);
    const auto num_reverse_duplicate_bps = sequence_size(reverse, duplicated_region);
    auto num_duplicate_bps = std::min(num_forward_duplicate_bps, num_reverse_duplicate_bps);
    bool select_first{true};
    for (; num_duplicate_bps > 0; --num_duplicate_bps) {
        assert(forward_qual_itr < std::rend(forward.base_qualities()));
        assert(reverse_qual_itr < std::end(reverse.base_qualities()));
        if (*forward_qual_itr == *reverse_qual_itr) {
            if (select_first) {
                *reverse_qual_itr++ = 0;
                select_first = false;
            } else {
                *forward_qual_itr++ = 0;
                select_first = true;
            }
        } else if (*forward_qual_itr < *reverse_qual_itr) {
            *forward_qual_itr++ = 0;
        } else {
            *reverse_qual_itr++ = 0;
        }
    }
}

void mask_strand_of_duplicated_bases(AlignedRead& forward, AlignedRead& reverse) {
    const auto duplicated_region = overlapped_region(forward, reverse);
    if (duplicated_region) {
        mask_strand_of_duplicated_bases(forward, reverse, *duplicated_region);
    }
}

} // namespace

void MaskStrandOfDuplicatedBases::operator()(ReadReferenceVector& read_template) const
{
    const auto template_size = read_template.size();
    if (template_size < 2) {
        return;
    } else if (template_size == 2) {
        if (read_template.front().get().is_marked_reverse_mapped()) {
            mask_strand_of_duplicated_bases(read_template.back(), read_template.front());
        } else {
            mask_strand_of_duplicated_bases(read_template.front(), read_template.back());
        }
    } else {
        // TODO
    }
}

namespace {

boost::optional<GenomicRegion>
duplicated_clipped_region(const AlignedRead& first, const AlignedRead& second, const GenomicRegion& duplicated_region)
{
    if (is_front_soft_clipped(first) && is_front_soft_clipped(second)) {
        const auto first_clipping  = splice_cigar(first, duplicated_region).front();
        const auto second_clipping = splice_cigar(second, duplicated_region).front();
        assert(is_clipping(first_clipping) && is_clipping(second_clipping));
        const auto num_duplicate_clipped_bases = std::min(first_clipping.size(), second_clipping.size());
        return expand_rhs(head_region(duplicated_region), num_duplicate_clipped_bases);
    } else if (is_back_soft_clipped(first) && is_back_soft_clipped(second)) {
        const auto first_clipping  = splice_cigar(first, duplicated_region).back();
        const auto second_clipping = splice_cigar(second, duplicated_region).back();
        assert(is_clipping(first_clipping) && is_clipping(second_clipping));
        const auto num_duplicate_clipped_bases = std::min(first_clipping.size(), second_clipping.size());
        return expand_lhs(tail_region(duplicated_region), num_duplicate_clipped_bases);
    } else {
        return boost::none;
    }
}

void mask(AlignedRead& read, const GenomicRegion& region)
{
    if (!is_empty(region) && overlaps(read, region)) {
        auto& base_qualities = read.base_qualities();
        if (contains(region, read)) {
            std::fill(std::begin(base_qualities), std::end(base_qualities), 0);
        } else if (begins_equal(read, region)) {
            const auto num_mask_bases = sequence_size(read, region);
            std::fill_n(std::begin(base_qualities), num_mask_bases, 0);
        } else if (ends_equal(read, region)) {
            const auto num_mask_bases = sequence_size(read, region);
            std::fill_n(std::rbegin(base_qualities), num_mask_bases, 0);
        } else {
            // Take the shortest side for optimisation only
            if (begin_distance(read, region) < end_distance(region, read)) {
                const auto num_lhs_mask_bases = sequence_size(read, left_overhang_region(read, region));
                const auto num_mask_bases = sequence_size(read, region);
                std::fill_n(std::next(std::begin(base_qualities), num_lhs_mask_bases), num_mask_bases, 0);
            } else {
                const auto num_rhs_mask_bases = sequence_size(read, right_overhang_region(read, region));
                const auto num_mask_bases = sequence_size(read, region);
                std::fill_n(std::next(std::rbegin(base_qualities), num_rhs_mask_bases), num_mask_bases, 0);
            }
        }
    }
}

void mask_both_strands(AlignedRead& forward, AlignedRead& reverse, const GenomicRegion& region)
{
    mask(forward, region);
    mask(reverse, region);
}

void mask_both_strands_of_clipped_duplicated_bases(AlignedRead& forward, AlignedRead& reverse)
{
    const auto duplicated_region = overlapped_region(forward, reverse);
    if (duplicated_region) {
        const auto mask_region = duplicated_clipped_region(forward, reverse, *duplicated_region);
        if (mask_region) {
            mask_both_strands(forward, reverse, *mask_region);
        }
    }
}

} // namespace

void MaskClippedDuplicatedBases::operator()(ReadReferenceVector& read_template) const
{
    const auto template_size = read_template.size();
    if (template_size < 2) {
        return;
    } else if (template_size == 2) {
        if (read_template.front().get().is_marked_reverse_mapped()) {
            mask_both_strands_of_clipped_duplicated_bases(read_template.back(), read_template.front());
        } else {
            mask_both_strands_of_clipped_duplicated_bases(read_template.front(), read_template.back());
        }
    } else {
        // TODO
    }
}

} // namespace readpipe
} // namespace octopus
