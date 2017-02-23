// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_transform.hpp"

#include <algorithm>
#include <iterator>
#include <tuple>

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
    auto& qualities = read.qualities();
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
    std::transform(first, last, first, [value] (auto v) noexcept { return v < value ? 0 : value; });
}

void mask_low_quality_front_bases(AlignedRead& read, std::size_t num_bases, AlignedRead::BaseQuality min_quality) noexcept
{
    auto& qualities = read.qualities();
    zero_if_less_than(std::begin(qualities), std::next(std::begin(qualities), std::min(num_bases, sequence_size(read))), min_quality);
}

void mask_low_quality_back_bases(AlignedRead& read, std::size_t num_bases, AlignedRead::BaseQuality min_quality) noexcept
{
    auto& qualities = read.qualities();
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

// template transforms

void mask_adapter_contamination(AlignedRead& first, AlignedRead& second) noexcept
{
    if (begins_before(second, first)) {
        const auto num_adapter_bases = begin_distance(second, first);
        zero_front_qualities(second, num_adapter_bases);
    }
    if (ends_before(second, first)) {
        const auto num_adapter_bases = end_distance(second, first);
        zero_back_qualities(first, num_adapter_bases);
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

void mask_duplicated_bases(AlignedRead& first, AlignedRead& second, const GenomicRegion& duplicated_region)
{
    auto first_qual_itr = std::rbegin(first.qualities());
    if (ends_before(second, first)) {
        first_qual_itr += end_distance(duplicated_region, first);
    }
    auto second_qual_itr = std::begin(second.qualities());
    if (begins_before(second, first)) {
        second_qual_itr += begin_distance(second, duplicated_region);
    }
    auto num_duplicate_bases = static_cast<int>(size(duplicated_region));
    bool select_first {true};
    for (; num_duplicate_bases > 0; --num_duplicate_bases) {
        if (*first_qual_itr == *second_qual_itr) {
            if (select_first) {
                *second_qual_itr++ = 0;
                select_first = false;
            } else {
                *first_qual_itr++ = 0;
                select_first = true;
            }
        } else if (*first_qual_itr < *second_qual_itr) {
            *first_qual_itr++ = 0;
        } else {
            *second_qual_itr++ = 0;
        }
    }
}

void mask_duplicated_bases(AlignedRead& first, AlignedRead& second)
{
    const auto duplicated_region = overlapped_region(first, second);
    if (duplicated_region) {
        mask_duplicated_bases(first, second, *duplicated_region);
    }
}

void MaskDuplicatedBases::operator()(ReadReferenceVector& read_template) const
{
    const auto template_size = read_template.size();
    if (template_size < 2) {
        return;
    } else if (template_size == 2) {
        if (read_template.front().get().is_marked_reverse_mapped()) {
            mask_duplicated_bases(read_template.back(), read_template.front());
        } else {
            mask_duplicated_bases(read_template.front(), read_template.back());
        }
    } else {
        // TODO
    }
}

} // namespace readpipe
} // namespace octopus
