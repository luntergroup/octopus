// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_transform.hpp"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <tuple>

namespace octopus { namespace readpipe
{

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
            set_back_qualities(read, overlapped_size);
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
                set_front_qualities(read, num_adapter_bases);
            } else {
                set_back_qualities(read, num_adapter_bases);
            }
        }
    }
}

MaskTail::MaskTail(Length num_bases) : num_bases_ {num_bases} {}

void MaskTail::operator()(AlignedRead& read) const noexcept
{
    if (read.is_marked_reverse_mapped()) {
        set_front_qualities(read, num_bases_);
    } else {
        set_back_qualities(read, num_bases_);
    }
}

void MaskSoftClipped::operator()(AlignedRead& read) const noexcept
{
    if (is_soft_clipped(read)) {
        const auto p = get_soft_clipped_sizes(read);
        set_front_qualities(read, p.first);
        set_back_qualities(read, p.second);
    }
}

MaskSoftClippedBoundraryBases::MaskSoftClippedBoundraryBases(Length num_bases) : num_bases_ {num_bases} {}

void MaskSoftClippedBoundraryBases::operator()(AlignedRead& read) const noexcept
{
    if (is_soft_clipped(read)) {
        Length num_front_bases, num_back_bases;
        std::tie(num_front_bases, num_back_bases) = get_soft_clipped_sizes(read);
        
        if (num_front_bases > 0) {
            set_front_qualities(read, num_front_bases + num_bases_);
        }
        if (num_back_bases > 0) {
            set_back_qualities(read, num_back_bases + num_bases_);
        }
    }
}

MaskLowQualitySoftClippedBases::MaskLowQualitySoftClippedBases(BaseQuality max) : max_ {max} {}

namespace {

template<typename InputIterator>
void zero_if_less_than(InputIterator first, InputIterator last,
                       typename std::iterator_traits<InputIterator>::value_type value) noexcept {
    std::transform(first, last, first, [value](auto v) { return v < value ? 0 : value; });
}

void mask_low_quality_front_bases(AlignedRead& read, std::size_t num_bases, AlignedRead::BaseQuality min_quality) noexcept
{
    auto& sequence = read.sequence();
    zero_if_less_than(std::begin(sequence), std::next(std::begin(sequence), num_bases), min_quality);
}

void mask_low_quality_back_bases(AlignedRead& read, std::size_t num_bases, AlignedRead::BaseQuality min_quality) noexcept
{
    auto& qualities = read.qualities();
    zero_if_less_than(std::rbegin(qualities), std::next(std::rbegin(qualities), num_bases), min_quality);
}


}

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

} // namespace readpipe
} // namespace octopus
