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

MaskSoftClippedBoundries::MaskSoftClippedBoundries(Length num_bases) : num_bases_ {num_bases} {}

void MaskSoftClippedBoundries::operator()(AlignedRead& read) const noexcept
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
    
}

void MaskLowQualitySoftClippedBases::operator()(AlignedRead& read) const noexcept
{
    if (is_soft_clipped(read)) {
        const auto p = get_soft_clipped_sizes(read);
        auto& qualities = read.qualities();
        zero_if_less_than(std::begin(qualities), std::next(std::begin(qualities), p.first), max_);
        zero_if_less_than(std::rbegin(qualities), std::next(std::rbegin(qualities), p.second), max_);
    }
}

void QualityAdjustedSoftClippedMasker::operator()(AlignedRead& read) const noexcept
{
    if (is_soft_clipped(read)) {
        using std::cbegin; using std::crbegin; using std::next;
        using std::accumulate; using std::min_element; using std::min;
        
        using S = AlignedRead::NucleotideSequence::size_type;
        S num_front_bases, num_back_bases;
        std::tie(num_front_bases, num_back_bases) = get_soft_clipped_sizes(read);
        const auto& qualities = read.qualities();
        
        using Q = AlignedRead::BaseQuality;
        
        if (num_front_bases > 0) {
            const auto sum = accumulate(cbegin(qualities), next(cbegin(qualities), num_front_bases), 0.0);
            const auto mean = static_cast<Q>(sum / num_front_bases);
            const auto min_quality = *min_element(cbegin(qualities), next(cbegin(qualities)));
            const auto mask_size = num_front_bases + min(static_cast<S>(mean - min_quality), num_front_bases);
            set_front_qualities(read, mask_size);
        }
        
        if (num_back_bases > 0) {
            const auto sum = accumulate(crbegin(qualities), next(crbegin(qualities), num_back_bases), 0.0);
            const auto mean = static_cast<Q>(sum / num_back_bases);
            const auto min_quality = *min_element(cbegin(qualities), next(cbegin(qualities)));
            const auto mask_size = num_back_bases + min(static_cast<S>(mean - min_quality), num_back_bases);
            set_back_qualities(read, mask_size);
        }
    }
}

} // namespace readpipe
} // namespace octopus
