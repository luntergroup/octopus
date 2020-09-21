// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_duplicates.hpp"

#include <cassert>

namespace octopus {

namespace {

bool are_duplicates(const AlignedRead::Segment& lhs, const AlignedRead::Segment& rhs) noexcept
{
    return lhs.contig_name() == rhs.contig_name()
           && lhs.is_marked_unmapped() == rhs.is_marked_unmapped()
           && lhs.is_marked_reverse_mapped() == rhs.is_marked_reverse_mapped()
           && lhs.inferred_template_length() == rhs.inferred_template_length();
}

bool next_segments_are_duplicates(const AlignedRead& lhs, const AlignedRead& rhs) noexcept
{
    if (lhs.has_other_segment()) {
        return rhs.has_other_segment() && are_duplicates(lhs.next_segment(), rhs.next_segment());
    } else {
        return !rhs.has_other_segment();
    }
}

} // namespace

bool FivePrimeDuplicateDefinition::unpaired_equal(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
{
    return are_same_strand(lhs, rhs)
        && is_same_contig(lhs, rhs)
		&& five_prime_mapping_position(lhs) == five_prime_mapping_position(rhs);
}

bool FivePrimeDuplicateDefinition::paired_equal(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
{
    return unpaired_equal(lhs, rhs) && next_segments_are_duplicates(lhs, rhs);
}

bool FivePrimeAndCigarDuplicateDefinition::unpaired_equal(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
{
    return are_same_strand(lhs, rhs)
        && is_same_region(lhs, rhs)
        && lhs.cigar() == rhs.cigar();
}

bool FivePrimeAndCigarDuplicateDefinition::paired_equal(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
{
    return unpaired_equal(lhs, rhs) && next_segments_are_duplicates(lhs, rhs);
}

namespace detail {

bool NextSegmentLess::operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
{
    assert(lhs.has_other_segment() && rhs.has_other_segment());
    return lhs.next_segment().inferred_template_length() < rhs.next_segment().inferred_template_length();
}

} // namespace

} // namespace octopus
