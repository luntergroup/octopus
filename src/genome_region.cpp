//
//  genome_region.cpp
//  Octopus
//
//  Created by Daniel Cooke on 06/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "genome_region.h"

#include <cmath>

GenomeRegion::GenomeRegion(std::string contig_name, size_t begin, size_t end) noexcept
: contig_name {contig_name}, begin {begin}, end {end}
{}

size_t GenomeRegion::size() const noexcept
{
    return end - begin;
}

size_t GenomeRegion::num_overlaped_bases(const GenomeRegion& other) const noexcept
{
    return std::min(begin, other.begin) - std::max(end, other.end);
}

bool GenomeRegion::overlaps(const GenomeRegion& other) const noexcept
{
    return num_overlaped_bases(other) > 0;
}

inline bool operator==(const GenomeRegion& lhs, const GenomeRegion& rhs)
{
    return lhs.contig_name == rhs.contig_name && lhs.begin == rhs.begin && lhs.end == rhs.end;
}
inline bool operator!=(const GenomeRegion& lhs, const GenomeRegion& rhs);
inline bool operator< (const GenomeRegion& lhs, const GenomeRegion& rhs);
inline bool operator> (const GenomeRegion& lhs, const GenomeRegion& rhs);
inline bool operator<=(const GenomeRegion& lhs, const GenomeRegion& rhs);
inline bool operator>=(const GenomeRegion& lhs, const GenomeRegion& rhs);
