//
//  genomic_region.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "genomic_region.h"

#include <cmath>

GenomicRegion::GenomicRegion(std::string contig_name, size_t begin, size_t end) noexcept
: contig_name {contig_name}, begin {begin}, end {end}
{}

size_t GenomicRegion::size() const noexcept
{
    return end - begin;
}

size_t GenomicRegion::num_overlaped_bases(const GenomicRegion& other) const noexcept
{
    return std::min(begin, other.begin) - std::max(end, other.end);
}

bool GenomicRegion::overlaps(const GenomicRegion& other) const noexcept
{
    return num_overlaped_bases(other) > 0;
}

inline bool operator==(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return lhs.contig_name == rhs.contig_name && lhs.begin == rhs.begin && lhs.end == rhs.end;
}
inline bool operator!=(const GenomicRegion& lhs, const GenomicRegion& rhs);
inline bool operator< (const GenomicRegion& lhs, const GenomicRegion& rhs);
inline bool operator> (const GenomicRegion& lhs, const GenomicRegion& rhs);
inline bool operator<=(const GenomicRegion& lhs, const GenomicRegion& rhs);
inline bool operator>=(const GenomicRegion& lhs, const GenomicRegion& rhs);
