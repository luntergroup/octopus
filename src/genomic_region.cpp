//
//  genomic_region.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "genomic_region.h"

GenomicRegion::GenomicRegion(std::string sequence_name, int_fast32_t begin, int_fast32_t end)
: the_sequence_name_ {sequence_name}, the_region_(begin, end)
{}

//GenomicRegion::GenomicRegion(std::string the_region)
//{
//    // parse region
//}

const std::string& GenomicRegion::get_sequence_name() const noexcept
{
    return the_sequence_name_;
}

const SequenceRegion& GenomicRegion::get_sequence_region() const noexcept
{
    return the_region_;
}

int_fast32_t GenomicRegion::get_begin_pos() const noexcept
{
    return the_region_.begin_pos;
}

int_fast32_t GenomicRegion::get_end_pos() const noexcept
{
    return the_region_.end_pos;
}

int_fast32_t size(const GenomicRegion& a_region) noexcept
{
    return size(a_region.get_sequence_region());
}

inline bool operator==(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return lhs.get_sequence_name() == rhs.get_sequence_name() && lhs.get_sequence_region() == rhs.get_sequence_region();
}
inline bool operator!=(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return !operator==(lhs, rhs);
}

std::string to_string(const GenomicRegion& a_region)
{
    return "";
}

int_fast32_t overlap_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    if (lhs.get_sequence_name() == rhs.get_sequence_name()) {
        return overlap_size(lhs.get_sequence_region(), rhs.get_sequence_region());
    } else {
        return 0;
    }
}

bool overlaps(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return overlap_size(lhs, rhs) > 0;
}
