//
//  genomic_region.h
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__genomic_region__
#define __Octopus__genomic_region__

#include <string>
#include <cstdint>
#include <ostream>

#include "sequence_region.h"
#include "comparable.h"

using std::uint_fast32_t;
using std::int_fast64_t;

/**
    Represents a continuous region of a sequence in a genome. The sequence
    name is the reference sequence name (usually a chromosome), and the
    begin and end positions are zero-indexed half open - [begin,end) - indexes.
 */
class GenomicRegion : Comparable<GenomicRegion>
{
public:
    GenomicRegion() = delete;
    GenomicRegion(std::string contig_name, uint_fast32_t begin, uint_fast32_t end);
    //GenomicRegion(std::string the_region);
    
    GenomicRegion(const GenomicRegion&)            = default;
    GenomicRegion& operator=(const GenomicRegion&) = default;
    GenomicRegion(GenomicRegion&&)                 = default;
    GenomicRegion& operator=(GenomicRegion&&)      = default;
    
    const std::string& get_contig_name() const noexcept;
    const SequenceRegion& get_contig_region() const noexcept;
    uint_fast32_t get_begin() const noexcept;
    uint_fast32_t get_end() const noexcept;

private:
    std::string contig_name_;
    SequenceRegion contig_region_;
};

inline
GenomicRegion::GenomicRegion(std::string contig_name, uint_fast32_t begin, uint_fast32_t end)
:contig_name_ {contig_name},
 contig_region_ {begin, end}
{}

//GenomicRegion::GenomicRegion(std::string the_region)
//{
//    // parse region
//}

inline const std::string& GenomicRegion::get_contig_name() const noexcept
{
    return contig_name_;
}

inline const SequenceRegion& GenomicRegion::get_contig_region() const noexcept
{
    return contig_region_;
}

inline uint_fast32_t GenomicRegion::get_begin() const noexcept
{
    return contig_region_.get_begin();
}

inline uint_fast32_t GenomicRegion::get_end() const noexcept
{
    return contig_region_.get_end();
}

inline uint_fast32_t size(const GenomicRegion& a_region) noexcept
{
    return size(a_region.get_contig_region());
}

inline bool is_same_contig(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return lhs.get_contig_name() == rhs.get_contig_name();
}

inline int_fast64_t overlap_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) ? overlap_size(lhs.get_contig_region(), rhs.get_contig_region()) : 0;
}

inline bool overlaps(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return overlap_size(lhs, rhs) > 0;
}

inline bool operator==(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return is_same_contig(lhs, rhs) && lhs.get_contig_region() == rhs.get_contig_region();
}

// It doesn't really make sense to define ordering operators for GenomicRegion
// (as oposed to SequenceRegion), as non-continuous sequences have no natural ordering.
inline bool operator<(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) return lhs.get_contig_region() < rhs.get_contig_region();
    throw std::logic_error {"Cannot compare regions on different contigs"};
}

inline std::string to_string(const GenomicRegion& a_region)
{
    return a_region.get_contig_name() + ':' + std::to_string(a_region.get_begin()) + '-'
            + std::to_string(a_region.get_end());
}

namespace std {
    template <> struct hash<GenomicRegion>
    {
        size_t operator()(const GenomicRegion& r) const
        {
            return hash<std::string>()(to_string(r));
        }
    };
}

inline
std::ostream& operator<<(std::ostream& os, const GenomicRegion& a_region)
{
    os << to_string(a_region);
    return os;
}

#endif /* defined(__Octopus__genomic_region__) */
