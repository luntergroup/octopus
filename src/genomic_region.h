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
#include <ostream>
#include <stdexcept>
#include <vector>

#include "sequence_region.h"
#include "comparable.h"

/**
    Represents a continuous region of a sequence in a genome. The sequence
    name is the reference sequence name (usually a chromosome), and the
    begin and end positions are zero-indexed half open - [begin,end) - indices.
 */
class GenomicRegion : public Comparable<GenomicRegion>
{
public:
    using StringType     = std::string;
    using SizeType       = SequenceRegion::SizeType;
    using DifferenceType = SequenceRegion::DifferenceType;
    
    GenomicRegion() = default;
    template <typename T>
    explicit GenomicRegion(T&& contig_name, SizeType begin, SizeType end);
    
    GenomicRegion(const GenomicRegion&)            = default;
    GenomicRegion& operator=(const GenomicRegion&) = default;
    GenomicRegion(GenomicRegion&&)                 = default;
    GenomicRegion& operator=(GenomicRegion&&)      = default;
    
    const StringType& get_contig_name() const noexcept;
    const SequenceRegion& get_contig_region() const noexcept;
    SizeType get_begin() const noexcept;
    SizeType get_end() const noexcept;

private:
    StringType contig_name_;
    SequenceRegion contig_region_;
};

template <typename T>
GenomicRegion::GenomicRegion(T&& contig_name, SizeType begin, SizeType end)
:contig_name_ {std::forward<T>(contig_name)},
 contig_region_ {begin, end}
{}

inline const GenomicRegion::StringType& GenomicRegion::get_contig_name() const noexcept
{
    return contig_name_;
}

inline const SequenceRegion& GenomicRegion::get_contig_region() const noexcept
{
    return contig_region_;
}

inline GenomicRegion::SizeType GenomicRegion::get_begin() const noexcept
{
    return contig_region_.get_begin();
}

inline GenomicRegion::SizeType GenomicRegion::get_end() const noexcept
{
    return contig_region_.get_end();
}

inline GenomicRegion::SizeType size(const GenomicRegion& a_region) noexcept
{
    return size(a_region.get_contig_region());
}

inline bool is_same_contig(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return lhs.get_contig_name() == rhs.get_contig_name();
}

inline GenomicRegion::DifferenceType overlap_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) ? overlap_size(lhs.get_contig_region(), rhs.get_contig_region()) : 0;
}

inline bool overlaps(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return overlap_size(lhs, rhs) > 0;
}

inline std::string to_string(const GenomicRegion& a_region)
{
    return a_region.get_contig_name() + ':' + std::to_string(a_region.get_begin()) + '-'
    + std::to_string(a_region.get_end());
}

inline bool operator==(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return is_same_contig(lhs, rhs) && lhs.get_contig_region() == rhs.get_contig_region();
}

inline bool operator<(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) return lhs.get_contig_region() < rhs.get_contig_region();
    throw std::runtime_error {"Cannot compare regions on different contigs: "  + to_string(lhs) +
        " & " + to_string(rhs)};
}

inline GenomicRegion shift(const GenomicRegion& a_region, GenomicRegion::DifferenceType n)
{
    if (n < 0 && a_region.get_begin() + n > a_region.get_begin()) {
        throw std::out_of_range {"Shifted past contig start"};
    }
    
    return GenomicRegion {
        a_region.get_contig_name(),
        static_cast<GenomicRegion::SizeType>(a_region.get_begin() + n),
        static_cast<GenomicRegion::SizeType>(a_region.get_end() + n)
    };
}

inline GenomicRegion get_intervening_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (lhs >= rhs || overlaps(lhs, rhs)) {
        throw std::runtime_error {"Cannot get intervening region between overlapping regions"};
    }
    
    return GenomicRegion {
        lhs.get_contig_name(),
        lhs.get_end(),
        rhs.get_begin()
    };
}

//inline std::vector<GenomicRegion> split(const GenomicRegion& a_region, GenomicRegion::SizeType max_size)
//{
//    std::vector<GenomicRegion> result {};
//    
//    return result;
//}

namespace std {
    template <> struct hash<GenomicRegion>
    {
        size_t operator()(const GenomicRegion& r) const
        {
            return hash<string>()(to_string(r));
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const GenomicRegion& a_region)
{
    os << to_string(a_region);
    return os;
}

#endif /* defined(__Octopus__genomic_region__) */
