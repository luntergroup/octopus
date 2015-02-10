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

#include "sequence_region.h"

using std::uint_fast32_t;
using std::int_fast64_t;

/*
    Represents a continuous region of a sequence in a genome. The sequence
    name is the reference sequence name (usually a chromosome), and the
    begin and end positions are zero-indexed half open - [begin,end) - indexes.
 */
class GenomicRegion
{
public:
    GenomicRegion() = delete;
    GenomicRegion(std::string contig_name, uint_fast32_t begin_pos, uint_fast32_t end_pos);
    //GenomicRegion(std::string the_region);
    
    const std::string& get_contig_name() const noexcept;
    const SequenceRegion& get_contig_region() const noexcept;
    uint_fast32_t get_begin_pos() const noexcept;
    uint_fast32_t get_end_pos() const noexcept;

private:
    const std::string contig_name_;
    const SequenceRegion region_;
};

inline
GenomicRegion::GenomicRegion(std::string contig_name, uint_fast32_t begin_pos, uint_fast32_t end_pos)
: contig_name_ {contig_name}, region_(begin_pos, end_pos)
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
    return region_;
}

inline uint_fast32_t GenomicRegion::get_begin_pos() const noexcept
{
    return region_.get_begin_pos();
}

inline uint_fast32_t GenomicRegion::get_end_pos() const noexcept
{
    return region_.get_end_pos();
}

inline uint_fast32_t size(const GenomicRegion& a_region) noexcept
{
    return size(a_region.get_contig_region());
}

inline std::string to_string(const GenomicRegion& a_region)
{
    return "";
}

inline int_fast64_t overlap_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return (lhs.get_contig_name() == rhs.get_contig_name()) ?
        overlap_size(lhs.get_contig_region(), rhs.get_contig_region()) : 0;
}

inline bool overlaps(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return overlap_size(lhs, rhs) > 0;
}

// It doesn't really make sense to define ordering operators for GenomicRegion
// (as oposed to SequenceRegion), as non-continuous sequences have no natural ordering.
inline bool operator==(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return lhs.get_contig_name() == rhs.get_contig_name() &&
            lhs.get_contig_region() == rhs.get_contig_region();
}
inline bool operator!=(const GenomicRegion& lhs, const GenomicRegion& rhs) {return !operator==(lhs, rhs);}

#endif /* defined(__Octopus__genomic_region__) */
