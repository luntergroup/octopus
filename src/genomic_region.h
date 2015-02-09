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
#include <cstddef>

using std::size_t;

struct GenomicRegion
{
    std::string contig_name;
    size_t begin, end;
    
    GenomicRegion() = delete;
    GenomicRegion(std::string contig_name, size_t begin, size_t end) noexcept;
    
    size_t size() const noexcept;
    size_t num_overlaped_bases(const GenomicRegion& other) const noexcept;
    bool overlaps(const GenomicRegion& other) const noexcept;
};

inline bool operator==(const GenomicRegion& lhs, const GenomicRegion& rhs);
inline bool operator!=(const GenomicRegion& lhs, const GenomicRegion& rhs);
inline bool operator< (const GenomicRegion& lhs, const GenomicRegion& rhs);
inline bool operator> (const GenomicRegion& lhs, const GenomicRegion& rhs);
inline bool operator<=(const GenomicRegion& lhs, const GenomicRegion& rhs);
inline bool operator>=(const GenomicRegion& lhs, const GenomicRegion& rhs);

#endif /* defined(__Octopus__genomic_region__) */
