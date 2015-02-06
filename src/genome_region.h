//
//  genome_region.h
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_genome_region_h
#define Octopus_genome_region_h

#include <string>
#include <cstddef>

using std::size_t;

struct GenomeRegion
{
    std::string contig_name;
    size_t begin, end;
    
    GenomeRegion() = delete;
    GenomeRegion(std::string contig_name, size_t begin, size_t end) noexcept;
    
    size_t size() const noexcept;
    size_t num_overlaped_bases(const GenomeRegion& other) const noexcept;
    bool overlaps(const GenomeRegion& other) const noexcept;
};

inline bool operator==(const GenomeRegion& lhs, const GenomeRegion& rhs);
inline bool operator!=(const GenomeRegion& lhs, const GenomeRegion& rhs);
inline bool operator< (const GenomeRegion& lhs, const GenomeRegion& rhs);
inline bool operator> (const GenomeRegion& lhs, const GenomeRegion& rhs);
inline bool operator<=(const GenomeRegion& lhs, const GenomeRegion& rhs);
inline bool operator>=(const GenomeRegion& lhs, const GenomeRegion& rhs);

#endif
