//
//  allele.h
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_allele_h
#define Octopus_allele_h

#include <string>
#include <ostream>
#include <utility> // std::forward
#include <boost/functional/hash.hpp> // boost::hash_combine

#include "genomic_region.h"
#include "comparable.h"
#include "mappable.h"
#include "reference_genome.h"

class Allele : public Comparable<Allele>, public Mappable<Allele>
{
public:
    using SizeType     = GenomicRegion::SizeType;
    using SequenceType = ReferenceGenome::SequenceType;
    
    Allele() = default;
    template <typename GenomicRegion_, typename SequenceType_>
    Allele(GenomicRegion_&& reference_region, SequenceType_&& sequence);
    ~Allele() = default;
    
    Allele(const Allele&)            = default;
    Allele& operator=(const Allele&) = default;
    Allele(Allele&&)                 = default;
    Allele& operator=(Allele&&)      = default;
    
    const GenomicRegion& get_region() const noexcept;
    const SequenceType& get_sequence() const noexcept;
    
private:
    GenomicRegion reference_region_;
    SequenceType sequence_;
};

template <typename GenomicRegion_, typename SequenceType_>
Allele::Allele(GenomicRegion_&& reference_region, SequenceType_&& sequence)
:
reference_region_ {std::forward<GenomicRegion_>(reference_region)},
sequence_ {std::forward<SequenceType_>(sequence)}
{}

bool is_reference(const Allele& an_allele, ReferenceGenome& the_reference);

Allele get_reference_allele(const GenomicRegion& region, ReferenceGenome& the_reference);

Allele::SequenceType get_subsequence(const Allele& an_allele, const GenomicRegion& region);

bool contains(const Allele& lhs, const Allele& rhs);

Allele splice(const Allele& an_allele, const GenomicRegion& region);

bool is_insertion(const Allele& an_allele);

bool is_deletion(const Allele& an_allele);

std::vector<Allele> decompose(const Allele& an_allele);

inline bool operator==(const Allele& lhs, const Allele& rhs)
{
    return lhs.get_region() == rhs.get_region() && lhs.get_sequence() == rhs.get_sequence();
}

inline bool operator<(const Allele& lhs, const Allele& rhs)
{
    return (lhs.get_region() == rhs.get_region()) ? lhs.get_sequence() < rhs.get_sequence() :
    lhs.get_region() < rhs.get_region();
}

namespace std {
    template <> struct hash<Allele>
    {
        size_t operator()(const Allele& a) const
        {
            size_t seed {};
            boost::hash_combine(seed, hash<GenomicRegion>()(a.get_region()));
            boost::hash_combine(seed, hash<Allele::SequenceType>()(a.get_sequence()));
            return seed;
        }
    };
}

std::ostream& operator<<(std::ostream& os, const Allele& an_allele);

#endif
