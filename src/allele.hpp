//
//  allele.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_allele_hpp
#define Octopus_allele_hpp

#include <string>
#include <ostream>
#include <utility> // std::forward
#include <boost/functional/hash.hpp> // boost::hash_combine

#include "genomic_region.hpp"
#include "comparable.hpp"
#include "mappable.hpp"
#include "reference_genome.hpp"

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

bool is_reference(const Allele& allele, ReferenceGenome& reference);

Allele get_reference_allele(const GenomicRegion& region, ReferenceGenome& reference);

std::vector<Allele> get_reference_alleles(const std::vector<GenomicRegion>& regions, ReferenceGenome& reference);

std::vector<Allele> get_positional_reference_alleles(const GenomicRegion& region, ReferenceGenome& reference);

//Allele::SequenceType get_subsequence(const Allele& allele, const GenomicRegion& region);

Allele splice(const Allele& allele, const GenomicRegion& region);

bool contains(const Allele& lhs, const Allele& rhs);

bool is_insertion(const Allele& allele);

bool is_deletion(const Allele& allele);

std::vector<Allele> decompose(const Allele& allele);

bool operator==(const Allele& lhs, const Allele& rhs);
bool operator<(const Allele& lhs, const Allele& rhs);

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
} // end namespace std

std::ostream& operator<<(std::ostream& os, const Allele& allele);

#endif
