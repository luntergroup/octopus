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
#include <algorithm> // std::min
#include <cstddef>   // std::size_t

#include "genomic_region.h"
#include "comparable.h"
#include "mappable.h"
#include "reference_genome.h"
#include "string_utils.h"

class Allele : public Comparable<Allele>, public Mappable<Allele>
{
public:
    using SizeType     = GenomicRegion::SizeType;
    using SequenceType = ReferenceGenome::SequenceType;
    using RealType     = double;
    
    Allele() = default;
    template <typename GenomicRegion_, typename SequenceType_>
    Allele(GenomicRegion_&& the_reference_region, SequenceType_&& the_sequence);
    template <typename GenomicRegion_, typename SequenceType_>
    Allele(GenomicRegion_&& the_reference_region, SequenceType_&& the_sequence, RealType probability);
    ~Allele() = default;
    
    Allele(const Allele&)            = default;
    Allele& operator=(const Allele&) = default;
    Allele(Allele&&)                 = default;
    Allele& operator=(Allele&&)      = default;
    
    const GenomicRegion& get_region() const noexcept;
    const SequenceType& get_sequence() const noexcept;
    
    void set_probability(RealType probability) noexcept;
    RealType get_probability() const noexcept;
    
private:
    GenomicRegion the_reference_region_;
    SequenceType the_sequence_;
    RealType probability_;
};

template <typename GenomicRegion_, typename SequenceType_>
Allele::Allele(GenomicRegion_&& the_reference_region, SequenceType_&& the_sequence)
:
the_reference_region_ {std::forward<GenomicRegion_>(the_reference_region)},
the_sequence_ {std::forward<SequenceType_>(the_sequence)},
probability_ {0}
{}

template <typename GenomicRegion_, typename SequenceType_>
Allele::Allele(GenomicRegion_&& the_reference_region, SequenceType_&& the_sequence, RealType probability)
:
the_reference_region_ {std::forward<GenomicRegion_>(the_reference_region)},
the_sequence_ {std::forward<SequenceType_>(the_sequence)},
probability_ {probability}
{}

inline const GenomicRegion& Allele::get_region() const noexcept
{
    return the_reference_region_;
}

inline const Allele::SequenceType& Allele::get_sequence() const noexcept
{
    return the_sequence_;
}

inline void Allele::set_probability(RealType probability) noexcept
{
    probability_ = probability;
}

inline Allele::RealType Allele::get_probability() const noexcept
{
    return probability_;
}

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
            return hash<GenomicRegion>()(a.get_region());
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const Allele& an_allele)
{
    os << an_allele.get_region() << " " << an_allele.get_sequence();
    return os;
}

inline bool is_reference_allele(const Allele& an_allele, ReferenceGenome& the_reference)
{
    return an_allele.get_sequence() == the_reference.get_sequence(an_allele.get_region());
}

inline Allele get_reference_allele(const GenomicRegion& a_region, ReferenceGenome& the_reference)
{
    return Allele {a_region, the_reference.get_sequence(a_region)};
}

inline Allele::SequenceType get_subsequence(const Allele& an_allele, const GenomicRegion& a_region)
{
    if (!contains(an_allele, a_region)) {
        return Allele::SequenceType {};
    } else {
        auto first = std::cbegin(an_allele.get_sequence()) + (get_begin(a_region) - get_begin(an_allele));
        
        // The minimum of the allele sequence size and region size is used as deletions will
        // result in a sequence size smaller than the region size
        return Allele::SequenceType {first, first +
            std::min(an_allele.get_sequence().size(), static_cast<std::size_t>(size(a_region)))};
    }
}

inline bool contains(const Allele& lhs, const Allele& rhs)
{
    if (!contains(lhs.get_region(), rhs.get_region())) {
        return false;
    } else if (empty(lhs.get_region())) {
        // If the alleles are both insertions then both regions will be the same so we can only test
        // if the inserted rhs sequence is a subsequence of the lhs sequence. The rhs sequence
        // is required to be non-empty otherwise it would be a subsequence of everything.
        return !rhs.get_sequence().empty() && ::contains(lhs.get_sequence(), rhs.get_sequence());
    } else {
        return get_subsequence(lhs, rhs.get_region()) == rhs.get_sequence();
    }
}

#endif
