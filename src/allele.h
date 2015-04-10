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
#include <cstdint>
#include <ostream>

#include "genomic_region.h"
#include "comparable.h"
#include "mappable.h"
#include "reference_genome.h"

using std::uint_fast32_t;

class Allele : public Comparable<Allele>, public Mappable<Allele>
{
public:
    using SizeType     = GenomicRegion::SizeType;
    using SequenceType = std::string;
    
    Allele() = default;
    template <typename T>
    Allele(const GenomicRegion& the_reference_region, T&& the_sequence);
    ~Allele() = default;
    
    const GenomicRegion& get_region() const noexcept;
    const SequenceType& get_sequence() const noexcept;
    
private:
    GenomicRegion the_reference_region_;
    SequenceType the_sequence_;
};

template <typename T>
Allele::Allele(const GenomicRegion& the_reference_region, T&& the_sequence)
:
the_reference_region_ {the_reference_region},
the_sequence_ {std::forward<T>(the_sequence)}
{}

const GenomicRegion& Allele::get_region() const noexcept
{
    return the_reference_region_;
}

const Allele::SequenceType& Allele::get_sequence() const noexcept
{
    return the_sequence_;
}

inline bool operator==(const Allele& lhs, const Allele& rhs)
{
    return lhs.get_region() == rhs.get_region() && lhs.get_sequence() == rhs.get_sequence();
}

inline bool operator<(const Allele& lhs, const Allele& rhs)
{
    return (lhs.get_region() == rhs.get_region()) ? lhs.get_sequence() == rhs.get_sequence() :
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

#endif
