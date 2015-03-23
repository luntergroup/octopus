//
//  haplotype.h
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__haplotype__
#define __Octopus__haplotype__

#include <queue>
#include <ostream>

#include "variant.h"
#include "equitable.h"

class ReferenceGenome;
class GenomicRegion;

class Haplotype : public Equitable<Haplotype>
{
public:
    using SequenceType = Variant::SequenceType;
    
    Haplotype() = delete;
    explicit Haplotype(ReferenceGenome& the_reference);
    ~Haplotype() = default;
    
    Haplotype(const Haplotype&)            = default;
    Haplotype& operator=(const Haplotype&) = default;
    Haplotype(Haplotype&&)                 = default;
    Haplotype& operator=(Haplotype&&)      = default;
    
    GenomicRegion get_region() const;
    SequenceType get_sequence() const;
    
    void emplace_back(const Variant& a_variant);
    template <typename T>
    void emplace_back(const GenomicRegion& the_allele_region, T&& the_allele_sequence);
    void emplace_front(const Variant& a_variant);
    template <typename T>
    void emplace_front(const GenomicRegion& the_allele_region, T&& the_allele_sequence);
    
private:
    struct Allele
    {
        template <typename T>
        Allele(const GenomicRegion& the_reference_region, T&& the_sequence)
        :
        the_reference_region {the_reference_region},
        the_sequence {std::forward<T>(the_sequence)}
        {}
        
        GenomicRegion the_reference_region;
        SequenceType the_sequence;
    };
    
    ReferenceGenome& the_reference_;
    std::deque<Allele> the_haplotype_;
};

template <typename T>
void Haplotype::emplace_back(const GenomicRegion& the_allele_region, T&& the_allele_sequence)
{
    the_haplotype_.emplace_back(the_allele_region, std::forward<T>(the_allele_sequence));
}

template <typename T>
void Haplotype::emplace_front(const GenomicRegion& the_allele_region, T&& the_allele_sequence)
{
    the_haplotype_.emplace_front(the_allele_region, std::forward<T>(the_allele_sequence));
}

inline bool operator==(const Haplotype& lhs, const Haplotype& rhs)
{
    return lhs.get_region() == rhs.get_region();
}

namespace std {
    template <> struct hash<Haplotype>
    {
        size_t operator()(const Haplotype& h) const
        {
            return hash<string>()(to_string(h.get_region()));
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const Haplotype& a_haplotype)
{
    os << a_haplotype.get_sequence();
    return os;
}

#endif /* defined(__Octopus__haplotype__) */
