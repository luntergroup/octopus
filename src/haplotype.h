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
#include <algorithm> // std::equal
#include <iterator>  // std::cbegin etc
#include <stdexcept> // std::runtime_error
#include <ostream>

#include "variant.h"
#include "equitable.h"
#include "mappable.h"

class ReferenceGenome;
class GenomicRegion;

class Haplotype : public Equitable<Haplotype>, public Mappable<Haplotype>
{
public:
    using SequenceType = Variant::SequenceType;
    using SizeType     = Variant::SizeType;
    
    Haplotype() = delete;
    explicit Haplotype(ReferenceGenome& the_reference);
    explicit Haplotype(ReferenceGenome& the_reference, const GenomicRegion& the_region);
    ~Haplotype() = default;
    
    Haplotype(const Haplotype&)            = default;
    Haplotype& operator=(const Haplotype&) = default;
    Haplotype(Haplotype&&)                 = default;
    Haplotype& operator=(Haplotype&&)      = default;
    
    template <typename T>
    void emplace_back(const GenomicRegion& the_allele_region, T&& the_allele_sequence);
    template <typename T>
    void emplace_front(const GenomicRegion& the_allele_region, T&& the_allele_sequence);
    
    bool contains(const GenomicRegion& the_allele_region, const SequenceType& the_allele_sequence) const;
    GenomicRegion get_region() const;
    SequenceType get_sequence() const;
    SequenceType get_sequence(const GenomicRegion& a_region) const;
    
    //void operator+=(const Haplotype& other);
    friend bool operator==(const Haplotype& lhs, const Haplotype& rhs);
private:
    struct Allele : public Mappable<Allele>
    {
        template <typename T>
        Allele(const GenomicRegion& the_reference_region, T&& the_sequence)
        :
        the_reference_region {the_reference_region},
        the_sequence {std::forward<T>(the_sequence)}
        {}
        
        GenomicRegion the_reference_region;
        SequenceType the_sequence;
        
        GenomicRegion get_region() const { return the_reference_region; }
    };
    
    using AlleleIterator = std::deque<Allele>::const_iterator;
    
    friend bool operator==(const Allele& lhs, const Allele& rhs);
    friend bool operator<(const GenomicRegion& lhs, const Allele& rhs);
    friend bool operator<(const Allele& lhs, const GenomicRegion& rhs);
    
    GenomicRegion get_region_bounded_by_alleles() const;
    SequenceType get_sequence_bounded_by_alleles(AlleleIterator first, AlleleIterator last) const;
    SequenceType get_sequence_bounded_by_alleles() const;
    
    ReferenceGenome& the_reference_;
    bool is_region_set_;
    GenomicRegion the_reference_region_;
    std::deque<Allele> the_haplotype_;
};

void add_to_back(const Variant& a_variant, Haplotype& a_haplotype);
void add_to_front(const Variant& a_variant, Haplotype& a_haplotype);
bool contains(const Haplotype& a_haplotype, const Variant& a_variant);

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

inline bool operator==(const Haplotype::Allele& lhs, const Haplotype::Allele& rhs)
{
    return (lhs.the_reference_region == rhs.the_reference_region && lhs.the_sequence == rhs.the_sequence);
}

inline bool operator==(const Haplotype& lhs, const Haplotype& rhs)
{
    if (lhs.the_haplotype_.size() != rhs.the_haplotype_.size()) return false;
    if (lhs.get_region() != rhs.get_region()) return false;
    return std::equal(std::cbegin(lhs.the_haplotype_), std::cend(lhs.the_haplotype_),
                      std::cbegin(rhs.the_haplotype_));
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
