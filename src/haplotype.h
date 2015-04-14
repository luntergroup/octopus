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

#include "allele.h"
#include "variant.h"
#include "equitable.h"
#include "mappable.h"

class ReferenceGenome;
class GenomicRegion;

class Haplotype : public Equitable<Haplotype>, public Mappable<Haplotype>
{
public:
    using SequenceType = Allele::SequenceType;
    using SizeType     = Allele::SizeType;
    
    Haplotype() = delete;
    explicit Haplotype(ReferenceGenome& the_reference);
    explicit Haplotype(ReferenceGenome& the_reference, const GenomicRegion& the_region);
    ~Haplotype() = default;
    
    Haplotype(const Haplotype&)            = default;
    Haplotype& operator=(const Haplotype&) = default;
    Haplotype(Haplotype&&)                 = default;
    Haplotype& operator=(Haplotype&&)      = default;
    
    template <typename T> void push_back(T&& an_allele);
    template <typename T> void push_front(T&& an_allele);
    
    bool contains(const Allele& an_allele) const;
    void set_region(const GenomicRegion& a_region);
    GenomicRegion get_region() const;
    SequenceType get_sequence() const;
    SequenceType get_sequence(const GenomicRegion& a_region) const;
    
    //void operator+=(const Haplotype& other);
private:
    using AlleleIterator = std::deque<Allele>::const_iterator;
    
    GenomicRegion get_region_bounded_by_explicit_alleles() const;
    SequenceType get_sequence_bounded_by_explicit_alleles(AlleleIterator first, AlleleIterator last) const;
    SequenceType get_sequence_bounded_by_explicit_alleles() const;
    
    ReferenceGenome& the_reference_;
    bool is_region_set_;
    GenomicRegion the_reference_region_;
    std::deque<Allele> the_explicit_alleles_;
    
    mutable SequenceType cached_sequence_;
    mutable bool is_cached_sequence_outdated_;
};

template <typename T>
void Haplotype::push_back(T&& an_allele)
{
    if (!the_explicit_alleles_.empty()) {
        if (!is_after(an_allele, the_explicit_alleles_.back())) {
            throw std::runtime_error {"Cannot append out of order allele to back of haplotype"};
        } else if (!are_adjacent(the_explicit_alleles_.back(), an_allele)) {
            auto intervening_region = get_intervening_region(the_explicit_alleles_.back(), an_allele);
            the_explicit_alleles_.push_back(get_reference_allele(intervening_region, the_reference_));
        }
    }
    
    the_explicit_alleles_.push_back(std::forward<T>(an_allele));
    is_cached_sequence_outdated_ = true;
}

template <typename T>
void Haplotype::push_front(T&& an_allele)
{
    if (!the_explicit_alleles_.empty()) {
        if (!is_before(an_allele, the_explicit_alleles_.front())) {
            throw std::runtime_error {"Cannot append out of order allele to front of haplotype"};
        } else if (!are_adjacent(an_allele, the_explicit_alleles_.front())) {
            auto intervening_region = get_intervening_region(the_explicit_alleles_.back(), an_allele);
            the_explicit_alleles_.push_front(get_reference_allele(intervening_region, the_reference_));
        }
    }
    
    the_explicit_alleles_.push_front(std::forward<T>(an_allele));
    is_cached_sequence_outdated_ = true;
}

inline bool operator==(const Haplotype& lhs, const Haplotype& rhs)
{
    return lhs.get_region() == rhs.get_region() && lhs.get_sequence() == rhs.get_sequence();
}

namespace std {
    template <> struct hash<Haplotype>
    {
        size_t operator()(const Haplotype& h) const
        {
            return hash<string>()(to_string(h.get_region())); //TODO: see if this can be improved
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const Haplotype& a_haplotype)
{
    os << a_haplotype.get_region() << " " << a_haplotype.get_sequence();
    return os;
}

void add_to_back(const Variant& a_variant, Haplotype& a_haplotype);
void add_to_front(const Variant& a_variant, Haplotype& a_haplotype);
bool contains(const Haplotype& a_haplotype, const Variant& a_variant);

#endif /* defined(__Octopus__haplotype__) */
