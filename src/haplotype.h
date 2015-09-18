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
#include <stdexcept> // std::runtime_error
#include <iostream>
#include <boost/functional/hash.hpp> // boost::hash_combine

#include "allele.h"
#include "variant.h"
#include "comparable.h"
#include "mappable.h"

class ReferenceGenome;
class GenomicRegion;

class Haplotype : public Comparable<Haplotype>, public Mappable<Haplotype>
{
public:
    using SequenceType = Allele::SequenceType;
    using SizeType     = Allele::SizeType;
    
    Haplotype() = default;
    explicit Haplotype(ReferenceGenome& reference);
    explicit Haplotype(ReferenceGenome& reference, const GenomicRegion& region);
    ~Haplotype() = default;
    
    Haplotype(const Haplotype&)            = default;
    Haplotype& operator=(const Haplotype&) = default;
    Haplotype(Haplotype&&)                 = default;
    Haplotype& operator=(Haplotype&&)      = default;
    
    explicit operator Allele() const;
    
    template <typename T> void push_back(T&& allele);
    template <typename T> void push_front(T&& allele);
    
    bool contains(const Allele& allele) const;
    void set_region(const GenomicRegion& region);
    GenomicRegion get_region() const;
    SequenceType get_sequence() const;
    SequenceType get_sequence(const GenomicRegion& region) const;
    
    unsigned num_transitions() const noexcept;
    unsigned num_transversions() const noexcept;
    
    void operator+=(const Haplotype& other);
    
    friend bool is_less_complex(const Haplotype& lhs, const Haplotype& rhs) noexcept;
    friend bool contains(const Haplotype& lhs, const Haplotype& rhs);
    friend Haplotype splice(const Haplotype& haplotype, const GenomicRegion& region);
    
    void print_explicit_alleles() const; // TEST
private:
    using AlleleIterator = std::deque<Allele>::const_iterator;
    
    GenomicRegion get_region_bounded_by_explicit_alleles() const;
    SequenceType get_sequence_bounded_by_explicit_alleles(AlleleIterator first, AlleleIterator last) const;
    SequenceType get_sequence_bounded_by_explicit_alleles() const;
    
    ReferenceGenome* reference_; // a non-owning pointer (rather than a reference) so Haplotype copyable
    GenomicRegion region_;
    mutable SequenceType cached_sequence_;
    std::deque<Allele> explicit_alleles_;
    bool is_region_set_;
    mutable bool is_cached_sequence_outdated_;
};

template <typename T>
void Haplotype::push_back(T&& allele)
{
    if (!explicit_alleles_.empty()) {
        if (!is_after(allele, explicit_alleles_.back())) {
            throw std::runtime_error {"cannot append out of order allele to back of haplotype"};
        } else if (!are_adjacent(explicit_alleles_.back(), allele)) {
            auto intervening_region = get_intervening(explicit_alleles_.back(), allele);
            explicit_alleles_.push_back(get_reference_allele(intervening_region, *reference_));
        }
        
        if (is_region_set_ && ends_before(region_, allele)) {
            region_ = get_encompassing(region_, allele);
        }
    } else if (is_region_set_ && begins_before(allele, region_)) {
        region_ = get_encompassing(allele, region_);
    }
    
    explicit_alleles_.push_back(std::forward<T>(allele));
    
    is_cached_sequence_outdated_ = true;
}

template <typename T>
void Haplotype::push_front(T&& allele)
{
    if (!explicit_alleles_.empty()) {
        if (!is_after(explicit_alleles_.front(), allele)) {
            throw std::runtime_error {"cannot append out of order allele to front of haplotype"};
        } else if (!are_adjacent(allele, explicit_alleles_.front())) {
            auto intervening_region = get_intervening(allele, explicit_alleles_.front());
            explicit_alleles_.push_front(get_reference_allele(intervening_region, *reference_));
        }
        
        if (is_region_set_ && begins_before(allele, region_)) {
            region_ = get_encompassing(allele, region_);
        }
    } else if (is_region_set_ && ends_before(region_, allele)) {
        region_ = get_encompassing(region_, allele);
    }
    
    explicit_alleles_.push_front(std::forward<T>(allele));
    
    is_cached_sequence_outdated_ = true;
}

// non-members

bool contains(const Haplotype& lhs, const Allele& rhs);
bool contains(const Haplotype& lhs, const Haplotype& rhs);

Haplotype splice(const Haplotype& haplotype, const GenomicRegion& region);

bool is_reference(const Haplotype& haplotype, ReferenceGenome& reference);

bool is_less_complex(const Haplotype& lhs, const Haplotype& rhs) noexcept;

void unique_least_complex(std::vector<Haplotype>& haplotypes);

bool operator==(const Haplotype& lhs, const Haplotype& rhs);
bool operator<(const Haplotype& lhs, const Haplotype& rhs);

namespace std
{
    template <> struct hash<Haplotype>
    {
        size_t operator()(const Haplotype& h) const
        {
            size_t seed {};
            boost::hash_combine(seed, hash<GenomicRegion>()(h.get_region()));
            boost::hash_combine(seed, hash<Haplotype::SequenceType>()(h.get_sequence()));
            return seed;
        }
    };
} // end namespace std

namespace boost
{
    template <> struct hash<Haplotype> : std::unary_function<Haplotype, std::size_t>
    {
        std::size_t operator()(const Haplotype& h) const
        {
            return std::hash<Haplotype>()(h);
        }
    };
} // end namespace boost

std::ostream& operator<<(std::ostream& os, const Haplotype& haplotype);

void add_to_back(const Variant& a_variant, Haplotype& haplotype);
void add_to_front(const Variant& a_variant, Haplotype& haplotype);
bool contains(const Haplotype& haplotype, const Variant& a_variant);

inline void Haplotype::print_explicit_alleles() const
{
    for (const auto& allele : explicit_alleles_) {
        std::cout << allele << std::endl;
    }
}

#endif /* defined(__Octopus__haplotype__) */
