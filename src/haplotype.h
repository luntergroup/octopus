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
#include <ostream>
#include <boost/functional/hash.hpp> // boost::hash_combine

#include "allele.h"
#include "variant.h"
#include "comparable.h"
#include "mappable.h"

#include <iostream> // TEST

class ReferenceGenome;
class GenomicRegion;

class Haplotype : public Comparable<Haplotype>, public Mappable<Haplotype>
{
public:
    using SequenceType = Allele::SequenceType;
    using SizeType     = Allele::SizeType;
    
    Haplotype() = default;
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
    
    ReferenceGenome* the_reference_; // a non-owning pointer (rather than a reference) so Haplotype copyable
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
            std::cout << the_explicit_alleles_.back() << std::endl;
            std::cout << an_allele << std::endl;
            throw std::runtime_error {"cannot append out of order allele to back of haplotype"};
        } else if (!are_adjacent(the_explicit_alleles_.back(), an_allele)) {
            auto intervening_region = get_intervening(the_explicit_alleles_.back(), an_allele);
            the_explicit_alleles_.push_back(get_reference_allele(intervening_region, *the_reference_));
        }
        
        if (is_region_set_ && ends_before(the_reference_region_, an_allele)) {
            the_reference_region_ = get_encompassing(the_reference_region_, an_allele);
        }
    } else if (is_region_set_ && begins_before(an_allele, the_reference_region_)) {
        the_reference_region_ = get_encompassing(an_allele, the_reference_region_);
    }
    
    the_explicit_alleles_.push_back(std::forward<T>(an_allele));
    
    is_cached_sequence_outdated_ = true;
}

template <typename T>
void Haplotype::push_front(T&& an_allele)
{
    if (!the_explicit_alleles_.empty()) {
        if (!is_after(the_explicit_alleles_.front(), an_allele)) {
            throw std::runtime_error {"cannot append out of order allele to front of haplotype"};
        } else if (!are_adjacent(an_allele, the_explicit_alleles_.front())) {
            auto intervening_region = get_intervening(an_allele, the_explicit_alleles_.front());
            the_explicit_alleles_.push_front(get_reference_allele(intervening_region, *the_reference_));
        }
        
        if (is_region_set_ && begins_before(an_allele, the_reference_region_)) {
            the_reference_region_ = get_encompassing(an_allele, the_reference_region_);
        }
    } else if (is_region_set_ && ends_before(the_reference_region_, an_allele)) {
        the_reference_region_ = get_encompassing(the_reference_region_, an_allele);
    }
    
    the_explicit_alleles_.push_front(std::forward<T>(an_allele));
    
    is_cached_sequence_outdated_ = true;
}

// non-members

bool contains(const Haplotype& lhs, const Haplotype& rhs);

Haplotype splice(const Haplotype& haplotype, const GenomicRegion& region);

bool is_reference(const Haplotype& a_haplotype, ReferenceGenome& the_reference);

bool is_less_complex(const Haplotype& lhs, const Haplotype& rhs) noexcept;

void unique_least_complex(std::vector<Haplotype>& haplotypes);

inline bool operator==(const Haplotype& lhs, const Haplotype& rhs)
{
    return lhs.get_region() == rhs.get_region() && lhs.get_sequence() == rhs.get_sequence();
}

inline bool operator<(const Haplotype& lhs, const Haplotype& rhs)
{
    return (lhs.get_region() == rhs.get_region()) ? lhs.get_sequence() < rhs.get_sequence() :
                                                    lhs.get_region() < rhs.get_region();
}

namespace std {
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
}

namespace boost {
    template <> struct hash<Haplotype> : std::unary_function<Haplotype, std::size_t>
    {
        std::size_t operator()(const Haplotype& h) const
        {
            return std::hash<Haplotype>()(h);
        }
    };
}

std::ostream& operator<<(std::ostream& os, const Haplotype& a_haplotype);

void add_to_back(const Variant& a_variant, Haplotype& a_haplotype);
void add_to_front(const Variant& a_variant, Haplotype& a_haplotype);
bool contains(const Haplotype& a_haplotype, const Variant& a_variant);

inline void Haplotype::print_explicit_alleles() const
{
    for (const auto& allele : the_explicit_alleles_) {
        std::cout << allele << std::endl;
    }
}

#endif /* defined(__Octopus__haplotype__) */
