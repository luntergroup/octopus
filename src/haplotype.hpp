//
//  haplotype.hpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__haplotype__
#define __Octopus__haplotype__

#include <queue>
#include <cstddef>   // size_t
#include <stdexcept> // std::runtime_error
#include <iostream>

#include "contig_region.hpp"
#include "mappable.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "comparable.hpp"
#include "mappable.hpp"

class ReferenceGenome;
class GenomicRegion;

class Haplotype;

namespace detail {
    Haplotype do_splice(const Haplotype& haplotype, const GenomicRegion& region, Haplotype);
    Allele do_splice(const Haplotype& haplotype, const GenomicRegion& region, Allele);
} // namespace detail

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
    
    std::vector<Variant> difference(const Haplotype& from) const;
    unsigned num_transitions() const noexcept;
    unsigned num_transversions() const noexcept;
    
    void operator+=(const Haplotype& other);
    
    size_t get_hash() const;
    
    friend struct IsLessComplex;
    
    friend bool contains(const Haplotype& lhs, const Haplotype& rhs);
    friend Haplotype detail::do_splice(const Haplotype& haplotype, const GenomicRegion& region, Haplotype);
    friend bool have_same_alleles(const Haplotype& lhs, const Haplotype& rhs);
    
    friend void print_alleles(const Haplotype& haplotype);
    friend void print_variant_alleles(const Haplotype& haplotype);
    
private:
    // TODO: store these instead of Allele's... don't need to keep the Contig name in each allele.. it will
    // save some space
    class ContigAllele : Mappable<ContigAllele>
    {
    public:
        using SequenceType = Allele::SequenceType;
        
        ContigAllele() = default;
        template <typename R, typename S> ContigAllele(R&& region, S&& sequence);
        //template <typename S> ContigAllele(const GenomicRegion& region, S&& sequence);
        
        const ContigRegion& get_region() const;
        const SequenceType& get_sequence() const;
        
    private:
        const ContigRegion region_;
        const SequenceType sequence_;
    };
    
    ReferenceGenome* reference_; // non-owning pointer rather than a reference so Haplotype copyable
    GenomicRegion region_;
    mutable SequenceType cached_sequence_;
    std::deque<Allele> explicit_alleles_;
    mutable size_t cached_hash_;
    bool is_region_set_;
    mutable bool is_cached_sequence_outdated_;
    
    using AlleleIterator = decltype(explicit_alleles_)::const_iterator;
    
    GenomicRegion get_region_bounded_by_explicit_alleles() const;
    SequenceType get_sequence_bounded_by_explicit_alleles(AlleleIterator first, AlleleIterator last) const;
    SequenceType get_sequence_bounded_by_explicit_alleles() const;
};

template <typename R, typename S>
Haplotype::ContigAllele::ContigAllele(R&& region, S&& sequence)
:
region_ {std::forward<R>(region)},
sequence_ {std::forward<S>(sequence)}
{}

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

template <typename MappableType>
MappableType splice(const Haplotype& haplotype, const GenomicRegion& region)
{
    return detail::do_splice(haplotype, region, MappableType());
}

bool is_reference(const Haplotype& haplotype, ReferenceGenome& reference);

struct IsLessComplex
{
    bool operator()(const Haplotype& lhs, const Haplotype& rhs) noexcept;
};

void unique_least_complex(std::vector<Haplotype>& haplotypes);

bool operator==(const Haplotype& lhs, const Haplotype& rhs);
bool operator<(const Haplotype& lhs, const Haplotype& rhs);

bool have_same_alleles(const Haplotype& lhs, const Haplotype& rhs);

namespace std
{
    template <> struct hash<Haplotype>
    {
        size_t operator()(const Haplotype& h) const
        {
            return h.get_hash();
        }
    };
} // namespace std

namespace boost
{
    template <> struct hash<Haplotype> : std::unary_function<Haplotype, size_t>
    {
        size_t operator()(const Haplotype& h) const
        {
            return std::hash<Haplotype>()(h);
        }
    };
} // namespace boost

std::ostream& operator<<(std::ostream& os, const Haplotype& haplotype);

void add_to_back(const Variant& a_variant, Haplotype& haplotype);
void add_to_front(const Variant& a_variant, Haplotype& haplotype);
bool contains(const Haplotype& haplotype, const Variant& a_variant);

void print_alleles(const Haplotype& haplotype);
void print_variant_alleles(const Haplotype& haplotype);

#endif /* defined(__Octopus__haplotype__) */
