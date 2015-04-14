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
    
    void push_back(const Allele& an_allele);
    void push_front(const Allele& an_allele);
    
    template <typename T>
    void emplace_back(const GenomicRegion& the_allele_region, T&& the_allele_sequence);
    template <typename T>
    void emplace_front(const GenomicRegion& the_allele_region, T&& the_allele_sequence);
    
    bool contains(const Allele& an_allele) const;
    bool contains(const GenomicRegion& the_allele_region, const SequenceType& the_allele_sequence) const;
    void set_region(const GenomicRegion& a_region);
    GenomicRegion get_region() const;
    SequenceType get_sequence() const;
    SequenceType get_sequence(const GenomicRegion& a_region) const;
    
    //void operator+=(const Haplotype& other);
    friend bool operator==(const Haplotype& lhs, const Haplotype& rhs);
private:
    using AlleleIterator = std::deque<Allele>::const_iterator;
    
    GenomicRegion get_region_bounded_by_explicit_alleles() const;
    SequenceType get_sequence_bounded_by_explicit_alleles(AlleleIterator first, AlleleIterator last) const;
    SequenceType get_sequence_bounded_by_explicit_alleles() const;
    
    ReferenceGenome& the_reference_;
    bool is_region_set_;
    GenomicRegion the_reference_region_;
    std::deque<Allele> the_explicit_alleles_;
};

template <typename T>
void Haplotype::emplace_back(const GenomicRegion& the_allele_region, T&& the_allele_sequence)
{
    the_explicit_alleles_.emplace_back(the_allele_region, std::forward<T>(the_allele_sequence));
}

template <typename T>
void Haplotype::emplace_front(const GenomicRegion& the_allele_region, T&& the_allele_sequence)
{
    the_explicit_alleles_.emplace_front(the_allele_region, std::forward<T>(the_allele_sequence));
}

inline bool operator==(const Haplotype& lhs, const Haplotype& rhs)
{
    //return (lhs.get_region() == rhs.get_region() && lhs.get_sequence() == rhs.get_sequence());
    
//    if (lhs.get_region() != rhs.get_region()) return false;
//    auto p = std::minmax(lhs, rhs, [] (const auto& lhs, const auto& rhs) {
//        return lhs.the_explicit_alleles_.size() < rhs.the_explicit_alleles_.size();
//    });
//    return std::all_of(std::cbegin(p.second.the_explicit_alleles_), std::cend(p.second.the_explicit_alleles_),
//                       [&p] (const auto& allele) {
//                           return p.first.contains(allele);
//                       });
    
    if (lhs.the_explicit_alleles_.size() != rhs.the_explicit_alleles_.size()) return false;
    if (lhs.get_region() != rhs.get_region()) return false;
    return std::equal(std::cbegin(lhs.the_explicit_alleles_), std::cend(lhs.the_explicit_alleles_),
                      std::cbegin(rhs.the_explicit_alleles_));
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

void add_to_back(const Variant& a_variant, Haplotype& a_haplotype);
void add_to_front(const Variant& a_variant, Haplotype& a_haplotype);
bool contains(const Haplotype& a_haplotype, const Variant& a_variant);

#endif /* defined(__Octopus__haplotype__) */
