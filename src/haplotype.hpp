//
//  haplotype.hpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__haplotype__
#define __Octopus__haplotype__

#include <deque>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>
#include <numeric>
#include <ostream>

#include <boost/functional/hash.hpp>

#include "contig_region.hpp"
#include "mappable.hpp"
#include "allele.hpp"
#include "comparable.hpp"

class ReferenceGenome;
class GenomicRegion;
class Variant;

class Haplotype;

namespace detail
{
    Haplotype do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::true_type);
    Allele do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::false_type);
}

class Haplotype : public Comparable<Haplotype>, public Mappable<Haplotype>
{
public:
    using SequenceType = Allele::SequenceType;
    using SizeType     = Allele::SizeType;
    
    class Builder;
    
    Haplotype() = default;
    
    template <typename R>
    explicit Haplotype(R&& region, const ReferenceGenome& reference);
    
    template <typename R, typename ForwardIt>
    explicit Haplotype(R&& region, ForwardIt first_allele, ForwardIt last_allele,
                       const ReferenceGenome& reference);
    
    ~Haplotype() = default;
    
    Haplotype(const Haplotype&)            = default;
    Haplotype& operator=(const Haplotype&) = default;
    Haplotype(Haplotype&&)                 = default;
    Haplotype& operator=(Haplotype&&)      = default;
    
    const GenomicRegion& get_region() const;
    
    bool contains(const ContigAllele& allele) const;
    bool contains_exact(const ContigAllele& allele) const;
    bool contains(const Allele& allele) const;
    bool contains_exact(const Allele& allele) const;
    
    SequenceType get_sequence(const ContigRegion& region) const;
    SequenceType get_sequence(const GenomicRegion& region) const;
    const SequenceType& get_sequence() const noexcept;
    
    std::vector<Variant> difference(const Haplotype& other) const;
    
    size_t get_hash() const noexcept;
    
    friend struct HaveSameAlleles;
    
    friend bool contains(const Haplotype& lhs, const Haplotype& rhs);
    friend Haplotype detail::do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::true_type);
    friend bool is_reference(const Haplotype& haplotype);
    
    friend void print_alleles(const Haplotype& haplotype);
    friend void print_variant_alleles(const Haplotype& haplotype);
    
private:
    GenomicRegion region_;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    std::vector<ContigAllele> explicit_alleles_;
    
    ContigRegion explicit_allele_region_;
    
    SequenceType cached_sequence_;
    size_t cached_hash_;
    
    using AlleleIterator = decltype(explicit_alleles_)::const_iterator;
    
    SequenceType get_reference_sequence(const GenomicRegion& region) const;
    SequenceType get_reference_sequence(const ContigRegion& region) const;
    SequenceType get_sequence_bounded_by_explicit_alleles(AlleleIterator first, AlleleIterator last) const;
    SequenceType get_sequence_bounded_by_explicit_alleles() const;
};

template <typename R>
Haplotype::Haplotype(R&& region, const ReferenceGenome& reference)
:
region_ {std::forward<R>(region)},
reference_ {reference},
explicit_alleles_ {},
explicit_allele_region_ {},
cached_sequence_ {get_reference_sequence(region_)},
cached_hash_ {std::hash<SequenceType>()(cached_sequence_)}
{}

template <typename R, typename ForwardIt>
Haplotype::Haplotype(R&& region, ForwardIt first_allele, ForwardIt last_allele,
                     const ReferenceGenome& reference)
:
region_ {std::forward<R>(region)},
reference_ {reference},
explicit_alleles_ {first_allele, last_allele},
explicit_allele_region_ {},
cached_sequence_ {},
cached_hash_ {0}
{
    if (!explicit_alleles_.empty()) {
        explicit_allele_region_ = encompassing_region(explicit_alleles_.front(), explicit_alleles_.back());
        
        auto num_bases = std::accumulate(std::cbegin(explicit_alleles_),
                                         std::cend(explicit_alleles_), 0,
                                         [] (const auto curr, const auto& allele) {
                                             return curr + sequence_size(allele);
                                         });
        
        const auto lhs_reference_region = left_overhang_region(region_.get_contig_region(),
                                                               explicit_allele_region_);
        
        const auto rhs_reference_region = right_overhang_region(region_.get_contig_region(),
                                                                explicit_allele_region_);
        
        num_bases += region_size(lhs_reference_region) + region_size(rhs_reference_region);
        
        cached_sequence_.reserve(num_bases);
        
        if (!is_empty_region(lhs_reference_region)) {
            cached_sequence_ += get_reference_sequence(lhs_reference_region);
        }
        
        cached_sequence_ += get_sequence_bounded_by_explicit_alleles();
        
        if (!is_empty_region(rhs_reference_region)) {
            cached_sequence_ += get_reference_sequence(rhs_reference_region);
        }
    } else {
        cached_sequence_ = get_reference_sequence(region_);
    }
    
    cached_hash_ = std::hash<SequenceType>()(cached_sequence_);
}

class Haplotype::Builder
{
public:
    Builder() = default;
    explicit Builder(const GenomicRegion& region, const ReferenceGenome& reference);
    ~Builder() = default;
    
    Builder(const Builder&)            = default;
    Builder& operator=(const Builder&) = default;
    Builder(Builder&&)                 = default;
    Builder& operator=(Builder&&)      = default;
    
    void push_back(const ContigAllele& allele);
    void push_front(const ContigAllele& allele);
    void push_back(ContigAllele&& allele);
    void push_front(ContigAllele&& allele);
    void push_back(const Allele& allele);
    void push_front(const Allele& allele);
    void push_back(Allele&& allele);
    void push_front(Allele&& allele);
    
    Haplotype build();
    
    friend Haplotype detail::do_splice(const Haplotype& haplotype, const GenomicRegion& region,
                                       std::true_type);
    
private:
    GenomicRegion region_;
    
    std::deque<ContigAllele> explicit_alleles_;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    ContigAllele get_intervening_reference_allele(const ContigAllele& lhs, const ContigAllele& rhs) const;
    void update_region(const ContigAllele& allele) noexcept;
    void update_region(const Allele& allele);
};

// non-members

Haplotype::SizeType sequence_size(const Haplotype& haplotype) noexcept;

bool is_empty_sequence(const Haplotype& haplotype) noexcept;

bool contains(const Haplotype& lhs, const Allele& rhs);
bool contains(const Haplotype& lhs, const Haplotype& rhs);

template <typename MappableType>
MappableType splice(const Haplotype& haplotype, const GenomicRegion& region)
{
    return detail::do_splice(haplotype, region, std::is_same<Haplotype, std::decay_t<MappableType>> {});
}

bool is_reference(const Haplotype& haplotype);

bool operator==(const Haplotype& lhs, const Haplotype& rhs);
bool operator<(const Haplotype& lhs, const Haplotype& rhs);

struct HaveSameAlleles
{
    bool operator()(const Haplotype& lhs, const Haplotype& rhs) const;
};

bool have_same_alleles(const Haplotype& lhs, const Haplotype& rhs);

bool are_equal_in_region(const Haplotype& lhs, const Haplotype& rhs, const GenomicRegion& region);

void add_ref_to_back(const Variant& variant, Haplotype& haplotype);
void add_ref_to_front(const Variant& variant, Haplotype& haplotype);
void add_alt_to_back(const Variant& variant, Haplotype& haplotype);
void add_alt_to_front(const Variant& variant, Haplotype& haplotype);

template <typename MappableType, typename Container,
          typename = std::enable_if_t<std::is_same<typename Container::value_type, Haplotype>::value>>
std::deque<MappableType> splice_all(const Container& haplotypes, const GenomicRegion& region)
{
    std::deque<MappableType> result {};
    
    for (const auto& haplotype : haplotypes) {
        result.push_back(splice<MappableType>(haplotype, region));
    }
    
    return result;
}

namespace std
{
    template <> struct hash<Haplotype>
    {
        size_t operator()(const Haplotype& haplotype) const
        {
            return haplotype.get_hash();
        }
    };
    
    template <> struct hash<reference_wrapper<const Haplotype>>
    {
        size_t operator()(const reference_wrapper<const Haplotype> haplotype) const
        {
            return hash<Haplotype>()(haplotype);
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

void print_alleles(const Haplotype& haplotype);
void print_variant_alleles(const Haplotype& haplotype);

#endif /* defined(__Octopus__haplotype__) */
