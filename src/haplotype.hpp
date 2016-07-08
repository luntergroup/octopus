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
#include <iosfwd>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>

#include "contig_region.hpp"
#include "mappable.hpp"
#include "allele.hpp"
#include "comparable.hpp"

class ReferenceGenome;
class GenomicRegion;
class Variant;

class Haplotype;

namespace debug
{
    template <typename S> void print_alleles(S&& stream, const Haplotype& haplotype);
    template <typename S> void print_variant_alleles(S&& stream, const Haplotype& haplotype);
} // namespace debug

namespace detail
{
    Haplotype do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::true_type);
    Allele do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::false_type);
} // namespace detail

class Haplotype : public Comparable<Haplotype>, public Mappable<Haplotype>
{
public:
    using SequenceType = Allele::SequenceType;
    using SizeType     = Allele::SizeType;
    
    class Builder;
    
    Haplotype() = delete;
    
    template <typename R>
    Haplotype(R&& region, const ReferenceGenome& reference);
    
    template <typename R, typename S>
    Haplotype(R&& region, S&& sequence, const ReferenceGenome& reference);
    
    template <typename R, typename ForwardIt>
    Haplotype(R&& region, ForwardIt first_allele, ForwardIt last_allele,
              const ReferenceGenome& reference);
    
    ~Haplotype() = default;
    
    Haplotype(const Haplotype&)            = default;
    Haplotype& operator=(const Haplotype&) = default;
    Haplotype(Haplotype&&)                 = default;
    Haplotype& operator=(Haplotype&&)      = default;
    
    const GenomicRegion& mapped_region() const;
    
    bool contains(const ContigAllele& allele) const;
    bool contains(const Allele& allele) const;
    bool contains_exact(const Allele& allele) const;
    bool contains_exact(const ContigAllele& allele) const;
    
    SequenceType sequence(const ContigRegion& region) const;
    SequenceType sequence(const GenomicRegion& region) const;
    const SequenceType& sequence() const noexcept;
    
    SizeType sequence_size(const ContigRegion& region) const;
    SizeType sequence_size(const GenomicRegion& region) const;
    
    std::vector<Variant> difference(const Haplotype& other) const;
    
    std::size_t get_hash() const noexcept;
    
    friend struct HaveSameAlleles;
    friend struct IsLessComplex;
    
    friend bool contains(const Haplotype& lhs, const Haplotype& rhs);
    friend Haplotype detail::do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::true_type);
    friend bool is_reference(const Haplotype& haplotype);
    friend Haplotype expand(const Haplotype& haplotype, Haplotype::SizeType n);
    
    template <typename S> friend void debug::print_alleles(S&&, const Haplotype&);
    template <typename S> friend void debug::print_variant_alleles(S&&, const Haplotype&);
    
private:
    GenomicRegion region_;
    
    std::vector<ContigAllele> explicit_alleles_;
    
    ContigRegion explicit_allele_region_;
    
    SequenceType sequence_;
    
    std::size_t cached_hash_;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    using AlleleIterator = decltype(explicit_alleles_)::const_iterator;
    
    void append(SequenceType& result, const ContigAllele& allele) const;
    void append(SequenceType& result, AlleleIterator first, AlleleIterator last) const;
    void append_reference(SequenceType& result, const ContigRegion& region) const;
    SequenceType fetch_reference_sequence(const ContigRegion& region) const;
};

template <typename R>
Haplotype::Haplotype(R&& region, const ReferenceGenome& reference)
:
region_ {std::forward<R>(region)},
explicit_alleles_ {},
explicit_allele_region_ {},
sequence_ {reference.fetch_sequence(region_)},
cached_hash_ {std::hash<SequenceType>()(sequence_)},
reference_ {reference}
{}

template <typename R, typename S>
Haplotype::Haplotype(R&& region, S&& sequence, const ReferenceGenome& reference)
:
region_ {std::forward<R>(region)},
explicit_alleles_ {},
explicit_allele_region_ {region_.contig_region()},
sequence_ {std::forward<S>(sequence)},
cached_hash_ {std::hash<SequenceType>()(sequence_)},
reference_ {reference}
{
    explicit_alleles_.reserve(1);
    explicit_alleles_.emplace_back(explicit_allele_region_, sequence_);
}

namespace detail
{
    template <typename T>
    void append(T& result, const ReferenceGenome& reference,
                const GenomicRegion::ContigNameType& contig, const ContigRegion& region)
    {
        result.append(reference.fetch_sequence(GenomicRegion {contig, region}));
    }
}

template <typename R, typename ForwardIt>
Haplotype::Haplotype(R&& region, ForwardIt first_allele, ForwardIt last_allele,
                     const ReferenceGenome& reference)
:
region_ {std::forward<R>(region)},
explicit_alleles_ {first_allele, last_allele},
explicit_allele_region_ {},
sequence_ {},
cached_hash_ {0},
reference_ {reference}
{
    if (!explicit_alleles_.empty()) {
        explicit_allele_region_ = encompassing_region(explicit_alleles_.front(), explicit_alleles_.back());
        
        auto num_bases = std::accumulate(std::cbegin(explicit_alleles_),
                                         std::cend(explicit_alleles_), 0,
                                         [] (const auto curr, const auto& allele) {
                                             return curr + ::sequence_size(allele);
                                         });
        
        const auto lhs_reference_region = left_overhang_region(region_.contig_region(),
                                                               explicit_allele_region_);
        
        const auto rhs_reference_region = right_overhang_region(region_.contig_region(),
                                                                explicit_allele_region_);
        
        num_bases += region_size(lhs_reference_region) + region_size(rhs_reference_region);
        
        sequence_.reserve(num_bases);
        
        const auto& contig = region_.contig_name();
        
        if (!is_empty(lhs_reference_region)) {
            detail::append(sequence_, reference, contig, lhs_reference_region);
        }
        
        append(sequence_, std::cbegin(explicit_alleles_), std::cend(explicit_alleles_));
        
        if (!is_empty(rhs_reference_region)) {
            detail::append(sequence_, reference, contig, rhs_reference_region);
        }
    } else {
        sequence_ = reference.fetch_sequence(region_);
    }
    
    cached_hash_ = std::hash<SequenceType>()(sequence_);
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

ContigAllele splice(const Haplotype& haplotype, const ContigRegion& region);

template <typename MappableType, typename Container,
          typename = std::enable_if_t<std::is_same<typename Container::value_type, Haplotype>::value>>
std::vector<MappableType> splice_all(const Container& haplotypes, const GenomicRegion& region)
{
    std::vector<MappableType> result {};
    result.reserve(haplotypes.size());
    
    std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::back_inserter(result),
                   [&region] (const auto& haplotype) { return splice<MappableType>(haplotype, region); });
    
    std::sort(std::begin(result), std::end(result));
    
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    
    return result;
}

template <typename Container,
          typename = std::enable_if_t<std::is_same<typename Container::value_type, Haplotype>::value>>
std::vector<ContigAllele> splice_all(const Container& haplotypes, const ContigRegion& region)
{
    std::vector<ContigAllele> result {};
    result.reserve(haplotypes.size());
    
    for (const auto& haplotype : haplotypes) {
        result.emplace_back(splice(haplotype, region));
    }
    
    std::sort(std::begin(result), std::end(result));
    
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    
    return result;
}

bool is_reference(const Haplotype& haplotype);

Haplotype expand(const Haplotype& haplotype, Haplotype::SizeType n);

std::vector<Variant> difference(const Haplotype& lhs, const Haplotype& rhs);

bool operator==(const Haplotype& lhs, const Haplotype& rhs);
bool operator<(const Haplotype& lhs, const Haplotype& rhs);

struct HaveSameAlleles
{
    bool operator()(const Haplotype& lhs, const Haplotype& rhs) const;
};

struct IsLessComplex
{
    IsLessComplex() = default;
    explicit IsLessComplex(boost::optional<Haplotype> reference);
    bool operator()(const Haplotype& lhs, const Haplotype& rhs) const noexcept;
private:
    boost::optional<Haplotype> reference_;
};

unsigned unique_least_complex(std::vector<Haplotype>& haplotypes, boost::optional<Haplotype> = boost::none);

bool have_same_alleles(const Haplotype& lhs, const Haplotype& rhs);

bool are_equal_in_region(const Haplotype& lhs, const Haplotype& rhs, const GenomicRegion& region);

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
    template <> struct hash<Haplotype> : std::unary_function<Haplotype, std::size_t>
    {
        std::size_t operator()(const Haplotype& h) const
        {
            return std::hash<Haplotype>()(h);
        }
    };
} // namespace boost

std::ostream& operator<<(std::ostream& os, const Haplotype& haplotype);

namespace debug
{
    template <typename S>
    void print_alleles(S&& stream, const Haplotype& haplotype)
    {
        stream << "< ";
        for (const auto& allele : haplotype.explicit_alleles_) {
            stream << "{" << allele << "} ";
        }
        stream << ">";
    }
    
    template <typename S>
    void print_variant_alleles(S&& stream, const Haplotype& haplotype)
    {
        if (is_reference(haplotype)) {
            stream << "< >";
        } else {
            const auto& contig = contig_name(haplotype);
            stream << "< ";
            for (const auto& contig_allele : haplotype.explicit_alleles_) {
                Allele allele {GenomicRegion {contig, contig_allele.mapped_region()}, contig_allele.sequence()};
                if (!is_reference(allele, haplotype.reference_)) stream << "{" << allele << "} ";
            }
            stream << ">";
        }
    }
    
    void print_alleles(const Haplotype& haplotype);
    void print_variant_alleles(const Haplotype& haplotype);
    
    Haplotype make_haplotype(const std::string& str, const GenomicRegion& region,
                             const ReferenceGenome& reference);
    Haplotype make_haplotype(const std::string& str, const std::string& region,
                             const ReferenceGenome& reference);
} // namespace debug

#endif /* defined(__Octopus__haplotype__) */
