// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef haplotype_hpp
#define haplotype_hpp

#include <deque>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>
#include <numeric>
#include <iosfwd>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>

#include "concepts/comparable.hpp"
#include "concepts/mappable.hpp"
#include "basics/contig_region.hpp"
#include "basics/cigar_string.hpp"
#include "allele.hpp"
#include "variant.hpp"

namespace octopus {

class GenomicRegion;
class ReferenceGenome;

/*
    A Haplotype is an ordered, non-overlapping, set of Alleles, and therefore implictly
    defines a sequence in a given GenomicRegion.
 */
class Haplotype;

namespace debug {

template <typename S> void print_alleles(S&& stream, const Haplotype& haplotype);
template <typename S> void print_variant_alleles(S&& stream, const Haplotype& haplotype);

} // namespace debug

namespace detail {

Haplotype do_copy(const Haplotype& haplotype, const GenomicRegion& region, std::true_type);
Allele do_copy(const Haplotype& haplotype, const GenomicRegion& region, std::false_type);

} // namespace detail

class Haplotype : public Comparable<Haplotype>, public Mappable<Haplotype>
{
public:
    using MappingDomain      = Allele::MappingDomain;
    using NucleotideSequence = Allele::NucleotideSequence;
    
    class Builder;
    
    Haplotype() = delete;
    
    template <typename R>
    Haplotype(R&& region, const ReferenceGenome& reference);
    
    template <typename R, typename S>
    Haplotype(R&& region, S&& sequence, const ReferenceGenome& reference);
    
    template <typename R, typename ForwardIt>
    Haplotype(R&& region, ForwardIt first_allele, ForwardIt last_allele,
              const ReferenceGenome& reference);
    
    Haplotype(const Haplotype&)            = default;
    Haplotype& operator=(const Haplotype&) = default;
    Haplotype(Haplotype&&)                 = default;
    Haplotype& operator=(Haplotype&&)      = default;
    
    ~Haplotype() = default;
    
    const GenomicRegion& mapped_region() const;
    
    bool contains(const ContigAllele& allele) const;
    bool contains(const Allele& allele) const;
    
    bool includes(const ContigAllele& allele) const;
    bool includes(const Allele& allele) const;
    
    NucleotideSequence sequence(const ContigRegion& region) const;
    NucleotideSequence sequence(const GenomicRegion& region) const;
    const NucleotideSequence& sequence() const noexcept;
    
    NucleotideSequence::size_type sequence_size(const ContigRegion& region) const;
    NucleotideSequence::size_type sequence_size(const GenomicRegion& region) const;
    
    std::vector<Variant> difference(const Haplotype& other) const; // w.r.t this
    CigarString cigar() const; // w.r.t reference
    
    std::size_t get_hash() const noexcept;
    
    friend struct StrictLess;
    friend struct HaveSameAlleles;
    friend struct IsLessComplex;
    
    friend bool contains(const Haplotype& lhs, const Haplotype& rhs);
    friend Haplotype detail::do_copy(const Haplotype& haplotype, const GenomicRegion& region, std::true_type);
    friend Haplotype copy(const Haplotype&, const std::vector<GenomicRegion>&);
    friend bool is_reference(const Haplotype& haplotype);
    friend Haplotype expand(const Haplotype& haplotype, MappingDomain::Position n);
    friend Haplotype remap(const Haplotype& haplotype, const GenomicRegion& region);
    
    template <typename S> friend void debug::print_alleles(S&&, const Haplotype&);
    template <typename S> friend void debug::print_variant_alleles(S&&, const Haplotype&);
    
private:
    GenomicRegion region_;
    std::vector<ContigAllele> explicit_alleles_;
    ContigRegion explicit_allele_region_;
    NucleotideSequence sequence_;
    std::size_t cached_hash_;
    std::reference_wrapper<const ReferenceGenome> reference_;

public:
    using AlleleIterator = decltype(explicit_alleles_)::const_iterator;
    
    std::pair<AlleleIterator, AlleleIterator> alleles() const noexcept;

private:
    
    void append(NucleotideSequence& result, const ContigAllele& allele) const;
    void append(NucleotideSequence& result, AlleleIterator first, AlleleIterator last) const;
    void append_reference(NucleotideSequence& result, const ContigRegion& region) const;
    NucleotideSequence fetch_reference_sequence(const ContigRegion& region) const;
};

template <typename R>
Haplotype::Haplotype(R&& region, const ReferenceGenome& reference)
: region_ {std::forward<R>(region)}
, explicit_alleles_ {}
, explicit_allele_region_ {}
, sequence_ {reference.fetch_sequence(region_)}
, cached_hash_ {std::hash<NucleotideSequence>()(sequence_)}
, reference_ {reference}
{}

template <typename R, typename S>
Haplotype::Haplotype(R&& region, S&& sequence, const ReferenceGenome& reference)
: region_ {std::forward<R>(region)}
, explicit_alleles_ {}
, explicit_allele_region_ {region_.contig_region()}
, sequence_ {std::forward<S>(sequence)}
, cached_hash_ {std::hash<NucleotideSequence>()(sequence_)}
, reference_ {reference}
{
    explicit_alleles_.reserve(1);
    explicit_alleles_.emplace_back(explicit_allele_region_, sequence_);
}

namespace detail {
    template <typename T>
    void append(T& result, const ReferenceGenome& reference,
                const GenomicRegion::ContigName& contig,
                const ContigRegion& region)
    {
        result.append(reference.fetch_sequence(GenomicRegion {contig, region}));
    }
}

template <typename R, typename ForwardIt>
Haplotype::Haplotype(R&& region, ForwardIt first_allele, ForwardIt last_allele,
                     const ReferenceGenome& reference)
: region_ {std::forward<R>(region)}
, explicit_alleles_ {first_allele, last_allele}
, explicit_allele_region_ {}
, sequence_ {}
, cached_hash_ {0}
, reference_ {reference}
{
    if (!explicit_alleles_.empty()) {
        explicit_allele_region_ = encompassing_region(explicit_alleles_.front(), explicit_alleles_.back());
        auto num_bases = std::accumulate(std::cbegin(explicit_alleles_), std::cend(explicit_alleles_),
                                         0, [] (const auto curr, const auto& allele) {
                                             return curr + ::octopus::sequence_size(allele);
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
    cached_hash_ = std::hash<NucleotideSequence>()(sequence_);
}

class Haplotype::Builder
{
public:
    Builder() = delete;
    
    explicit Builder(const GenomicRegion& region, const ReferenceGenome& reference);
    
    Builder(const Builder&)            = default;
    Builder& operator=(const Builder&) = default;
    Builder(Builder&&)                 = default;
    Builder& operator=(Builder&&)      = default;
    
    ~Builder() = default;
    
    bool can_push_back(const ContigAllele& allele) const noexcept;
    bool can_push_back(const Allele& allele) const noexcept;
    bool can_push_front(const ContigAllele& allele) const noexcept;
    bool can_push_front(const Allele& allele) const noexcept;
    
    void push_back(const ContigAllele& allele);
    void push_front(const ContigAllele& allele);
    void push_back(ContigAllele&& allele);
    void push_front(ContigAllele&& allele);
    void push_back(const Allele& allele);
    void push_front(const Allele& allele);
    void push_back(Allele&& allele);
    void push_front(Allele&& allele);
    
    Haplotype build();
    
    friend Haplotype detail::do_copy(const Haplotype& haplotype, const GenomicRegion& region, std::true_type);
    friend Haplotype copy(const Haplotype&, const std::vector<GenomicRegion>&);
    
private:
    GenomicRegion region_;
    std::deque<ContigAllele> explicit_alleles_;
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    ContigAllele get_intervening_reference_allele(const ContigAllele& lhs, const ContigAllele& rhs) const;
    void update_region(const ContigAllele& allele) noexcept;
    void update_region(const Allele& allele);
};

// non-members

Haplotype::NucleotideSequence::size_type sequence_size(const Haplotype& haplotype) noexcept;

bool is_sequence_empty(const Haplotype& haplotype) noexcept;

bool contains(const Haplotype& lhs, const Allele& rhs);
bool contains(const Haplotype& lhs, const Haplotype& rhs);
bool includes(const Haplotype& lhs, const Allele& rhs);

template <typename MappableType>
MappableType copy(const Haplotype& haplotype, const GenomicRegion& region)
{
    return detail::do_copy(haplotype, region, std::is_same<Haplotype, std::decay_t<MappableType>> {});
}

ContigAllele copy(const Haplotype& haplotype, const ContigRegion& region);

Haplotype copy(const Haplotype& haplotype, const std::vector<GenomicRegion>& regions);

template <typename MappableType, typename Container,
          typename = std::enable_if_t<std::is_same<typename Container::value_type, Haplotype>::value>>
std::vector<MappableType> copy_unique(const Container& haplotypes, const GenomicRegion& region)
{
    std::vector<MappableType> result {};
    result.reserve(haplotypes.size());
    std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::back_inserter(result),
                   [&region] (const auto& haplotype) { return copy<MappableType>(haplotype, region); });
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

template <typename Container,
          typename = std::enable_if_t<std::is_same<typename Container::value_type, Haplotype>::value>>
std::vector<ContigAllele> copy_unique(const Container& haplotypes, const ContigRegion& region)
{
    std::vector<ContigAllele> result {};
    result.reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) {
        result.emplace_back(copy(haplotype, region));
    }
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

bool is_reference(const Haplotype& haplotype);

Haplotype expand(const Haplotype& haplotype, Haplotype::MappingDomain::Size n);
Haplotype remap(const Haplotype& haplotype, const GenomicRegion& region);

std::vector<Variant> difference(const Haplotype& lhs, const Haplotype& rhs);

bool operator==(const Haplotype& lhs, const Haplotype& rhs) noexcept;
bool operator<(const Haplotype& lhs, const Haplotype& rhs);

struct HaveSameAlleles
{
    bool operator()(const Haplotype& lhs, const Haplotype& rhs) const;
};

bool have_same_alleles(const Haplotype& lhs, const Haplotype& rhs);

bool are_equal_in_region(const Haplotype& lhs, const Haplotype& rhs, const GenomicRegion& region);

struct IsLessComplex
{
    IsLessComplex() = default;
    explicit IsLessComplex(boost::optional<Haplotype> reference);
    bool operator()(const Haplotype& lhs, const Haplotype& rhs) const;
private:
    boost::optional<Haplotype> reference_;
};

struct StrictLess
{
    bool operator()(const Haplotype& lhs, const Haplotype& rhs) const;
};

// Removes all duplicates haplotypes (w.r.t operator==) keeping the duplicate which is considered least complex w.r.t cmp.
template <typename RandomIt, typename Compare>
RandomIt remove_duplicates(RandomIt first_itr, RandomIt last_itr, const Compare& cmp)
{
    std::sort(first_itr, last_itr, StrictLess {});
    auto first_dup_itr  = first_itr;
    while (true) {
        first_dup_itr = std::adjacent_find(first_dup_itr, last_itr);
        if (first_dup_itr == last_itr) break;
        auto dup_keep_itr = (cmp(*first_dup_itr, *std::next(first_dup_itr))) ? first_dup_itr : std::next(first_dup_itr);
        auto last_dup_itr = std::next(first_dup_itr, 2);
        for (; last_dup_itr != last_itr; ++last_dup_itr) {
            if (*last_dup_itr != *first_dup_itr) {
                break;
            }
            if (cmp(*last_dup_itr, *dup_keep_itr)) {
                dup_keep_itr = last_dup_itr;
            }
        }
        std::iter_swap(first_dup_itr, dup_keep_itr);
        first_dup_itr = last_dup_itr;
    }
    return std::unique(first_itr, last_itr);
}

template <typename Container, typename Compare>
unsigned remove_duplicates(Container& haplotypes, const Compare& cmp)
{
    const auto erasable_itr = remove_duplicates(std::begin(haplotypes), std::end(haplotypes), cmp);
    const auto result = std::distance(erasable_itr, std::end(haplotypes));
    haplotypes.erase(erasable_itr, std::end(haplotypes));
    return result;
}

template <typename Container>
unsigned remove_duplicates(Container& haplotypes)
{
    return remove_duplicates(haplotypes, IsLessComplex {});
}
template <typename Container>
unsigned remove_duplicates(Container& haplotypes, Haplotype reference)
{
    IsLessComplex cmp {std::move(reference)};
    return remove_duplicates(haplotypes, cmp);
}

std::ostream& operator<<(std::ostream& os, const Haplotype& haplotype);

struct HaplotypeHash
{
    std::size_t operator()(const Haplotype& haplotype) const noexcept
    {
        return haplotype.get_hash();
    }
};

namespace debug {

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

} // namespace debug
} // namespace octopus

namespace std {

template <> struct hash<octopus::Haplotype>
{
    size_t operator()(const octopus::Haplotype& haplotype) const
    {
        return octopus::HaplotypeHash()(haplotype);
    }
};

template <> struct hash<reference_wrapper<const octopus::Haplotype>>
{
    size_t operator()(const reference_wrapper<const octopus::Haplotype> haplotype) const
    {
        return hash<octopus::Haplotype>()(haplotype);
    }
};

} // namespace std

namespace boost {

template <> struct hash<octopus::Haplotype>
{
    std::size_t operator()(const octopus::Haplotype& h) const
    {
        return std::hash<octopus::Haplotype>()(h);
    }
};

} // namespace boost

#endif
