// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef allele_hpp
#define allele_hpp

#include <string>
#include <utility>
#include <algorithm>
#include <iterator>
#include <cstddef>
#include <ostream>
#include <type_traits>
#include <functional>
#include <cassert>

#include <boost/functional/hash.hpp>

#include "concepts/comparable.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "io/reference/reference_genome.hpp"

namespace octopus {

/**
 An Allele is simply a mapped chuck of sequence. Although the mapping is only fully defined
 by GenomicRegion, when the user knows the contig mapping it may be more efficient to store
 ContigAllele which uses the less informative ContigRegion mapping.
 */
template <typename RegionTp>
class BasicAllele : public Comparable<BasicAllele<RegionTp>>, public Mappable<BasicAllele<RegionTp>>
{
public:
    using MappingDomain      = RegionTp;
    using NucleotideSequence = ReferenceGenome::GeneticSequence;
    
    //static_assert(is_region<RegionType>, "Not a region");
    
    BasicAllele() = default;
    
    template <typename R, typename S>
    BasicAllele(R&& region, S&& sequence);
    template <typename S> // RegionTp=ContigRegion
    BasicAllele(ContigRegion::Position begin, S&& sequence);
    template <typename T, typename S> // RegionTp=GenomicRegion
    BasicAllele(T&& contig_name, GenomicRegion::Position begin, S&& sequence);
    
    BasicAllele(const BasicAllele&)            = default;
    BasicAllele& operator=(const BasicAllele&) = default;
    BasicAllele(BasicAllele&&)                 = default;
    BasicAllele& operator=(BasicAllele&&)      = default;
    
    ~BasicAllele() = default;
    
    const RegionTp& mapped_region() const noexcept;
    
    const NucleotideSequence& sequence() const noexcept;
    NucleotideSequence& sequence() noexcept;
    
    friend BasicAllele<ContigRegion> demote(BasicAllele<GenomicRegion>&&);
    
private:
    NucleotideSequence sequence_;
    RegionTp region_;
};

// concrete types

using ContigAllele = BasicAllele<ContigRegion>;
using Allele       = BasicAllele<GenomicRegion>;

// template base member methods

template <typename RegionTp>
template <typename R, typename S>
BasicAllele<RegionTp>::BasicAllele(R&& region, S&& sequence)
: sequence_ {std::forward<S>(sequence)}
, region_ {std::forward<R>(region)}
{}

template <typename RegionTp>
template <typename S>
BasicAllele<RegionTp>::BasicAllele(const ContigRegion::Position begin, S&& sequence)
: sequence_ {std::forward<S>(sequence)}
, region_ {begin, static_cast<ContigRegion::Position>(begin + sequence_.size())}
{
    static_assert(std::is_same<RegionTp, ContigRegion>::value,
                  "This constructor is only for ContigAllele");
}

template <typename RegionTp>
template <typename T, typename S>
BasicAllele<RegionTp>::BasicAllele(T&& contig_name, const GenomicRegion::Position begin, S&& sequence)
: sequence_ {std::forward<S>(sequence)}
, region_ {
    std::forward<T>(contig_name), begin,
    static_cast<GenomicRegion::Position>(begin + sequence_.size())
}
{
    static_assert(std::is_same<RegionTp, GenomicRegion>::value, "This constructor is only for Allele");
}

template <typename RegionTp>
const RegionTp& BasicAllele<RegionTp>::mapped_region() const noexcept
{
    return region_;
}

template <typename RegionTp>
const typename BasicAllele<RegionTp>::NucleotideSequence& BasicAllele<RegionTp>::sequence() const noexcept
{
    return sequence_;
}

template <typename RegionTp>
typename BasicAllele<RegionTp>::NucleotideSequence& BasicAllele<RegionTp>::sequence() noexcept
{
    return sequence_;
}

// template base non-member methods

ContigAllele demote(const Allele& allele);
ContigAllele demote(Allele&& allele);

template <typename RegionTp>
auto sequence_size(const BasicAllele<RegionTp>& allele) noexcept
{
    return allele.sequence().size();
}

template <typename RegionTp>
bool is_sequence_empty(const BasicAllele<RegionTp>& allele) noexcept
{
    return allele.sequence().empty();
}

namespace detail {

template <typename Sequence>
bool is_subsequence(const Sequence& lhs, const Sequence& rhs)
{
    using std::cbegin; using std::cend;
    return std::search(cbegin(lhs), cend(lhs), cbegin(rhs), cend(rhs)) != cend(lhs);
}

template <typename RegionTp>
auto copy_sequence(const BasicAllele<RegionTp>& allele, const RegionTp& region)
{
    using NucleotideSequence = typename BasicAllele<RegionTp>::NucleotideSequence;
    assert(contains(allele, region));
    const auto& sequence = allele.sequence();
    if (mapped_region(allele) == region) return sequence;
    const auto region_offset = static_cast<std::size_t>(begin_distance(allele, region));
    auto first_base_itr = std::cbegin(sequence), last_base_itr = std::cend(sequence);
    if (is_deletion(allele)) {
        if (!is_sequence_empty(allele)) {
            first_base_itr = std::cbegin(allele.sequence());
            const auto num_deleted_based = region_size(allele) - sequence_size(allele);
            if (size(region) <= num_deleted_based) {
                last_base_itr = first_base_itr;
            } else {
                last_base_itr = std::next(first_base_itr, size(region) - num_deleted_based);
            }
        }
    } else {
        first_base_itr = std::next(std::cbegin(allele.sequence()), region_offset);
        if (is_simple_insertion(allele)) {
            const auto num_trailing_bases = static_cast<std::size_t>(end_distance(region, allele));
            const auto num_subsequence_bases = sequence_size(allele) - region_offset - num_trailing_bases;
            last_base_itr = std::next(first_base_itr, num_subsequence_bases);
        } else {
            last_base_itr = std::next(first_base_itr, size(region));
        }
    }
    assert(first_base_itr <= last_base_itr);
    assert(first_base_itr >= std::cbegin(sequence));
    assert(last_base_itr <= std::cend(sequence));
    return NucleotideSequence {first_base_itr, last_base_itr};
}

template <typename RegionTp>
auto copy_sequence(BasicAllele<RegionTp>&& allele, const RegionTp& region)
{
    assert(contains(allele, region));
    auto sequence = std::move(allele.sequence());
    if (mapped_region(allele) == region) return sequence;
    const auto region_offset = static_cast<std::size_t>(begin_distance(allele, region));
    auto first_base_itr = std::cbegin(sequence), last_base_itr = std::cend(sequence);
    const auto region_size = static_cast<std::size_t>(size(region));
    if (is_deletion(allele)) {
        if (!is_sequence_empty(allele)) {
            const auto base_offset = std::min(region_offset, sequence_size(allele));
            first_base_itr = std::next(std::cbegin(allele.sequence()), base_offset);
            const auto num_remaining_bases = std::min(region_size, sequence_size(allele) - base_offset);
            last_base_itr = std::next(first_base_itr, num_remaining_bases);
        }
    } else {
        first_base_itr = std::next(std::cbegin(allele.sequence()), region_offset);
        if (is_insertion(allele)) {
            const auto num_trailing_bases = static_cast<std::size_t>(end_distance(region, allele));
            const auto num_subsequence_bases = sequence_size(allele) - region_offset - num_trailing_bases;
            last_base_itr = std::next(first_base_itr, num_subsequence_bases);
        } else {
            last_base_itr = std::next(first_base_itr, region_size);
        }
    }
    assert(first_base_itr <= last_base_itr);
    assert(first_base_itr >= std::cbegin(sequence));
    assert(last_base_itr <= std::cend(sequence));
    if (first_base_itr == last_base_itr) {
        sequence.clear();
    } else {
        if (last_base_itr != std::cend(sequence)) {
            sequence.erase(std::next(last_base_itr), std::cend(sequence));
        }
        if (first_base_itr != std::cbegin(sequence)) {
            sequence.erase(std::cbegin(sequence), first_base_itr);
        }
    }
    return sequence;
}

} // namespace detail

// Note this hides Mappable::contains
template <typename RegionTp>
bool contains(const BasicAllele<RegionTp>& lhs, const BasicAllele<RegionTp>& rhs)
{
    if (!contains(mapped_region(lhs), mapped_region(rhs))) {
        return false;
    }
    if (is_empty_region(lhs)) {
        // If the alleles are both insertions then both regions will be the same so we can only test
        // if the inserted rhs sequence is a subsequence of the lhs sequence. The rhs sequence
        // is required to be non-empty otherwise it would be a subsequence of everything.
        return !rhs.sequence().empty() && detail::is_subsequence(lhs.sequence(), rhs.sequence());
    }
    return detail::copy_sequence(lhs, rhs.mapped_region()) == rhs.sequence();
}

template <typename RegionTp>
BasicAllele<RegionTp> copy(const BasicAllele<RegionTp>& allele, const RegionTp& region)
{
    if (!contains(allele, region)) {
        throw std::logic_error {"Allele: trying to copy an uncontained region"};
    }
    return BasicAllele<RegionTp> {region, detail::copy_sequence(allele, region)};
}

template <typename RegionTp>
BasicAllele<RegionTp> copy(BasicAllele<RegionTp>&& allele, const RegionTp& region)
{
    if (!contains(allele, region)) {
        throw std::logic_error {"Allele: trying to copy an uncontained region"};
    }
    return BasicAllele<RegionTp> {region, detail::copy_sequence(std::move(allele), region)};
}

template <typename RegionTp>
bool is_insertion(const BasicAllele<RegionTp>& allele) noexcept
{
    return allele.sequence().size() > region_size(allele);
}

template <typename RegionTp>
bool is_deletion(const BasicAllele<RegionTp>& allele) noexcept
{
    return allele.sequence().size() < region_size(allele);
}

template <typename RegionTp>
bool is_indel(const BasicAllele<RegionTp>& allele) noexcept
{
    return is_insertion(allele) || is_deletion(allele);
}

template <typename RegionTp>
bool is_simple_insertion(const BasicAllele<RegionTp>& allele) noexcept
{
    return is_insertion(allele) && is_empty_region(allele);
}

template <typename RegionTp>
bool is_simple_deletion(const BasicAllele<RegionTp>& allele) noexcept
{
    return is_deletion(allele) && is_sequence_empty(allele);
}

template <typename RegionTp>
bool is_simple_indel(const BasicAllele<RegionTp>& allele) noexcept
{
    return is_simple_insertion(allele) || is_simple_deletion(allele);
}

template <typename RegionTp>
typename BasicAllele<RegionTp>::MappingDomain::Size reference_distance(const BasicAllele<RegionTp>& allele)
{
    if (is_insertion(allele)) {
        return sequence_size(allele) - region_size(allele);
    } else if (is_deletion(allele)) {
        return region_size(allele) - sequence_size(allele);
    } else {
        return 0;
    }
}

template <typename RegionTp>
bool operator==(const BasicAllele<RegionTp>& lhs, const BasicAllele<RegionTp>& rhs)
{
    return lhs.mapped_region() == rhs.mapped_region() && lhs.sequence() == rhs.sequence();
}

template <typename RegionTp>
bool operator<(const BasicAllele<RegionTp>& lhs, const BasicAllele<RegionTp>& rhs)
{
    if (lhs.mapped_region() == rhs.mapped_region()) {
        return lhs.sequence() < rhs.sequence();
    }
    return lhs.mapped_region() < rhs.mapped_region();
}

template <typename RegionTp>
std::ostream& operator<<(std::ostream& os, const BasicAllele<RegionTp>& allele)
{
    os << allele.mapped_region() << " " << allele.sequence();
    return os;
}

// Allele specific methods

bool is_reference(const Allele& allele, const ReferenceGenome& reference);

Allele make_reference_allele(const GenomicRegion& region, const ReferenceGenome& reference);

std::vector<Allele> make_reference_alleles(const std::vector<GenomicRegion>& regions,
                                           const ReferenceGenome& reference);

std::vector<Allele> make_positional_reference_alleles(const GenomicRegion& region,
                                                      const ReferenceGenome& reference);

struct AlleleHash
{
    template <typename RegionTp>
    std::size_t operator()(const BasicAllele<RegionTp>& allele) const
    {
        using boost::hash_combine;
        std::size_t result {};
        hash_combine(result, std::hash<RegionTp>()(allele.mapped_region()));
        using T = typename BasicAllele<RegionTp>::NucleotideSequence;
        hash_combine(result, std::hash<T>()(allele.sequence()));
        return result;
    }
};

} // namespace octopus

namespace std {

template <typename RegionTp> struct hash<octopus::BasicAllele<RegionTp>>
{
    size_t operator()(const octopus::BasicAllele<RegionTp>& allele) const
    {
        return octopus::AlleleHash()(allele);
    }
};

} // namespace std

namespace boost {

template <typename RegionTp> struct hash<octopus::BasicAllele<RegionTp>>
{
    std::size_t operator()(const octopus::BasicAllele<RegionTp>& a) const
    {
        return std::hash<octopus::BasicAllele<RegionTp>>()(a);
    }
};

} // namespace boost

#endif
