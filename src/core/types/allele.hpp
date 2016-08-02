//
//  allele.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_allele_hpp
#define Octopus_allele_hpp

#include <string>
#include <utility>
#include <algorithm>
#include <iterator>
#include <cstddef>
#include <ostream>
#include <type_traits>
#include <functional>

#include <boost/functional/hash.hpp>

#include <basics/genomic_region.hpp>
#include <interfaces/comparable.hpp>
#include <interfaces/mappable.hpp>
#include <io/reference/reference_genome.hpp>

namespace octopus {

template <typename RegionTp, typename = EnableIfRegion<RegionTp>> class BasicAllele;

template <typename RegionTp>
class BasicAllele<RegionTp> : public Comparable<BasicAllele<RegionTp>>, public Mappable<BasicAllele<RegionTp>>
{
public:
    using RegionType         = RegionTp;
    using NucleotideSequence = ReferenceGenome::GeneticSequence;
    
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
:
sequence_ {std::forward<S>(sequence)},
region_ {std::forward<R>(region)}
{}

template <typename RegionTp>
template <typename S>
BasicAllele<RegionTp>::BasicAllele(const ContigRegion::Position begin, S&& sequence)
:
sequence_ {std::forward<S>(sequence)},
region_ {begin, static_cast<ContigRegion::Position>(begin + sequence_.size())}
{
    static_assert(std::is_same<RegionTp, ContigRegion>::value,
                  "This constructor is only for ContigAllele");
}

template <typename RegionTp>
template <typename T, typename S>
BasicAllele<RegionTp>::BasicAllele(T&& contig_name, const GenomicRegion::Position begin, S&& sequence)
:
sequence_ {std::forward<S>(sequence)},
region_ {
    std::forward<T>(contig_name), begin,
    static_cast<GenomicRegion::Position>(begin + sequence_.size())
}
{
    static_assert(std::is_same<RegionTp, GenomicRegion>::value,
                  "This constructor is only for Allele");
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

// template base non-member methods

ContigAllele demote(const Allele& allele);
ContigAllele demote(Allele&& allele);

template <typename RegionTp>
auto sequence_size(const BasicAllele<RegionTp>& allele) noexcept
{
    return allele.sequence().size();
}

template <typename RegionTp>
bool is_empty_sequence(const BasicAllele<RegionTp>& allele) noexcept
{
    return allele.sequence().empty();
}

namespace detail
{
    template <typename Sequence>
    bool is_subsequence(const Sequence& lhs, const Sequence& rhs)
    {
        using std::cbegin; using std::cend;
        return std::search(cbegin(lhs), cend(lhs), cbegin(rhs), cend(rhs)) != cend(lhs);
    }
    
    template <typename RegionTp>
    auto subsequence(const BasicAllele<RegionTp>& allele, const RegionTp& region)
    {
        using NucleotideSequence = typename BasicAllele<RegionTp>::NucleotideSequence;
        
        if (!contains(allele, region)) {
            return NucleotideSequence {};
        }
        
        const auto& sequence = allele.sequence();
        
        if (mapped_region(allele) == region) {
            return sequence;
        }
        
        if (begins_equal(region, allele) && is_empty(region) && is_insertion(allele)) {
            auto first = std::cbegin(sequence);
            return NucleotideSequence {first, std::next(first, sequence.size() - region_size(allele))};
        }
        
        auto first = std::cbegin(allele.sequence()) + begin_distance(allele, region);
        // The minimum of the allele sequence size and region size is used as deletions will
        // result in a sequence size smaller than the region size
        return NucleotideSequence {
            first, std::next(first, std::min(sequence.size(), static_cast<std::size_t>(region_size(region))))
        };
    }
}

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
    
    return detail::subsequence(lhs, rhs.mapped_region()) == rhs.sequence();
}

template <typename RegionTp>
BasicAllele<RegionTp> splice(const BasicAllele<RegionTp>& allele, const RegionTp& region)
{
    if (!contains(allele, region)) {
        throw std::logic_error {"Allele: trying to splice an uncontained region"};
    }
    return BasicAllele<RegionTp> {region, detail::subsequence(allele, region)};
}

template <typename RegionTp>
bool is_insertion(const BasicAllele<RegionTp>& allele)
{
    return allele.sequence().size() > region_size(allele);
}

template <typename RegionTp>
bool is_deletion(const BasicAllele<RegionTp>& allele)
{
    return allele.sequence().size() < region_size(allele);
}

template <typename RegionTp>
bool is_indel(const BasicAllele<RegionTp>& allele)
{
    return is_insertion(allele) || is_deletion(allele);
}

template <typename RegionTp>
auto decompose(const BasicAllele<RegionTp>& allele)
{
    std::vector<BasicAllele<RegionTp>> result {};
    
    if (is_insertion(allele)) {
        const auto& sequence = allele.sequence();
        auto insertion_size = sequence.size();
        result.reserve(insertion_size);
        
        for (unsigned i {0}; i < insertion_size; ++i) {
            result.emplace_back(mapped_region(allele), sequence.substr(i, 1));
        }
    } else if (is_deletion(allele)) {
        result.reserve(region_size(allele));
        auto decomposed_regions = decompose(mapped_region(allele));
        
        std::transform(std::begin(decomposed_regions), std::end(decomposed_regions), std::back_inserter(result),
                       [] (const auto& region) { return BasicAllele<RegionTp> {region, ""}; } );
    } else {
        result.reserve(region_size(allele));
        const auto& sequence = allele.sequence();
        unsigned i {};
        
        for (const auto& region : decompose(mapped_region(allele))) {
            result.emplace_back(region, sequence.substr(i, 1));
            ++i;
        }
    }
    
    return result;
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

template <typename T>
Allele make_allele(const std::string& region, T&& sequence, const ReferenceGenome& reference)
{
    return Allele {parse_region(region, reference), std::forward<T>(sequence)};
}

Allele make_reference_allele(const GenomicRegion& region, const ReferenceGenome& reference);
Allele make_reference_allele(const std::string& region, const ReferenceGenome& reference);

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
        size_t result {0};
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
    : std::unary_function<octopus::BasicAllele<RegionTp>, std::size_t>
    {
        std::size_t operator()(const octopus::BasicAllele<RegionTp>& a) const
        {
            return std::hash<octopus::BasicAllele<RegionTp>>()(a);
        }
    };
} // namespace boost

#endif
