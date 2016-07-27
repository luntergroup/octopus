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

#include <boost/functional/hash.hpp>

#include "genomic_region.hpp"
#include "comparable.hpp"
#include "mappable.hpp"
#include "reference_genome.hpp"

template <typename RegionTp, typename = EnableIfRegion<RegionTp>> class BaseAllele;

template <typename RegionTp>
class BaseAllele<RegionTp> : public Comparable<BaseAllele<RegionTp>>, public Mappable<BaseAllele<RegionTp>>
{
public:
    using SizeType     = typename RegionTp::SizeType;
    using SequenceType = ReferenceGenome::SequenceType;
    
    BaseAllele() = default;
    
    template <typename R, typename S>
    BaseAllele(R&& region, S&& sequence);
    template <typename S> // RegionTp=ContigRegion
    BaseAllele(SizeType begin_pos, S&& sequence);
    template <typename T, typename S> // RegionTp=GenomicRegion
    BaseAllele(T&& contig_name, SizeType begin_pos, S&& sequence);
    
    BaseAllele(const BaseAllele&)            = default;
    BaseAllele& operator=(const BaseAllele&) = default;
    BaseAllele(BaseAllele&&)                 = default;
    BaseAllele& operator=(BaseAllele&&)      = default;
    
    ~BaseAllele() = default;
    
    const RegionTp& mapped_region() const noexcept;
    
    const SequenceType& sequence() const noexcept;
    
    friend BaseAllele<ContigRegion> demote(BaseAllele<GenomicRegion>&&);
    
private:
    SequenceType sequence_;
    RegionTp region_;
};

// concrete types

using ContigAllele = BaseAllele<ContigRegion>;
using Allele       = BaseAllele<GenomicRegion>;

// template base member methods

template <typename RegionTp>
template <typename R, typename S>
BaseAllele<RegionTp>::BaseAllele(R&& region, S&& sequence)
:
sequence_ {std::forward<S>(sequence)},
region_ {std::forward<R>(region)}
{}

template <typename RegionTp>
template <typename S>
BaseAllele<RegionTp>::BaseAllele(const SizeType begin_pos, S&& sequence)
:
sequence_ {std::forward<S>(sequence)},
region_ {begin_pos, static_cast<SizeType>(begin_pos + sequence_.size())}
{}

template <typename RegionTp>
template <typename T, typename S>
BaseAllele<RegionTp>::BaseAllele(T&& contig_name, const SizeType begin_pos, S&& sequence)
:
sequence_ {std::forward<S>(sequence)},
region_ {std::forward<T>(contig_name), begin_pos, static_cast<SizeType>(begin_pos + sequence_.size())}
{}

template <typename RegionTp>
const RegionTp& BaseAllele<RegionTp>::mapped_region() const noexcept
{
    return region_;
}

template <typename RegionTp>
const typename BaseAllele<RegionTp>::SequenceType& BaseAllele<RegionTp>::sequence() const noexcept
{
    return sequence_;
}

// template base non-member methods

ContigAllele demote(const Allele& allele);
ContigAllele demote(Allele&& allele);

template <typename RegionTp>
auto sequence_size(const BaseAllele<RegionTp>& allele) noexcept
{
    return static_cast<typename BaseAllele<RegionTp>::SizeType>(allele.sequence().size());
}

template <typename RegionTp>
auto is_empty_sequence(const BaseAllele<RegionTp>& allele) noexcept
{
    return allele.sequence().empty();
}

namespace detail
{
    template <typename SequenceType>
    bool is_subsequence(const SequenceType& lhs, const SequenceType& rhs)
    {
        using std::cbegin; using std::cend;
        return std::search(cbegin(lhs), cend(lhs), cbegin(rhs), cend(rhs)) != cend(lhs);
    }
    
    template <typename RegionTp>
    auto subsequence(const BaseAllele<RegionTp>& allele, const RegionTp& region)
    {
        using ResultType = typename BaseAllele<RegionTp>::SequenceType;
        
        if (!contains(allele, region)) {
            return ResultType {};
        }
        
        const auto& sequence = allele.sequence();
        
        if (mapped_region(allele) == region) {
            return sequence;
        }
        
        if (begins_equal(region, allele) && is_empty(region) && is_insertion(allele)) {
            auto first = std::cbegin(sequence);
            return ResultType {first, first + sequence.size() - region_size(allele)};
        }
        
        auto first = std::cbegin(allele.sequence()) + begin_distance(allele, region);
        // The minimum of the allele sequence size and region size is used as deletions will
        // result in a sequence size smaller than the region size
        return ResultType {first, first + std::min(sequence.size(), static_cast<size_t>(region_size(region)))};
    }
}

template <typename RegionTp>
bool contains(const BaseAllele<RegionTp>& lhs, const BaseAllele<RegionTp>& rhs)
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
BaseAllele<RegionTp> splice(const BaseAllele<RegionTp>& allele, const RegionTp& region)
{
    if (!contains(allele, region)) {
        throw std::logic_error {"Allele: trying to splice an uncontained region"};
    }
    return BaseAllele<RegionTp> {region, detail::subsequence(allele, region)};
}

template <typename RegionTp>
bool is_insertion(const BaseAllele<RegionTp>& allele)
{
    return allele.sequence().size() > region_size(allele);
}

template <typename RegionTp>
bool is_deletion(const BaseAllele<RegionTp>& allele)
{
    return allele.sequence().size() < region_size(allele);
}

template <typename RegionTp>
bool is_indel(const BaseAllele<RegionTp>& allele)
{
    return is_insertion(allele) || is_deletion(allele);
}

template <typename RegionTp>
auto decompose(const BaseAllele<RegionTp>& allele)
{
    std::vector<BaseAllele<RegionTp>> result {};
    
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
                       [] (const auto& region) { return BaseAllele<RegionTp> {region, ""}; } );
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
bool operator==(const BaseAllele<RegionTp>& lhs, const BaseAllele<RegionTp>& rhs)
{
    return lhs.mapped_region() == rhs.mapped_region() && lhs.sequence() == rhs.sequence();
}

template <typename RegionTp>
bool operator<(const BaseAllele<RegionTp>& lhs, const BaseAllele<RegionTp>& rhs)
{
    if (lhs.mapped_region() == rhs.mapped_region()) {
        return lhs.sequence() < rhs.sequence();
    }
    return lhs.mapped_region() < rhs.mapped_region();
}

namespace std {
    template <typename RegionTp> struct hash<BaseAllele<RegionTp>>
    {
        size_t operator()(const BaseAllele<RegionTp>& allele) const
        {
            using boost::hash_combine;
            size_t result {0};
            hash_combine(result, hash<RegionTp>()(allele.mapped_region()));
            hash_combine(result, hash<typename BaseAllele<RegionTp>::SequenceType>()(allele.sequence()));
            return result;
        }
    };
} // namespace std

namespace boost
{
    template <typename RegionTp> struct hash<BaseAllele<RegionTp>> : std::unary_function<BaseAllele<RegionTp>, std::size_t>
    {
        std::size_t operator()(const BaseAllele<RegionTp>& a) const
        {
            return std::hash<BaseAllele<RegionTp>>()(a);
        }
    };
} // namespace boost

template <typename RegionTp>
std::ostream& operator<<(std::ostream& os, const BaseAllele<RegionTp>& allele)
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

#endif
