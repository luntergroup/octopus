//
//  allele.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "allele.hpp"

#include <algorithm>
#include <iterator>
#include <cstddef>

#include "string_utils.hpp"
#include "mappable_algorithms.hpp"

// public methods

const GenomicRegion& Allele::get_region() const noexcept
{
    return reference_region_;
}

const Allele::SequenceType& Allele::get_sequence() const noexcept
{
    return sequence_;
}

// non-member functions

Allele::SizeType sequence_size(const Allele& allele) noexcept
{
    return static_cast<Allele::SizeType>(allele.get_sequence().size());
}

bool is_reference(const Allele& allele, const ReferenceGenome& reference)
{
    return allele.get_sequence() == reference.get_sequence(allele.get_region());
}

Allele get_reference_allele(const GenomicRegion& region, const ReferenceGenome& reference)
{
    return Allele {region, reference.get_sequence(region)};
}

std::vector<Allele> get_reference_alleles(const std::vector<GenomicRegion>& regions,
                                          const ReferenceGenome& reference)
{
    std::vector<Allele> result {};
    result.reserve(regions.size());
    
    std::transform(std::cbegin(regions), std::cend(regions), std::back_inserter(result),
                   [&reference] (const auto& region) { return get_reference_allele(region, reference); });
    
    return result;
}

std::vector<Allele> get_positional_reference_alleles(const GenomicRegion& region,
                                                     const ReferenceGenome& reference)
{
    auto sequence  = reference.get_sequence(region);
    auto positions = decompose(region);
    
    std::vector<Allele> result {};
    result.reserve(positions.size());
    
    std::transform(std::cbegin(positions), std::cend(positions), std::cbegin(sequence), std::back_inserter(result),
                   [&reference] (const auto& region, auto base) {
                       return Allele {region, base};
                   });
    
    return result;
}

Allele::SequenceType get_subsequence(const Allele& allele, const GenomicRegion& region)
{
    if (!contains(allele, region)) {
        return Allele::SequenceType {};
    }
    
    const auto& sequence = allele.get_sequence();
    
    if (get_region(allele) == region) {
        return sequence;
    }
    
    if (begins_equal(region, allele) && is_empty(region) && is_insertion(allele)) {
        auto first = std::cbegin(sequence);
        return Allele::SequenceType {first, first + sequence.size() - region_size(allele)};
    }
    
    auto first = std::cbegin(allele.get_sequence()) + get_begin(region) - get_begin(allele);
    // The minimum of the allele sequence size and region size is used as deletions will
    // result in a sequence size smaller than the region size
    return Allele::SequenceType {first, first + std::min(sequence.size(), static_cast<size_t>(region_size(region)))};
}

bool contains(const Allele& lhs, const Allele& rhs)
{
    if (!contains(get_region(lhs), get_region(rhs))) {
        return false;
    } else if (is_empty(lhs)) {
        // If the alleles are both insertions then both regions will be the same so we can only test
        // if the inserted rhs sequence is a subsequence of the lhs sequence. The rhs sequence
        // is required to be non-empty otherwise it would be a subsequence of everything.
        return !rhs.get_sequence().empty() && Octopus::contains(lhs.get_sequence(), rhs.get_sequence());
    } else {
        return get_subsequence(lhs, rhs.get_region()) == rhs.get_sequence();
    }
}

Allele splice(const Allele& allele, const GenomicRegion& region)
{
    if (!contains(allele, region)) {
        throw std::logic_error {"Allele: trying to splice an uncontained region"};
    }
    return Allele {region, get_subsequence(allele, region)};
}

bool is_insertion(const Allele& allele)
{
    return allele.get_sequence().size() > region_size(allele);
}

bool is_deletion(const Allele& allele)
{
    return allele.get_sequence().size() < region_size(allele);
}

bool is_indel(const Allele& allele)
{
    return is_insertion(allele) || is_deletion(allele);
}

std::vector<Allele> decompose(const Allele& allele)
{
    std::vector<Allele> result {};
    
    if (is_insertion(allele)) {
        const auto& sequence = allele.get_sequence();
        auto insertion_size = sequence.size();
        result.reserve(insertion_size);
        
        for (unsigned i {0}; i < insertion_size; ++i) {
            result.emplace_back(get_region(allele), sequence.substr(i, 1));
        }
    } else if (is_deletion(allele)) {
        result.reserve(region_size(allele));
        auto decomposed_regions = decompose(get_region(allele));
        
        std::transform(std::begin(decomposed_regions), std::end(decomposed_regions), std::back_inserter(result),
                       [] (const auto& region) { return Allele {region, ""}; } );
    } else {
        result.reserve(region_size(allele));
        const auto& sequence = allele.get_sequence();
        unsigned i {};
        
        for (const auto& region : decompose(get_region(allele))) {
            result.emplace_back(region, sequence.substr(i, 1));
            ++i;
        }
    }
    
    return result;
}

bool operator==(const Allele& lhs, const Allele& rhs)
{
    return lhs.get_region() == rhs.get_region() && lhs.get_sequence() == rhs.get_sequence();
}

bool operator<(const Allele& lhs, const Allele& rhs)
{
    return (lhs.get_region() == rhs.get_region()) ? lhs.get_sequence() < rhs.get_sequence() :
    lhs.get_region() < rhs.get_region();
}

std::ostream& operator<<(std::ostream& os, const Allele& allele)
{
    os << allele.get_region() << " " << allele.get_sequence();
    return os;
}
