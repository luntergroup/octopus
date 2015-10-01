//
//  allele.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "allele.hpp"

#include <algorithm> // std::min, std::transform
#include <iterator>  // std::back_inserter
#include <cstddef>   // size_t

#include "string_utils.hpp"
#include "mappable_algorithms.hpp"

const GenomicRegion& Allele::get_region() const noexcept
{
    return reference_region_;
}

const Allele::SequenceType& Allele::get_sequence() const noexcept
{
    return sequence_;
}

// non-member functions

bool is_reference(const Allele& allele, ReferenceGenome& reference)
{
    return allele.get_sequence() == reference.get_sequence(allele.get_region());
}

Allele get_reference_allele(const GenomicRegion& region, ReferenceGenome& reference)
{
    return Allele {region, reference.get_sequence(region)};
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
    
    if (begins_equal(region, allele) && empty(region) && is_insertion(allele)) {
        auto first = std::cbegin(sequence);
        return Allele::SequenceType {first, first + sequence.size() - size(allele)};
    }
    
    auto first = std::cbegin(allele.get_sequence()) + get_begin(region) - get_begin(allele);
    // The minimum of the allele sequence size and region size is used as deletions will
    // result in a sequence size smaller than the region size
    return Allele::SequenceType {first, first + std::min(sequence.size(), static_cast<size_t>(size(region)))};
}

Allele splice(const Allele& allele, const GenomicRegion& region)
{
    if (!contains(allele, region)) {
        throw std::runtime_error {"cannot splice region " + to_string(region) +
            " from Allele as the region is not contained"};
    }
    return Allele {region, get_subsequence(allele, region)};
}

bool contains(const Allele& lhs, const Allele& rhs)
{
    if (!contains(get_region(lhs), get_region(rhs))) {
        return false;
    } else if (empty(lhs)) {
        // If the alleles are both insertions then both regions will be the same so we can only test
        // if the inserted rhs sequence is a subsequence of the lhs sequence. The rhs sequence
        // is required to be non-empty otherwise it would be a subsequence of everything.
        return !rhs.get_sequence().empty() && ::contains(lhs.get_sequence(), rhs.get_sequence());
    } else {
        return get_subsequence(lhs, rhs.get_region()) == rhs.get_sequence();
    }
}

bool is_insertion(const Allele& allele)
{
    return allele.get_sequence().size() > size(allele);
}

bool is_deletion(const Allele& allele)
{
    return allele.get_sequence().size() < size(allele);
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
        result.reserve(size(allele));
        auto decomposed_regions = decompose(get_region(allele));
        
        std::transform(std::begin(decomposed_regions), std::end(decomposed_regions), std::back_inserter(result),
                       [] (const auto& region) { return Allele {region, ""}; } );
    } else {
        result.reserve(size(allele));
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
