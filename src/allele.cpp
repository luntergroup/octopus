//
//  allele.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "allele.h"

#include <algorithm> // std::min, std::transform
#include <iterator>  // std::back_inserter

#include "string_utils.h"
#include "mappable_algorithms.h"

const GenomicRegion& Allele::get_region() const noexcept
{
    return the_reference_region_;
}

const Allele::SequenceType& Allele::get_sequence() const noexcept
{
    return the_sequence_;
}

bool is_reference(const Allele& an_allele, ReferenceGenome& the_reference)
{
    return an_allele.get_sequence() == the_reference.get_sequence(an_allele.get_region());
}

Allele get_reference_allele(const GenomicRegion& a_region, ReferenceGenome& the_reference)
{
    return Allele {a_region, the_reference.get_sequence(a_region)};
}

Allele::SequenceType get_subsequence(const Allele& an_allele, const GenomicRegion& a_region)
{
    if (!contains(an_allele, a_region)) {
        return Allele::SequenceType {};
    } if (an_allele.get_region() == a_region) {
        return an_allele.get_sequence();
    } else {
        auto first = std::cbegin(an_allele.get_sequence()) + (get_begin(a_region) - get_begin(an_allele));
        
        // The minimum of the allele sequence size and region size is used as deletions will
        // result in a sequence size smaller than the region size
        return Allele::SequenceType {first, first +
            std::min(an_allele.get_sequence().size(), static_cast<std::size_t>(size(a_region)))};
    }
}

Allele splice(const Allele& an_allele, const GenomicRegion& a_region)
{
    if (!contains(an_allele, a_region)) {
        throw std::runtime_error {"cannot splice Allele region that is not contained"};
    }
    return Allele {a_region, get_subsequence(an_allele, a_region)};
}

bool contains(const Allele& lhs, const Allele& rhs)
{
    if (!contains(lhs.get_region(), rhs.get_region())) {
        return false;
    } else if (empty(lhs.get_region())) {
        // If the alleles are both insertions then both regions will be the same so we can only test
        // if the inserted rhs sequence is a subsequence of the lhs sequence. The rhs sequence
        // is required to be non-empty otherwise it would be a subsequence of everything.
        return !rhs.get_sequence().empty() && ::contains(lhs.get_sequence(), rhs.get_sequence());
    } else {
        return get_subsequence(lhs, rhs.get_region()) == rhs.get_sequence();
    }
}

bool is_insertion(const Allele& an_allele)
{
    return empty(an_allele.get_region()) && an_allele.get_sequence().size() > 0;
}

bool is_deletion(const Allele& an_allele)
{
    return !empty(an_allele.get_region()) && an_allele.get_sequence().size() == 0;
}

std::vector<Allele> decompose(const Allele& an_allele)
{
    std::vector<Allele> result {};
    
    if (is_insertion(an_allele)) {
        const auto& sequence = an_allele.get_sequence();
        auto insertion_size = sequence.size();
        result.reserve(insertion_size);
        
        for (unsigned i {0}; i < insertion_size; ++i) {
            result.emplace_back(an_allele.get_region(), sequence.substr(i, 1));
        }
    } else if (is_deletion(an_allele)) {
        result.reserve(size(an_allele));
        auto decomposed_regions = decompose(an_allele.get_region());
        
        std::transform(decomposed_regions.begin(), decomposed_regions.end(), std::back_inserter(result),
                       [] (const auto& region) { return Allele {region, ""}; } );
    } else {
        result.reserve(size(an_allele));
        const auto& sequence = an_allele.get_sequence();
        unsigned i {};
        
        for (auto& region : decompose(an_allele.get_region())) {
            result.emplace_back(region, sequence.substr(i, 1));
            ++i;
        }
    }
    
    return result;
}
