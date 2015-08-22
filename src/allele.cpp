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
    return reference_region_;
}

const Allele::SequenceType& Allele::get_sequence() const noexcept
{
    return sequence_;
}

// non-member functions

bool is_reference(const Allele& an_allele, ReferenceGenome& the_reference)
{
    return an_allele.get_sequence() == the_reference.get_sequence(an_allele.get_region());
}

Allele get_reference_allele(const GenomicRegion& region, ReferenceGenome& the_reference)
{
    return Allele {region, the_reference.get_sequence(region)};
}

Allele::SequenceType get_subsequence(const Allele& an_allele, const GenomicRegion& region)
{
    if (!contains(an_allele, region)) {
        return Allele::SequenceType {};
    } if (get_region(an_allele) == region) {
        return an_allele.get_sequence();
    } else {
        auto first = std::cbegin(an_allele.get_sequence()) + (get_begin(region) - get_begin(an_allele));
        
        // The minimum of the allele sequence size and region size is used as deletions will
        // result in a sequence size smaller than the region size
        return Allele::SequenceType {first, first +
            std::min(an_allele.get_sequence().size(), static_cast<std::size_t>(size(region)))};
    }
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

Allele splice(const Allele& an_allele, const GenomicRegion& region)
{
    if (!contains(an_allele, region)) {
        throw std::runtime_error {"cannot splice Allele region that is not contained"};
    }
    return Allele {region, get_subsequence(an_allele, region)};
}

bool is_insertion(const Allele& an_allele)
{
    return empty(an_allele) && an_allele.get_sequence().size() > 0;
}

bool is_deletion(const Allele& an_allele)
{
    return !empty(an_allele) && an_allele.get_sequence().size() == 0;
}

std::vector<Allele> decompose(const Allele& an_allele)
{
    std::vector<Allele> result {};
    
    if (is_insertion(an_allele)) {
        const auto& sequence = an_allele.get_sequence();
        auto insertion_size = sequence.size();
        result.reserve(insertion_size);
        
        for (unsigned i {0}; i < insertion_size; ++i) {
            result.emplace_back(get_region(an_allele), sequence.substr(i, 1));
        }
    } else if (is_deletion(an_allele)) {
        result.reserve(size(an_allele));
        auto decomposed_regions = decompose(get_region(an_allele));
        
        std::transform(std::begin(decomposed_regions), std::end(decomposed_regions), std::back_inserter(result),
                       [] (const auto& region) { return Allele {region, ""}; } );
    } else {
        result.reserve(size(an_allele));
        const auto& sequence = an_allele.get_sequence();
        unsigned i {};
        
        for (const auto& region : decompose(get_region(an_allele))) {
            result.emplace_back(region, sequence.substr(i, 1));
            ++i;
        }
    }
    
    return result;
}

std::ostream& operator<<(std::ostream& os, const Allele& an_allele)
{
    os << an_allele.get_region() << " " << an_allele.get_sequence();
    return os;
}
