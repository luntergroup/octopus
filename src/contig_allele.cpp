//
//  contig_allele.cpp
//  Octopus
//
//  Created by Daniel Cooke on 04/12/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "contig_allele.hpp"

#include <stdexcept>
#include <algorithm>
#include <iterator>

#include "string_utils.hpp"

ContigAllele::ContigAllele(const Allele& allele)
:
region_ {allele.get_region().get_contig_region()},
sequence_ {allele.get_sequence()}
{}

const ContigRegion& ContigAllele::get_region() const noexcept
{
    return region_;
}

const ContigAllele::SequenceType& ContigAllele::get_sequence() const noexcept
{
    return sequence_;
}

// non-member methods

ContigAllele::SizeType sequence_size(const ContigAllele& allele) noexcept
{
    return static_cast<ContigAllele::SizeType>(allele.get_sequence().size());
}

ContigAllele::SequenceType get_subsequence(const ContigAllele& allele, const ContigRegion& region)
{
    if (!::contains(allele, region)) {
        return ContigAllele::SequenceType {};
    }
    
    const auto& sequence = allele.get_sequence();
    
    if (allele.get_region() == region) {
        return sequence;
    }
    
    if (begins_equal(region, allele) && is_empty_region(region) && is_insertion(allele)) {
        const auto first = std::cbegin(sequence);
        return ContigAllele::SequenceType {first, first + sequence.size() - region_size(allele)};
    }
    
    auto first = std::cbegin(allele.get_sequence()) + get_begin(region) - get_begin(allele);
    // The minimum of the allele sequence size and region size is used as deletions will
    // result in a sequence size smaller than the region size
    return ContigAllele::SequenceType {first, first + std::min(sequence.size(), static_cast<size_t>(region_size(region)))};
}

bool contains(const ContigAllele& lhs, const ContigAllele& rhs)
{
    if (!contains(get_region(lhs), get_region(rhs))) {
        return false;
    } else if (is_empty_region(lhs)) {
        return !rhs.get_sequence().empty() && Octopus::contains(lhs.get_sequence(), rhs.get_sequence());
    } else {
        return get_subsequence(lhs, get_region(rhs)) == rhs.get_sequence();
    }
}

ContigAllele splice(const ContigAllele& allele, const ContigRegion& region)
{
    if (!contains(allele, region)) {
        throw std::logic_error {"ContigAllele: trying to splice an uncontained region"};
    }
    return ContigAllele {region, get_subsequence(allele, region)};
}

bool is_insertion(const ContigAllele& allele)
{
    return allele.get_sequence().size() > region_size(allele);
}

bool is_deletion(const ContigAllele& allele)
{
    return allele.get_sequence().size() < region_size(allele);
}

bool is_indel(const ContigAllele& allele)
{
    return is_insertion(allele) || is_deletion(allele);
}

bool operator==(const ContigAllele& lhs, const ContigAllele& rhs)
{
    return lhs.get_region() == rhs.get_region() && lhs.get_sequence() == rhs.get_sequence();
}

bool operator<(const ContigAllele& lhs, const ContigAllele& rhs)
{
    return (lhs.get_region() == rhs.get_region()) ? lhs.get_sequence() < rhs.get_sequence() :
    lhs.get_region() < rhs.get_region();
}

std::ostream& operator<<(std::ostream& os, const ContigAllele& allele)
{
    os << allele.get_region() << " " << allele.get_sequence();
    return os;
}
