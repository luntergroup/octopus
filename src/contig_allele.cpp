//
//  contig_allele.cpp
//  Octopus
//
//  Created by Daniel Cooke on 04/12/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "contig_allele.hpp"

ContigAllele::ContigAllele(const Allele& allele)
:
region_ {allele.get_region().get_contig_region()},
sequence_ {allele.get_sequence()}
{}

//ContigAllele::ContigAllele(Allele&& allele)
//:
//region_ {allele.get_region().get_contig_region()},
//sequence_ {allele.get_sequence()}
//{}

const ContigRegion& ContigAllele::get_region() const noexcept
{
    return region_;
}

const ContigAllele::SequenceType& ContigAllele::get_sequence() const noexcept
{
    return sequence_;
}

// non-member methods

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
