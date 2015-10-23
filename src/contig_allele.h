//
//  contig_allele.h
//  Octopus
//
//  Created by Daniel Cooke on 21/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef contig_allele_h
#define contig_allele_h

#include "comparable.hpp"
#include "contig_region.hpp"
#include "genomic_region.hpp"
#include "mappable.hpp"

/*
 ContigAllele is like Allele but uses ContigRegion rather than GenomicRegion, this makes it much
 less flexible, but more efficient to use within other classes as a set with a common contig.
 */
class ContigAllele : public Comparable<Allele>, public Mappable<Allele>
{
public:
    using SizeType     = ContigRegion::SizeType;
    using SequenceType = ReferenceGenome::SequenceType;
    
    ContigAllele() = default;
    template <typename SequenceType_> ContigAllele(ContigRegion region, SequenceType_&& sequence);
    template <typename SequenceType_> ContigAllele(const GenomicRegion& region, SequenceType_&& sequence);
    ~ContigAllele() = default;
    
    ContigAllele(const ContigAllele&)            = default;
    ContigAllele& operator=(const ContigAllele&) = default;
    ContigAllele(ContigAllele&&)                 = default;
    ContigAllele& operator=(ContigAllele&&)      = default;
    
    const GenomicRegion& get_region() const noexcept;
    const SequenceType& get_sequence() const noexcept;
    
private:
    SequenceType sequence_;
    ContigRegion region_;
};

template <typename SequenceType_>
ContigAllele::ContigAllele(ContigRegion region, SequenceType_&& sequence)
: region_ {region}, sequence_ {std::forward<SequenceType_>(sequence)}
{}

template <typename SequenceType_>
ContigAllele::ContigAllele(const GenomicRegion& region, SequenceType_&& sequence)
: region_ {region.get_contig_region()}, sequence_ {std::forward<SequenceType_>(sequence)}
{}

#endif /* contig_allele_h */
