//
//  haplotype.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype.h"

#include "reference_genome.h"
#include "genomic_region.h"

Haplotype::Haplotype(ReferenceGenome& the_reference)
:
the_reference_ {the_reference}
{}

GenomicRegion Haplotype::get_region() const
{
    return GenomicRegion {
        the_haplotype_.front().the_reference_region.get_contig_name(),
        the_haplotype_.front().the_reference_region.get_begin(),
        the_haplotype_.back().the_reference_region.get_end()
    };
}

Haplotype::SequenceType Haplotype::get_sequence() const
{
    SequenceType result {};
    
    GenomicRegion last_region {the_haplotype_.front().the_reference_region};
    
    for (const auto& allele : the_haplotype_) {
        auto this_region = allele.the_reference_region;
        
        if (!overlaps(last_region, this_region)) {
            auto intervening_reference_region   = get_intervening_region(last_region, this_region);
            auto intervening_reference_sequence = the_reference_.get_sequence(intervening_reference_region);
            result += intervening_reference_sequence;
        }
        
        result += allele.the_sequence;
        
        last_region = this_region;
    }
    
    return result;
}

void Haplotype::emplace_back(const Variant& a_variant)
{
    the_haplotype_.emplace_back(a_variant.get_reference_allele_region(),
                                a_variant.get_alternative_allele());
}

void Haplotype::emplace_front(const Variant& a_variant)
{
    the_haplotype_.emplace_front(a_variant.get_reference_allele_region(),
                                 a_variant.get_alternative_allele());
}
