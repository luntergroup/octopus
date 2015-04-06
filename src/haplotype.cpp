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
the_reference_ {the_reference},
is_region_set_ {false},
the_reference_region_ {},
the_haplotype_ {}
{}

Haplotype::Haplotype(ReferenceGenome& the_reference, const GenomicRegion& a_region)
:
the_reference_ {the_reference},
is_region_set_ {true},
the_reference_region_ {a_region},
the_haplotype_ {}
{}

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

bool Haplotype::contains(const Variant& a_variant) const
{
    auto er = std::equal_range(std::cbegin(the_haplotype_), std::cend(the_haplotype_), a_variant);
    return (std::distance(er.first, er.second) == 0) ? false : er.first->the_sequence == a_variant.get_alternative_allele();
}

GenomicRegion Haplotype::get_region() const
{
    return (is_region_set_) ? the_reference_region_ : get_region_bounded_by_alleles();
}

Haplotype::SequenceType Haplotype::get_sequence() const
{
    return (is_region_set_) ? get_sequence(the_reference_region_) : get_sequence_bounded_by_alleles();
}

Haplotype::SequenceType Haplotype::get_sequence(const GenomicRegion& a_region) const
{
    if (the_haplotype_.empty()) {
        return the_reference_.get_sequence(a_region);
    }
    
    auto the_region_bounded_by_alleles = get_region_bounded_by_alleles();
    
    GenomicRegion reference_left_part {
        a_region.get_contig_name(),
        a_region.get_begin(),
        the_region_bounded_by_alleles.get_begin()
    };
    GenomicRegion reference_right_part {
        a_region.get_contig_name(),
        the_region_bounded_by_alleles.get_end(),
        a_region.get_end()
    };
    
    SequenceType result {};
    
    if (size(reference_left_part) > 0) {
        result += the_reference_.get_sequence(reference_left_part);
    }
    
    result += get_sequence_bounded_by_alleles();
    
    if (size(reference_right_part) > 0) {
        result += the_reference_.get_sequence(reference_right_part);
    }
    
    return result;
}

GenomicRegion Haplotype::get_region_bounded_by_alleles() const
{
    if (the_haplotype_.empty()) throw std::runtime_error {"Cannot get region from empty allele list"};
    return GenomicRegion {
        the_haplotype_.front().the_reference_region.get_contig_name(),
        the_haplotype_.front().the_reference_region.get_begin(),
        the_haplotype_.back().the_reference_region.get_end()
    };
}

Haplotype::SequenceType Haplotype::get_sequence_bounded_by_alleles() const
{
    SequenceType result {};
    
    if (the_haplotype_.empty()) return result;
    
    GenomicRegion previous_region {the_haplotype_.front().the_reference_region};
    
    for (const auto& allele : the_haplotype_) {
        auto this_region = allele.the_reference_region;
        
        if (previous_region != this_region && !overlaps(previous_region, this_region)) {
            auto intervening_reference_region   = get_intervening_region(previous_region, this_region);
            auto intervening_reference_sequence = the_reference_.get_sequence(intervening_reference_region);
            result += intervening_reference_sequence;
        }
        
        result += allele.the_sequence;
        
        previous_region = this_region;
    }
    
    return result;
}

bool operator<(const Variant& lhs, const Haplotype::Allele& rhs)
{
    return lhs.get_reference_allele_region() < rhs.the_reference_region;
}

bool operator<(const Haplotype::Allele& lhs, const Variant& rhs)
{
    return lhs.the_reference_region < rhs.get_reference_allele_region();
}
