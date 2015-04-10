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

/**
 I currently only store sequence data explictly added through the emplace methods. Any intervening
 or flanking regions are assumed to be reference and are retreived as needed. This saves memory,
 but it might cause too many file reads if a lot of complex containment queries are made.
 */

Haplotype::Haplotype(ReferenceGenome& the_reference)
:
the_reference_ {the_reference},
is_region_set_ {false},
the_reference_region_ {},
the_explicit_alleles_ {}
{}

Haplotype::Haplotype(ReferenceGenome& the_reference, const GenomicRegion& a_region)
:
the_reference_ {the_reference},
is_region_set_ {true},
the_reference_region_ {a_region},
the_explicit_alleles_ {}
{}

bool Haplotype::contains(const GenomicRegion& the_allele_region, const SequenceType& the_allele_sequence) const
{
    if (the_explicit_alleles_.empty()) {
        return (!is_region_set_ || !::contains(get_region(), the_allele_region)) ? false :
                            the_allele_sequence == the_reference_.get_sequence(the_allele_region);
    } else if (get_region() == the_allele_region) {
        return the_allele_sequence == get_sequence();
    }
    
    auto er = std::equal_range(std::cbegin(the_explicit_alleles_), std::cend(the_explicit_alleles_), the_allele_region);
    
    if (std::distance(er.first, er.second) != 0) {
        return er.first->get_sequence() == the_allele_sequence;
    } else {
        auto the_region_bounded_by_alleles = get_region_bounded_by_alleles();
        
        auto reference_region_left_part = get_left_overhang(the_allele_region, the_region_bounded_by_alleles);
        
        if (::contains(reference_region_left_part, the_allele_region)) {
            return the_reference_.get_sequence(the_allele_region) == the_allele_sequence;
        }
        
        auto reference_region_right_part = get_right_overhang(the_allele_region, the_region_bounded_by_alleles);
        
        if (::contains(reference_region_right_part, the_allele_region)) {
            return the_reference_.get_sequence(the_allele_region) == the_allele_sequence;
        }
        
        SequenceType haplotype_sub_sequence {};
        AlleleIterator first_allele_it;
        
        if (begins_before(the_allele_region, the_region_bounded_by_alleles)) {
            haplotype_sub_sequence += the_reference_.get_sequence(reference_region_left_part);
            first_allele_it = std::cbegin(the_explicit_alleles_);
        } else {
            first_allele_it = std::lower_bound(std::cbegin(the_explicit_alleles_), std::cend(the_explicit_alleles_), the_allele_region);
        }
        
        if (begins_before(the_allele_region, first_allele_it->get_region())) {
            if (overlaps(the_allele_region, first_allele_it->get_region())) {
                auto reference_region_before_next_allele = get_left_overhang(the_allele_region,
                                                                             first_allele_it->get_region());
                haplotype_sub_sequence += the_reference_.get_sequence(reference_region_before_next_allele);
            } else {
                return the_allele_sequence == the_reference_.get_sequence(the_allele_region);
            }
        }
        
        if (ends_before(the_region_bounded_by_alleles, the_allele_region)) {
            haplotype_sub_sequence += get_sequence_bounded_by_alleles();
            haplotype_sub_sequence += the_reference_.get_sequence(reference_region_right_part);
        } else {
            auto last_allele_it = std::upper_bound(first_allele_it, std::cend(the_explicit_alleles_), the_allele_region);
            
            if (are_adjacent(the_allele_region, last_allele_it->get_region())) {
                haplotype_sub_sequence += get_sequence_bounded_by_alleles(first_allele_it, last_allele_it);
            } else {
                --last_allele_it;
                haplotype_sub_sequence += get_sequence_bounded_by_alleles(first_allele_it, last_allele_it);
                auto reference_region_after_last_allele = get_right_overhang(the_allele_region,
                                                                             last_allele_it->get_region());
                haplotype_sub_sequence += the_reference_.get_sequence(reference_region_after_last_allele);
            }
        }
        
        return haplotype_sub_sequence == the_allele_sequence;
    }
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
    if (the_explicit_alleles_.empty()) {
        return the_reference_.get_sequence(a_region);
    }
    
    auto the_region_bounded_by_alleles = get_region_bounded_by_alleles();
    
    SequenceType result {};
    
    auto reference_left_part  = get_left_overhang(a_region, the_region_bounded_by_alleles);
    if (size(reference_left_part) > 0) {
        result += the_reference_.get_sequence(reference_left_part);
    }
    
    result += get_sequence_bounded_by_alleles();
    
    auto reference_right_part = get_right_overhang(a_region, the_region_bounded_by_alleles);
    if (size(reference_right_part) > 0) {
        result += the_reference_.get_sequence(reference_right_part);
    }
    
    return result;
}

//void Haplotype::operator+=(const Haplotype& other)
//{
//    the_explicit_alleles_.insert(std::end(the_explicit_alleles_), std::cbegin(other.the_explicit_alleles_), std::cend(other.the_explicit_alleles_));
//}

GenomicRegion Haplotype::get_region_bounded_by_alleles() const
{
    if (the_explicit_alleles_.empty()) throw std::runtime_error {"Cannot get region from empty allele list"};
    
    return GenomicRegion {
        the_explicit_alleles_.front().get_region().get_contig_name(),
        the_explicit_alleles_.front().get_region().get_begin(),
        the_explicit_alleles_.back().get_region().get_end()
    };
}

Haplotype::SequenceType Haplotype::get_sequence_bounded_by_alleles(AlleleIterator first, AlleleIterator last) const
{
    SequenceType result {};
    
    if (first == last) return result;
    
    GenomicRegion previous_region {first->get_region()};
    
    std::for_each(first, last, [this, &previous_region, &result] (const auto& allele) {
        auto this_region = allele.get_region();
        
        if (previous_region != this_region && !overlaps(previous_region, this_region)) {
            auto intervening_reference_region   = get_intervening_region(previous_region, this_region);
            auto intervening_reference_sequence = the_reference_.get_sequence(intervening_reference_region);
            result += intervening_reference_sequence;
        }
        
        result += allele.get_sequence();
        
        previous_region = this_region;
    });
    
    return result;
}

Haplotype::SequenceType Haplotype::get_sequence_bounded_by_alleles() const
{
    return get_sequence_bounded_by_alleles(std::cbegin(the_explicit_alleles_), std::cend(the_explicit_alleles_));
}

void add_to_back(const Variant& a_variant, Haplotype& a_haplotype)
{
    a_haplotype.emplace_back(a_variant.get_region(), a_variant.get_alternative_allele_sequence());
}

void add_to_front(const Variant& a_variant, Haplotype& a_haplotype)
{
    a_haplotype.emplace_front(a_variant.get_region(), a_variant.get_alternative_allele_sequence());
}

bool contains(const Haplotype& a_haplotype, const Variant& a_variant)
{
    return a_haplotype.contains(a_variant.get_region(), a_variant.get_alternative_allele_sequence());
}
