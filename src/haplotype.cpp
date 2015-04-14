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
#include "region_utils.h"

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
the_explicit_alleles_ {},
cached_sequence_ {},
is_cached_sequence_outdated_ {false}
{}

Haplotype::Haplotype(ReferenceGenome& the_reference, const GenomicRegion& a_region)
:
the_reference_ {the_reference},
is_region_set_ {true},
the_reference_region_ {a_region},
the_explicit_alleles_ {},
cached_sequence_ {},
is_cached_sequence_outdated_ {true}
{}

bool Haplotype::contains(const Allele& an_allele) const
{
    if (the_explicit_alleles_.empty() && !is_region_set_) {
        return false;
    }
    
    return ::contains(get_region(), an_allele) && get_sequence(an_allele.get_region()) == an_allele.get_sequence();
}

void Haplotype::set_region(const GenomicRegion& a_region)
{
    the_reference_region_ = a_region;
    is_region_set_ = true;
    is_cached_sequence_outdated_ = true;
}

GenomicRegion Haplotype::get_region() const
{
    return (is_region_set_) ? the_reference_region_ : get_region_bounded_by_explicit_alleles();
}

Haplotype::SequenceType Haplotype::get_sequence() const
{
    if (!is_cached_sequence_outdated_) {
        return cached_sequence_;
    } else {
        cached_sequence_ = (is_region_set_) ? get_sequence(the_reference_region_) : get_sequence_bounded_by_explicit_alleles();
        is_cached_sequence_outdated_ = false;
        return cached_sequence_;
    }
}

Haplotype::SequenceType Haplotype::get_sequence(const GenomicRegion& a_region) const
{
    auto haplotype_region = get_region();
    
    SequenceType result {};
    
    if (begins_before(a_region, haplotype_region)) {
        result += the_reference_.get_sequence(get_left_overhang(a_region, haplotype_region));
    }
    
    if (!is_cached_sequence_outdated_) {
        auto first_cached_base = std::cbegin(cached_sequence_) + (get_begin(a_region) - get_begin(haplotype_region));
        result.insert(std::end(result), first_cached_base, first_cached_base +
                                        std::min(size(a_region), size(haplotype_region)));
        
    } else {
        auto the_region_bounded_by_alleles = get_region_bounded_by_explicit_alleles();
        
        if (begins_before(haplotype_region, the_region_bounded_by_alleles)) {
            result += the_reference_.get_sequence(get_left_overhang(haplotype_region, the_region_bounded_by_alleles));
        }
        
        auto overlapped_explicit_allele_range = overlap_range(std::cbegin(the_explicit_alleles_),
                                                              std::cend(the_explicit_alleles_), a_region);
        
        result += get_sequence_bounded_by_explicit_alleles(overlapped_explicit_allele_range.first,
                                                           overlapped_explicit_allele_range.second);
        
        if (ends_before(the_region_bounded_by_alleles, haplotype_region)) {
            result += the_reference_.get_sequence(get_right_overhang(haplotype_region, the_region_bounded_by_alleles));
        }
    }
    
    if (ends_before(haplotype_region, a_region)) {
        result += the_reference_.get_sequence(get_right_overhang(a_region, haplotype_region));
    }
    
    return result;
    
//    SequenceType result {};
//    
//    auto reference_left_part  = get_left_overhang(a_region, the_region_bounded_by_alleles);
//    if (size(reference_left_part) > 0) {
//        result += the_reference_.get_sequence(reference_left_part);
//    }
//    
//    result += get_sequence_bounded_by_explicit_alleles();
//    
//    auto reference_right_part = get_right_overhang(a_region, the_region_bounded_by_alleles);
//    if (size(reference_right_part) > 0) {
//        result += the_reference_.get_sequence(reference_right_part);
//    }
//    
//    return result;
}

//void Haplotype::operator+=(const Haplotype& other)
//{
//    the_explicit_alleles_.insert(std::end(the_explicit_alleles_), std::cbegin(other.the_explicit_alleles_), std::cend(other.the_explicit_alleles_));
//}

GenomicRegion Haplotype::get_region_bounded_by_explicit_alleles() const
{
    if (the_explicit_alleles_.empty()) throw std::runtime_error {"Cannot get region from empty allele list"};
    
    return GenomicRegion {
        the_explicit_alleles_.front().get_region().get_contig_name(),
        the_explicit_alleles_.front().get_region().get_begin(),
        the_explicit_alleles_.back().get_region().get_end()
    };
}

Haplotype::SequenceType Haplotype::get_sequence_bounded_by_explicit_alleles(AlleleIterator first, AlleleIterator last) const
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

Haplotype::SequenceType Haplotype::get_sequence_bounded_by_explicit_alleles() const
{
    return get_sequence_bounded_by_explicit_alleles(std::cbegin(the_explicit_alleles_),
                                                    std::cend(the_explicit_alleles_));
}

void add_to_back(const Variant& a_variant, Haplotype& a_haplotype)
{
    a_haplotype.push_back(a_variant.get_alternative_allele());
}

void add_to_front(const Variant& a_variant, Haplotype& a_haplotype)
{
    a_haplotype.push_front(a_variant.get_alternative_allele());
}

bool contains(const Haplotype& a_haplotype, const Variant& a_variant)
{
    return a_haplotype.contains(a_variant.get_alternative_allele());
}
