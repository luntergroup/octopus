//
//  haplotype.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype.h"

#include <algorithm> // std::for_each
#include <iterator>  // std::cbegin etc

#include "reference_genome.h"
#include "genomic_region.h"
#include "region_utils.h"

#include <iostream> // TEST

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
    
    if (::contains(get_region(), an_allele)) {
        // These binary searches are just optimisations and should not be required
        if (std::binary_search(std::cbegin(the_explicit_alleles_),
                               std::cend(the_explicit_alleles_), an_allele)) {
            return true;
        } else if (std::binary_search(std::cbegin(the_explicit_alleles_),
                                      std::cend(the_explicit_alleles_), an_allele.get_region())) {
            // If the allele is not explcitly contained but the region is then it must be a different
            // allele
            return false;
        }
        
        return get_sequence(an_allele.get_region()) == an_allele.get_sequence();
    } else {
        return false;
    }
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
    if (the_explicit_alleles_.empty()) {
        return the_reference_.get_sequence(a_region);
    }
    
    auto the_region_bounded_by_alleles = get_region_bounded_by_explicit_alleles();
    
    SequenceType result {};
    
    if (begins_before(a_region, the_region_bounded_by_alleles)) {
        result += the_reference_.get_sequence(get_left_overhang(a_region, the_region_bounded_by_alleles));
    }
    
    auto overlapped_explicit_allele_range = overlap_range(std::cbegin(the_explicit_alleles_),
                                                          std::cend(the_explicit_alleles_), a_region);
    
    result += get_sequence_bounded_by_explicit_alleles(overlapped_explicit_allele_range.first,
                                                       overlapped_explicit_allele_range.second);
    
    if (ends_before(the_region_bounded_by_alleles, a_region)) {
        result += the_reference_.get_sequence(get_right_overhang(a_region, the_region_bounded_by_alleles));
    }
    
    return result;
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

Haplotype::SequenceType Haplotype::get_sequence_bounded_by_explicit_alleles(AlleleIterator first,
                                                                            AlleleIterator last) const
{
    SequenceType result {};
    
    std::for_each(first, last, [this, &result] (const auto& allele) {
        result += allele.get_sequence();
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
