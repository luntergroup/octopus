//
//  haplotype.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype.h"

#include <algorithm> // std::for_each, std::binary_search, std::equal_range, std::sort,
                     // std::nth_element, std::find_if_not, std::adjacent_find, std::unique
#include <iterator>  // std::cbegin, std::cend, std::distance, std::next

#include "reference_genome.h"
#include "genomic_region.h"
#include "region_algorithms.h"

Haplotype::Haplotype(ReferenceGenome& the_reference)
:
the_reference_ {&the_reference},
is_region_set_ {false},
the_reference_region_ {},
the_explicit_alleles_ {},
cached_sequence_ {},
is_cached_sequence_outdated_ {false}
{}

Haplotype::Haplotype(ReferenceGenome& the_reference, const GenomicRegion& a_region)
:
the_reference_ {&the_reference},
is_region_set_ {true},
the_reference_region_ {a_region},
the_explicit_alleles_ {},
cached_sequence_ {},
is_cached_sequence_outdated_ {true}
{}

bool Haplotype::contains(const Allele& an_allele) const
{
    if (the_explicit_alleles_.empty() && !is_region_set_) return false;
    
    if (::contains(get_region(), an_allele)) {
        // These binary searches are just optimisations and should not be required
        if (std::binary_search(std::cbegin(the_explicit_alleles_),
                               std::cend(the_explicit_alleles_), an_allele)) {
            return true;
        } else if (std::binary_search(std::cbegin(the_explicit_alleles_),
                                      std::cend(the_explicit_alleles_), an_allele.get_region())) {
            // If the allele is not explcitly contained but the region is then it must be a different
            // allele, unless it is an insertion, in which case we must check the sequence
            if (empty(an_allele.get_region())) {
                const auto& haplotype_allele = *std::equal_range(std::cbegin(the_explicit_alleles_),
                                                                 std::cend(the_explicit_alleles_),
                                                                 an_allele.get_region()).first;
                return ::contains(haplotype_allele, an_allele);
            } else {
                return false;
            }
        }
        
        auto overlapped_range = overlap_range(std::cbegin(the_explicit_alleles_),
                                              std::cend(the_explicit_alleles_), an_allele);
        if (std::distance(overlapped_range.begin(), overlapped_range.end()) == 1 &&
            ::contains(*overlapped_range.begin(), an_allele)) {
            return an_allele.get_sequence() ==
                    get_subsequence(*overlapped_range.begin(), get_overlapped(*overlapped_range.begin(), an_allele));
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
        cached_sequence_ = (is_region_set_) ? get_sequence(the_reference_region_) :
                                                get_sequence_bounded_by_explicit_alleles();
        is_cached_sequence_outdated_ = false;
        return cached_sequence_;
    }
}

Haplotype::SequenceType Haplotype::get_sequence(const GenomicRegion& a_region) const
{
    if (the_explicit_alleles_.empty()) {
        return the_reference_->get_sequence(a_region);
    }
    
    auto the_region_bounded_by_alleles = get_region_bounded_by_explicit_alleles();
    
    SequenceType result {};
    result.reserve(size(a_region));
    
    if (begins_before(a_region, the_region_bounded_by_alleles)) {
        result += the_reference_->get_sequence(get_left_overhang(a_region, the_region_bounded_by_alleles));
    }
    
    // we know the alleles are bidirectionally sorted as it is a condition of them being on a single haplotype
    auto overlapped_explicit_alleles = bases(overlap_range(std::cbegin(the_explicit_alleles_),
                                                           std::cend(the_explicit_alleles_), a_region,
                                                           MappableRangeOrder::BidirectionallySorted));
    
    if (overlapped_explicit_alleles.begin() != std::cend(the_explicit_alleles_)) {
        if (::contains(*overlapped_explicit_alleles.begin(), a_region)) {
            result += get_subsequence(*overlapped_explicit_alleles.begin(), a_region);
            return result;
        } else if (begins_before(*overlapped_explicit_alleles.begin(), a_region)) {
            result += get_subsequence(*overlapped_explicit_alleles.begin(),
                                      get_overlapped(*overlapped_explicit_alleles.begin(), a_region));
            overlapped_explicit_alleles.advance_begin(1);
        }
    }
    
    bool region_ends_before_last_overlapped_allele {false};
    
    if (!empty(overlapped_explicit_alleles) &&
        overlapped_explicit_alleles.end() != std::cend(the_explicit_alleles_) &&
        ends_before(a_region, *overlapped_explicit_alleles.end())) {
        overlapped_explicit_alleles.advance_end(-1);
        region_ends_before_last_overlapped_allele = true;
    }
    
    result += get_sequence_bounded_by_explicit_alleles(overlapped_explicit_alleles.begin(),
                                                       overlapped_explicit_alleles.end());
    
    if (region_ends_before_last_overlapped_allele) {
        result += get_subsequence(*overlapped_explicit_alleles.end(),
                                  get_overlapped(*overlapped_explicit_alleles.end(), a_region));
    } else if (ends_before(the_region_bounded_by_alleles, a_region)) {
        result += the_reference_->get_sequence(get_right_overhang(a_region, the_region_bounded_by_alleles));
    }
    
    result.shrink_to_fit();
    
    return result;
}

void Haplotype::operator+=(const Haplotype& other)
{
    the_explicit_alleles_.insert(std::end(the_explicit_alleles_), std::cbegin(other.the_explicit_alleles_),
                                 std::cend(other.the_explicit_alleles_));
}

GenomicRegion Haplotype::get_region_bounded_by_explicit_alleles() const
{
    if (the_explicit_alleles_.empty()) throw std::runtime_error {"Cannot get region from empty allele list"};
    
    return get_encompassing(the_explicit_alleles_.front(), the_explicit_alleles_.back());
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

// non-member methods

bool is_reference(const Haplotype& a_haplotype, ReferenceGenome& the_reference)
{
    return a_haplotype.get_sequence() == the_reference.get_sequence(a_haplotype.get_region());
}

void unique_least_complex(std::vector<Haplotype>& haplotypes)
{
    std::sort(haplotypes.begin(), haplotypes.end());
    
    auto first_equal = haplotypes.begin();
    auto last_equal  = haplotypes.begin();
    auto last        = haplotypes.end();
    
    while (true) {
        first_equal = std::adjacent_find(first_equal, last);
        
        if (first_equal == last) break;
        
        // skips 2 as std::next(first_equal, 1) is a duplicate
        last_equal = std::find_if_not(std::next(first_equal, 2), last,
                                      [first_equal] (const auto& haplotype) {
                                          return haplotype == *first_equal;
                                      });
        
        std::nth_element(first_equal, first_equal, last_equal, is_less_complex);
        
        first_equal = last_equal;
    }
    
    first_equal = std::unique(haplotypes.begin(), haplotypes.end());
    
    haplotypes.erase(first_equal, last);
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
