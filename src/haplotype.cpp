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
#include "mappable_algorithms.h"

// public methods

Haplotype::Haplotype(ReferenceGenome& reference)
:
reference_ {&reference},
is_region_set_ {false},
reference_region_ {},
explicit_alleles_ {},
cached_sequence_ {},
is_cached_sequence_outdated_ {false}
{}

Haplotype::Haplotype(ReferenceGenome& reference, const GenomicRegion& region)
:
reference_ {&reference},
is_region_set_ {true},
reference_region_ {region},
explicit_alleles_ {},
cached_sequence_ {},
is_cached_sequence_outdated_ {true}
{}

bool Haplotype::contains(const Allele& allele) const
{
    if (explicit_alleles_.empty() && !is_region_set_) return false;
    
    if (::contains(get_region(), allele)) {
        // these binary searches are just optimisations
        if (std::binary_search(std::cbegin(explicit_alleles_), std::cend(explicit_alleles_), allele)) {
            return true;
        } else if (std::binary_search(std::cbegin(explicit_alleles_),
                                      std::cend(explicit_alleles_), allele.get_region())) {
            // If the allele is not explcitly contained but the region is then it must be a different
            // allele, unless it is an insertion, in which case we must check the sequence
            if (empty(allele.get_region())) {
                const auto& haplotype_allele = *std::equal_range(std::cbegin(explicit_alleles_),
                                                                 std::cend(explicit_alleles_),
                                                                 allele.get_region()).first;
                return ::contains(haplotype_allele, allele);
            } else {
                return false;
            }
        }
        
        auto overlapped_range = overlap_range(std::cbegin(explicit_alleles_),
                                              std::cend(explicit_alleles_), allele);
        if (std::distance(overlapped_range.begin(), overlapped_range.end()) == 1 &&
            ::contains(*overlapped_range.begin(), allele)) {
            return allele.get_sequence() ==
                    get_subsequence(*overlapped_range.begin(), get_overlapped(*overlapped_range.begin(), allele));
        }
        
        return get_sequence(allele.get_region()) == allele.get_sequence();
    } else {
        return false;
    }
}

void Haplotype::set_region(const GenomicRegion& region)
{
    reference_region_            = region;
    is_region_set_               = true;
    is_cached_sequence_outdated_ = true;
}

GenomicRegion Haplotype::get_region() const
{
    return (is_region_set_) ? reference_region_ : get_region_bounded_by_explicit_alleles();
}

Haplotype::SequenceType Haplotype::get_sequence() const
{
    if (!is_cached_sequence_outdated_) {
        return cached_sequence_;
    } else {
        cached_sequence_ = (is_region_set_) ? get_sequence(reference_region_) :
                                                get_sequence_bounded_by_explicit_alleles();
        is_cached_sequence_outdated_ = false;
        return cached_sequence_;
    }
}

Haplotype::SequenceType Haplotype::get_sequence(const GenomicRegion& region) const
{
    if (explicit_alleles_.empty()) {
        return reference_->get_sequence(region);
    }
    
    auto region_bounded_by_alleles = get_region_bounded_by_explicit_alleles();
    
    SequenceType result {};
    result.reserve(size(region));
    
    if (begins_before(region, region_bounded_by_alleles)) {
        result += reference_->get_sequence(get_left_overhang(region, region_bounded_by_alleles));
    }
    
    // we know the alleles are bidirectionally sorted as it is a condition of them being on a single haplotype
    auto overlapped_explicit_alleles = bases(overlap_range(std::cbegin(explicit_alleles_),
                                                           std::cend(explicit_alleles_), region,
                                                           MappableRangeOrder::BidirectionallySorted));
    
    if (overlapped_explicit_alleles.begin() != std::cend(explicit_alleles_)) {
        if (::contains(*overlapped_explicit_alleles.begin(), region)) {
            result += get_subsequence(*overlapped_explicit_alleles.begin(), region);
            return result;
        } else if (begins_before(*overlapped_explicit_alleles.begin(), region)) {
            result += get_subsequence(*overlapped_explicit_alleles.begin(),
                                      get_overlapped(*overlapped_explicit_alleles.begin(), region));
            overlapped_explicit_alleles.advance_begin(1);
        }
    }
    
    bool region_ends_before_last_overlapped_allele {false};
    
    if (!empty(overlapped_explicit_alleles) &&
        overlapped_explicit_alleles.end() != std::cend(explicit_alleles_) &&
        ends_before(region, *overlapped_explicit_alleles.end())) {
        overlapped_explicit_alleles.advance_end(-1);
        region_ends_before_last_overlapped_allele = true;
    }
    
    result += get_sequence_bounded_by_explicit_alleles(overlapped_explicit_alleles.begin(),
                                                       overlapped_explicit_alleles.end());
    
    if (region_ends_before_last_overlapped_allele) {
        result += get_subsequence(*overlapped_explicit_alleles.end(),
                                  get_overlapped(*overlapped_explicit_alleles.end(), region));
    } else if (ends_before(region_bounded_by_alleles, region)) {
        result += reference_->get_sequence(get_right_overhang(region, region_bounded_by_alleles));
    }
    
    result.shrink_to_fit();
    
    return result;
}

unsigned Haplotype::num_transitions() const noexcept
{
    return static_cast<unsigned>(std::count_if(std::cbegin(explicit_alleles_), std::cend(explicit_alleles_),
                                               [] (const auto& allele) {
                                                   return true;
                                               }));
}

unsigned Haplotype::num_transversions() const noexcept
{
    return static_cast<unsigned>(std::count_if(std::cbegin(explicit_alleles_), std::cend(explicit_alleles_),
                                               [] (const auto& allele) {
                                                   return true;
                                               }));
}

void Haplotype::operator+=(const Haplotype& other)
{
    explicit_alleles_.insert(std::end(explicit_alleles_), std::cbegin(other.explicit_alleles_),
                                 std::cend(other.explicit_alleles_));
}

// private methods

GenomicRegion Haplotype::get_region_bounded_by_explicit_alleles() const
{
    if (explicit_alleles_.empty()) throw std::runtime_error {"Cannot get region from empty allele list"};
    
    return get_encompassing(explicit_alleles_.front(), explicit_alleles_.back());
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
    return get_sequence_bounded_by_explicit_alleles(std::cbegin(explicit_alleles_),
                                                    std::cend(explicit_alleles_));
}

// non-member methods

bool contains(const Haplotype& lhs, const Haplotype& rhs)
{
    if (!contains(get_region(lhs), get_region(rhs))) return false;
    
    auto rhs_explicit_allele_region = rhs.get_region_bounded_by_explicit_alleles();
    
    return lhs.get_sequence(rhs_explicit_allele_region) == rhs.get_sequence_bounded_by_explicit_alleles();
}

Haplotype splice(const Haplotype& haplotype, const GenomicRegion& region)
{
    Haplotype result {*haplotype.reference_, region};
    
    auto contained = bases(contained_range(std::cbegin(haplotype.explicit_alleles_),
                                           std::cend(haplotype.explicit_alleles_), region));
    
    switch (size(contained)) {
        case 0:
            break;
        case 1:
            result.push_back(contained.front());
        default:
            result.push_back(contained.front());
            result.explicit_alleles_.insert(std::end(result.explicit_alleles_),
                                                std::next(contained.begin()), std::prev(contained.end()));
            result.push_back(contained.back());
            break;
    }
    
    return result;
}

bool is_reference(const Haplotype& haplotype, ReferenceGenome& reference)
{
    return haplotype.get_sequence() == reference.get_sequence(haplotype.get_region());
}

bool is_less_complex(const Haplotype& lhs, const Haplotype& rhs) noexcept
{
    return lhs.explicit_alleles_.size() < rhs.explicit_alleles_.size();
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
    
    haplotypes.erase(std::unique(haplotypes.begin(), haplotypes.end()), last);
}

bool operator==(const Haplotype& lhs, const Haplotype& rhs)
{
    return lhs.get_region() == rhs.get_region() && lhs.get_sequence() == rhs.get_sequence();
}

bool operator<(const Haplotype& lhs, const Haplotype& rhs)
{
    return (lhs.get_region() == rhs.get_region()) ? lhs.get_sequence() < rhs.get_sequence() :
    lhs.get_region() < rhs.get_region();
}

std::ostream& operator<<(std::ostream& os, const Haplotype& haplotype)
{
    os << haplotype.get_region() << " " << haplotype.get_sequence();
    return os;
}

void add_to_back(const Variant& a_variant, Haplotype& haplotype)
{
    haplotype.push_back(a_variant.get_alternative_allele());
}

void add_to_front(const Variant& a_variant, Haplotype& haplotype)
{
    haplotype.push_front(a_variant.get_alternative_allele());
}

bool contains(const Haplotype& haplotype, const Variant& a_variant)
{
    return haplotype.contains(a_variant.get_alternative_allele());
}
