//
//  haplotype.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype.hpp"

#include <algorithm>
#include <iterator>

#include <boost/functional/hash.hpp>

#include "reference_genome.hpp"
#include "genomic_region.hpp"
#include "mappable_algorithms.hpp"
#include "variant.hpp"

#include <iostream> // TEST

// public methods

Haplotype::Haplotype(ReferenceGenome& reference, const GenomicRegion& region)
:
reference_ {reference},
region_ {region},
explicit_alleles_ {},
cached_sequence_ {}
{}

const GenomicRegion& Haplotype::get_region() const
{
    return region_;
}

bool Haplotype::contains(const Allele& allele) const
{
    using std::cbegin; using std::cend; using std::binary_search;
    
    if (::contains(get_region(), allele)) {
        // these binary searches are just optimisations
        if (binary_search(cbegin(explicit_alleles_), cend(explicit_alleles_), allele)) {
            return true;
        } else if (binary_search(cbegin(explicit_alleles_), cend(explicit_alleles_), ::get_region(allele))) {
            // If the allele is not explcitly contained but the region is then it must be a different
            // allele, unless it is an insertion, in which case we must check the sequence
            if (is_insertion(allele)) {
                auto it = std::lower_bound(cbegin(explicit_alleles_), cend(explicit_alleles_),
                                           allele.get_region());
                //return *it == allele;
                return ::contains(*it, allele);
            } else {
                return false;
            }
        }
        
        auto overlapped = bases(overlap_range(cbegin(explicit_alleles_), cend(explicit_alleles_), allele,
                                              MappableRangeOrder::BidirectionallySorted));
        
        if (overlapped.size() == 1 && ::contains(overlapped.front(), allele)) {
            return allele == splice(overlapped.front(), get_overlapped(overlapped.front(), allele));
        }
        
        return get_sequence(allele.get_region()) == allele.get_sequence();
    } else {
        return false;
    }
}

bool Haplotype::contains_exact(const Allele& allele) const
{
    return has_exact_overlap(std::cbegin(explicit_alleles_), std::cend(explicit_alleles_), allele,
                             MappableRangeOrder::BidirectionallySorted);
}

Haplotype::SequenceType Haplotype::get_sequence() const
{
    if (is_cached_sequence_good()) {
        return cached_sequence_;
    } else {
        cached_sequence_ = get_sequence(region_);
        return cached_sequence_;
    }
}

void append(Haplotype::SequenceType& sequence, const Allele& allele)
{
    sequence += allele.get_sequence();
}

Haplotype::SequenceType Haplotype::get_sequence(const GenomicRegion& region) const
{
    using std::cbegin; using std::cend;
    
    if (explicit_alleles_.empty()) {
        return get_reference_sequence(region);
    }
    
    auto region_bounded_by_alleles = get_region_bounded_by_explicit_alleles();
    
    SequenceType result {};
    result.reserve(size(region)); // may be more or less depending on indels
    
    if (begins_before(region, region_bounded_by_alleles)) {
        result += get_reference_sequence(get_left_overhang(region, region_bounded_by_alleles));
        if (is_before(region, region_bounded_by_alleles)) {
            if (!explicit_alleles_.empty() && ::contains(region, explicit_alleles_.front())) {
                result += explicit_alleles_.front().get_sequence();
            }
            result.shrink_to_fit();
            return result;
        }
    } else if (is_after(region, region_bounded_by_alleles)) {
        result += get_reference_sequence(region);
    }
    
    // we know the alleles are bidirectionally sorted as it is a condition of them being on a haplotype
    auto overlapped_explicit_alleles = bases(overlap_range(cbegin(explicit_alleles_), cend(explicit_alleles_), region,
                                                           MappableRangeOrder::BidirectionallySorted));
    
    // captures insertion at the end... this breaks the rules of region overlaps as normally
    // insertions are only consdiered overlapped at the start of a region.
//    if (overlapped_explicit_alleles.size() > 1 && is_insertion(overlapped_explicit_alleles.front())) {
//        
//    }
    if (overlapped_explicit_alleles.end() != cend(explicit_alleles_) && ::contains(region, *overlapped_explicit_alleles.end())) {
        overlapped_explicit_alleles.advance_end(1);
    }
    
    if (empty(overlapped_explicit_alleles)) { // can this ever happen?
        result.shrink_to_fit();
        return result;
    }
    
    if (::contains(overlapped_explicit_alleles.front(), region)) {
        append(result, splice(overlapped_explicit_alleles.front(), region));
        overlapped_explicit_alleles.advance_begin(1);
        if (!empty(overlapped_explicit_alleles) && is_insertion(overlapped_explicit_alleles.front())) {
            append(result, overlapped_explicit_alleles.front());
        }
        result.shrink_to_fit();
        return result;
    } else if (begins_before(overlapped_explicit_alleles.front(), region)) {
        append(result, splice(overlapped_explicit_alleles.front(), get_overlapped(overlapped_explicit_alleles.front(), region)));
        overlapped_explicit_alleles.advance_begin(1);
    }
    
    bool region_ends_before_last_overlapped_allele {ends_before(region, overlapped_explicit_alleles.back())};
    
    if (region_ends_before_last_overlapped_allele) {
        overlapped_explicit_alleles.advance_end(-1); // as we don't want all of the last allele
        region_ends_before_last_overlapped_allele = true;
    }
    
    result += get_sequence_bounded_by_explicit_alleles(overlapped_explicit_alleles.begin(),
                                                       overlapped_explicit_alleles.end());
    
    if (region_ends_before_last_overlapped_allele) {
        overlapped_explicit_alleles.advance_end(1); // as we previously removed this allele
        append(result, splice(overlapped_explicit_alleles.back(), get_overlapped(overlapped_explicit_alleles.back(), region)));
    } else if (ends_before(region_bounded_by_alleles, region)) {
        result += get_reference_sequence(get_right_overhang(region, region_bounded_by_alleles));
    }
    
    result.shrink_to_fit();
    
    return result;
}

std::vector<Variant> Haplotype::difference(const Haplotype& other) const
{
    std::vector<Variant> result {};
    result.reserve(explicit_alleles_.size());
    
    for (const auto& allele : explicit_alleles_) {
        if (!other.contains(allele)) {
            result.emplace_back(allele.get_region(), other.get_sequence(allele.get_region()), allele.get_sequence());
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

size_t Haplotype::get_hash() const
{
    using boost::hash_combine;
    
    if (cached_hash_ == 0) {
        size_t result {};
        
        hash_combine(result, std::hash<GenomicRegion>()(region_));
        hash_combine(result, std::hash<SequenceType>()(get_sequence()));
        
        if (result == 0) {
            ++result; // 0 is reserved
        }
        
        cached_hash_ = result;
    }
    
    return cached_hash_;
}

// private methods

Haplotype::SequenceType Haplotype::get_reference_sequence(const GenomicRegion& region) const
{
    return reference_.get().get_sequence(region);
}

Allele Haplotype::get_intervening_reference_allele(const Allele& lhs, const Allele& rhs) const
{
    return get_reference_allele(get_intervening(lhs, rhs), reference_);
}

GenomicRegion Haplotype::get_region_bounded_by_explicit_alleles() const
{
    if (explicit_alleles_.empty()) throw std::runtime_error {"cannot get region from empty allele list"};
    
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

bool Haplotype::is_cached_sequence_good() const noexcept
{
    return !cached_sequence_.empty() || (explicit_alleles_.empty() && empty(region_));
}

void Haplotype::reset_cached_sequence()
{
    cached_sequence_.clear();
}

// non-member methods

bool contains(const Haplotype& lhs, const Allele& rhs)
{
    return lhs.contains(rhs);
}

bool contains(const Haplotype& lhs, const Haplotype& rhs)
{
    if (!contains(get_region(lhs), get_region(rhs))) return false;
    
    auto rhs_explicit_allele_region = rhs.get_region_bounded_by_explicit_alleles();
    
    return lhs.get_sequence(rhs_explicit_allele_region) == rhs.get_sequence_bounded_by_explicit_alleles();
}

namespace detail {
Haplotype do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::true_type)
{
    // TODO: this is buggy.. test more
    
    if (get_region(haplotype) == region) return haplotype;
    
    Haplotype result {haplotype.reference_, region};
    
    if (haplotype.explicit_alleles_.empty()) return result;
    
    auto contained = bases(contained_range(std::cbegin(haplotype.explicit_alleles_),
                                           std::cend(haplotype.explicit_alleles_), region));
    
    switch (size(contained)) {
        case 0:
            break;
        case 1:
            result.push_back(contained.front());
            break;
        default:
            result.push_back(contained.front());
            result.explicit_alleles_.insert(std::end(result.explicit_alleles_), std::next(contained.begin()),
                                            std::prev(contained.end()));
            result.push_back(contained.back());
            break;
    }
    
    return result;
}

Allele do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::false_type)
{
    return Allele {region, haplotype.get_sequence(region)};
}
} // namespace detail

bool is_reference(const Haplotype& haplotype, ReferenceGenome& reference)
{
    return haplotype.get_sequence() == reference.get_sequence(haplotype.get_region());
}

bool IsLessComplex::operator()(const Haplotype& lhs, const Haplotype& rhs) noexcept
{
    return lhs.explicit_alleles_.size() < rhs.explicit_alleles_.size();
}

void unique_least_complex(std::vector<Haplotype>& haplotypes)
{
    using std::begin; using std::end;
    
    std::sort(begin(haplotypes), end(haplotypes));
    
    auto first_equal = begin(haplotypes);
    auto last_equal  = begin(haplotypes);
    auto last        = end(haplotypes);
    
    while (true) {
        first_equal = std::adjacent_find(first_equal, last);
        
        if (first_equal == last) break;
        
        // skips 2 as std::next(first_equal, 1) is a duplicate
        last_equal = std::find_if_not(std::next(first_equal, 2), last,
                                      [first_equal] (const auto& haplotype) {
                                          return haplotype == *first_equal;
                                      });
        
        std::nth_element(first_equal, first_equal, last_equal, IsLessComplex());
        
        first_equal = last_equal;
    }
    
    haplotypes.erase(std::unique(begin(haplotypes), end(haplotypes)), last);
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

bool HaveSameAlleles::operator()(const Haplotype &lhs, const Haplotype &rhs) const
{
    return lhs.explicit_alleles_ == rhs.explicit_alleles_;
}

bool have_same_alleles(const Haplotype& lhs, const Haplotype& rhs)
{
    return HaveSameAlleles()(lhs, rhs);
}

bool are_equal_in_region(const Haplotype& lhs, const Haplotype& rhs, const GenomicRegion& region)
{
    return splice<Haplotype>(lhs, region) == splice<Haplotype>(rhs, region);
}

std::ostream& operator<<(std::ostream& os, const Haplotype& haplotype)
{
    os << haplotype.get_region() << " " << haplotype.get_sequence();
    return os;
}

void print_alleles(const Haplotype& haplotype)
{
    std::cout << "< ";
    for (const auto& allele : haplotype.explicit_alleles_) {
        std::cout << "{" << allele << "} ";
    }
    std::cout << ">";
}

void print_variant_alleles(const Haplotype& haplotype)
{
    if (is_reference(haplotype, haplotype.reference_)) {
        std::cout << "< >";
    } else {
        std::cout << "< ";
        for (const auto& allele : haplotype.explicit_alleles_) {
            if (!is_reference(allele, haplotype.reference_)) std::cout << "{" << allele << "} ";
        }
        std::cout << ">";
    }
}
