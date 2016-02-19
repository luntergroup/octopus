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
#include <utility>
#include <stdexcept>

#include <boost/functional/hash.hpp>

#include "reference_genome.hpp"
#include "genomic_region.hpp"
#include "mappable_algorithms.hpp"
#include "mappable_ranges.hpp"
#include "variant.hpp"

#include <iostream> // TEST

// helper

template <typename T, typename M>
auto haplotype_overlap_range(const T& alleles, const M& mappable)
{
    return bases(overlap_range(std::cbegin(alleles), std::cend(alleles),
                               mappable, MappableRangeOrder::BidirectionallySorted));
}

template <typename T, typename M>
auto haplotype_contained_range(const T& alleles, const M& mappable)
{
    return bases(contained_range(std::cbegin(alleles), std::cend(alleles), mappable));
}

// public methods

Haplotype::Haplotype(const GenomicRegion& region, const ReferenceGenome& reference)
:
region_ {region},
explicit_alleles_ {},
cached_sequence_ {},
reference_ {reference}
{}

const GenomicRegion& Haplotype::get_region() const
{
    return region_;
}

void Haplotype::push_back(const ContigAllele& allele)
{
    if (!explicit_alleles_.empty()) {
        if (!is_after(allele, explicit_alleles_.back())) {
            throw std::logic_error {"Haplotype::push_back called with out-of-order Allele"};
        }
        if (!are_adjacent(explicit_alleles_.back(), allele)) {
            explicit_alleles_.emplace_back(get_intervening_reference_allele(explicit_alleles_.back(), allele));
        }
    }
    update_region(allele);
    explicit_alleles_.emplace_back(allele);
    clear_cached_sequence();
}

void Haplotype::push_back(ContigAllele&& allele)
{
    if (!explicit_alleles_.empty()) {
        if (!is_after(allele, explicit_alleles_.back())) {
            throw std::logic_error {"Haplotype::push_back called with out-of-order Allele"};
        }
        if (!are_adjacent(explicit_alleles_.back(), allele)) {
            explicit_alleles_.emplace_back(get_intervening_reference_allele(explicit_alleles_.back(), allele));
        }
    }
    update_region(allele);
    explicit_alleles_.emplace_back(std::move(allele));
    clear_cached_sequence();
}

void Haplotype::push_front(const ContigAllele& allele)
{
    if (!explicit_alleles_.empty()) {
        if (!is_after(explicit_alleles_.front(), allele)) {
            throw std::logic_error {"Haplotype::push_front called with out-of-order Allele"};
        }
        if (!are_adjacent(allele, explicit_alleles_.front())) {
            explicit_alleles_.emplace_front(get_intervening_reference_allele(allele, explicit_alleles_.front()));
        }
    }
    update_region(allele);
    explicit_alleles_.emplace_front(allele);
    clear_cached_sequence();
}

void Haplotype::push_front(ContigAllele&& allele)
{
    if (!explicit_alleles_.empty()) {
        if (!is_after(explicit_alleles_.front(), allele)) {
            throw std::logic_error {"Haplotype::push_front called with out-of-order Allele"};
        }
        if (!are_adjacent(allele, explicit_alleles_.front())) {
            explicit_alleles_.emplace_front(get_intervening_reference_allele(allele, explicit_alleles_.front()));
        }
    }
    update_region(allele);
    explicit_alleles_.emplace_front(std::move(allele));
    clear_cached_sequence();
}

void Haplotype::push_back(const Allele& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::logic_error {"Haplotype::push_back called with Allele on different contig"};
    }
    push_back(ContigAllele {contig_region(allele), allele.get_sequence()});
}

void Haplotype::push_front(const Allele& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::logic_error {"Haplotype::push_front called with Allele on different contig"};
    }
    push_front(ContigAllele {contig_region(allele), allele.get_sequence()});
}

void Haplotype::push_back(Allele&& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::logic_error {"Haplotype::push_back called with Allele on different contig"};
    }
    push_back(ContigAllele {contig_region(allele), allele.get_sequence()});
}

void Haplotype::push_front(Allele&& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::logic_error {"Haplotype::push_front called with Allele on different contig"};
    }
    push_front(ContigAllele {contig_region(allele), allele.get_sequence()});
}

bool Haplotype::contains(const ContigAllele& allele) const
{
    using std::cbegin; using std::cend;
    
    if (::contains(region_.get_contig_region(), allele)) {
        const auto it = std::lower_bound(cbegin(explicit_alleles_), cend(explicit_alleles_),
                                         allele.get_region());
        
        if (it != cend(explicit_alleles_)) {
            if (*it == allele) return true;
            
            if (is_same_region(*it, allele)) {
                // If the allele is not explcitly contained but the region is then it must be a different
                // allele, unless it is an insertion, in which case we must check the sequence
                if (is_insertion(allele)) {
                    const auto it = std::lower_bound(cbegin(explicit_alleles_), cend(explicit_alleles_),
                                                     allele.get_region());
                    return ::contains(*it, allele);
                } else {
                    return false;
                }
            }
        }
        
        const auto overlapped = haplotype_overlap_range(explicit_alleles_, allele);
        
        if (overlapped.size() == 1 && ::contains(overlapped.front(), allele)) {
            return allele == splice(overlapped.front(), contig_region(allele));
        }
        
        return get_sequence(allele.get_region()) == allele.get_sequence();
    }
    
    return false;
}

bool Haplotype::contains_exact(const ContigAllele& allele) const
{
    return has_exact_overlap(explicit_alleles_, allele, MappableRangeOrder::BidirectionallySorted);
}

bool Haplotype::contains(const Allele& allele) const
{
    if (!is_same_contig(allele, region_)) return false;
    return contains(demote(allele));
}

bool Haplotype::contains_exact(const Allele& allele) const
{
    if (!is_same_contig(allele, region_)) return false;
    return contains_exact(demote(allele));
}

void append(Haplotype::SequenceType& sequence, const ContigAllele& allele)
{
    sequence += allele.get_sequence();
}

Haplotype::SequenceType Haplotype::get_sequence(const ContigRegion& region) const
{
    using std::cbegin; using std::cend;
    
    if (!::contains(region_.get_contig_region(), region)) {
        throw std::out_of_range {"Haplotype: attempting to get_sequence from region not contained by Haplotype region"};
    }
    
    if (explicit_alleles_.empty()) {
        if (is_cached_sequence_good()) {
            return cached_sequence_.substr(begin_distance(region, region_.get_contig_region()),
                                           region_size(region));
        } else {
            return get_reference_sequence(region);
        }
    }
    
    const auto region_bounded_by_alleles = get_region_bounded_by_explicit_alleles();
    
    if (is_before(region, region_bounded_by_alleles) || is_after(region, region_bounded_by_alleles)) {
        return get_reference_sequence(region);
    }
    
    SequenceType result {};
    result.reserve(region_size(region)); // may be more or less depending on indels
    
    if (begins_before(region, region_bounded_by_alleles)) {
        result += get_reference_sequence(left_overhang_region(region, region_bounded_by_alleles));
    }
    
    auto overlapped_explicit_alleles = haplotype_overlap_range(explicit_alleles_, region);
    
    // known that !overlapped_explicit_alleles.empty()
    
    if (::contains(overlapped_explicit_alleles.front(), region)) {
        append(result, splice(overlapped_explicit_alleles.front(), region));
        overlapped_explicit_alleles.advance_begin(1);
        if (!overlapped_explicit_alleles.empty() && is_insertion(overlapped_explicit_alleles.front())) {
            append(result, overlapped_explicit_alleles.front());
        }
        result.shrink_to_fit();
        return result;
    } else if (begins_before(overlapped_explicit_alleles.front(), region)) {
        append(result, splice(overlapped_explicit_alleles.front(),
                              overlapped_region(overlapped_explicit_alleles.front(), region)));
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
        append(result, splice(overlapped_explicit_alleles.back(),
                              overlapped_region(overlapped_explicit_alleles.back(), region)));
    } else if (ends_before(region_bounded_by_alleles, region)) {
        result += get_reference_sequence(right_overhang_region(region, region_bounded_by_alleles));
    }
    
    result.shrink_to_fit();
    
    return result;
}

Haplotype::SequenceType Haplotype::get_sequence(const GenomicRegion& region) const
{
    if (!is_same_contig(region, region_)) {
        throw std::logic_error {"Haplotype: cannot get_sequence from different contig"};
    }
    return get_sequence(region.get_contig_region());
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

std::vector<Variant> Haplotype::difference(const Haplotype& other) const
{
    std::vector<Variant> result {};
    result.reserve(explicit_alleles_.size());
    
    const auto& contig = region_.get_contig_name();
    
    for (const auto& allele : explicit_alleles_) {
        if (!other.contains(allele)) {
            result.emplace_back(GenomicRegion {contig, allele.get_region()},
                                other.get_sequence(allele.get_region()),
                                allele.get_sequence());
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

size_t Haplotype::get_hash() const
{
    using boost::hash_combine;
    
    if (cached_hash_ == 0) {
        size_t result {0};
        
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

Haplotype::SequenceType Haplotype::get_reference_sequence(const ContigRegion& region) const
{
    return get_reference_sequence(GenomicRegion {region_.get_contig_name(), region});
}

ContigAllele Haplotype::get_intervening_reference_allele(const ContigAllele& lhs, const ContigAllele& rhs) const
{
    const auto region = intervening_region(lhs, rhs);
    return ContigAllele {region, reference_.get().get_sequence(GenomicRegion {region_.get_contig_name(), region})};
}

ContigRegion Haplotype::get_region_bounded_by_explicit_alleles() const
{
    if (explicit_alleles_.empty()) {
        throw std::runtime_error {"Haplotype: trying to get region from empty allele list"};
    }
    return encompassing_region(explicit_alleles_.front(), explicit_alleles_.back());
}

Haplotype::SequenceType Haplotype::get_sequence_bounded_by_explicit_alleles(AlleleIterator first,
                                                                            AlleleIterator last) const
{
    SequenceType result {};
    result.reserve(64);
    
    std::for_each(first, last, [this, &result] (const auto& allele) {
        result += allele.get_sequence();
    });
    
    result.shrink_to_fit();
    
    return result;
}

Haplotype::SequenceType Haplotype::get_sequence_bounded_by_explicit_alleles() const
{
    return get_sequence_bounded_by_explicit_alleles(std::cbegin(explicit_alleles_),
                                                    std::cend(explicit_alleles_));
}

void Haplotype::update_region(const ContigAllele& allele) noexcept
{
    const auto new_contig_region = encompassing_region(region_.get_contig_region(), allele);
    region_ = GenomicRegion {region_.get_contig_name(), new_contig_region};
}

void Haplotype::update_region(const Allele& allele)
{
    region_ = encompassing_region(region_, allele);
}

bool Haplotype::is_cached_sequence_good() const noexcept
{
    return !cached_sequence_.empty() || (explicit_alleles_.empty() && is_empty_region(region_));
}

void Haplotype::clear_cached_sequence()
{
    cached_sequence_.clear();
}

// non-member methods

Haplotype::SizeType sequence_size(const Haplotype& haplotype) noexcept
{
    return static_cast<Haplotype::SizeType>(haplotype.get_sequence().size());
}

bool is_empty_sequence(const Haplotype& haplotype) noexcept
{
    return haplotype.get_sequence().empty();
}

bool contains(const Haplotype& lhs, const Allele& rhs)
{
    return lhs.contains(rhs);
}

bool contains(const Haplotype& lhs, const Haplotype& rhs)
{
    if (!contains(mapped_region(lhs), mapped_region(rhs))) return false;
    
    auto rhs_explicit_allele_region = rhs.get_region_bounded_by_explicit_alleles();
    
    return lhs.get_sequence(rhs_explicit_allele_region) == rhs.get_sequence_bounded_by_explicit_alleles();
}

namespace detail
{
Haplotype do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::true_type)
{
    using std::end; using std::cbegin; using std::cend; using std::prev;
    
    if (!is_same_contig(haplotype, region)) {
        throw std::logic_error {"Haplotype: trying to splice region from different contig"};
    }
    
    if (is_same_region(haplotype, region)) return haplotype;
    
    Haplotype result {region, haplotype.reference_};
    
    if (haplotype.explicit_alleles_.empty()) return result;
    
    const auto& contig_region = region.get_contig_region();
    
    const auto explicit_allele_region = haplotype.get_region_bounded_by_explicit_alleles();
    
    if (contains(contig_region, explicit_allele_region)) {
        result.explicit_alleles_ = haplotype.explicit_alleles_;
        return result;
    }
    
    if (!overlaps(contig_region, explicit_allele_region)) return result;
    
    auto overlapped = haplotype_overlap_range(haplotype.explicit_alleles_, region.get_contig_region());
    
    // known that !overlapped.empty()
    
    if (is_empty_region(contig_region)) {
        if (!is_empty_region(overlapped.front())) {
            overlapped.advance_begin(1);
        }
        
        if (is_empty_region(overlapped.front())) {
            result.push_back(overlapped.front());
        } else {
            result.push_back(ContigAllele {contig_region, ""});
        }
        
        return result;
    }
    
    if (!contains(contig_region, overlapped.front())) {
        result.push_front(splice(overlapped.front(), overlapped_region(overlapped.front(), contig_region)));
        overlapped.advance_begin(1);
    }
    
    if (!overlapped.empty()) {
        if (contains(contig_region, overlapped.back())) {
            result.explicit_alleles_.insert(end(result.explicit_alleles_),
                                            cbegin(overlapped), cend(overlapped));
        } else {
            result.explicit_alleles_.insert(end(result.explicit_alleles_),
                                            cbegin(overlapped), prev(cend(overlapped)));
            result.push_back(splice(overlapped.back(), overlapped_region(overlapped.back(), contig_region)));
        }
    }
    
    return result;
}

Allele do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::false_type)
{
    return Allele {region, haplotype.get_sequence(region)};
}
} // namespace detail

bool is_reference(const Haplotype& haplotype)
{
    if (haplotype.explicit_alleles_.empty()) return true;
    return haplotype.get_sequence() == haplotype.reference_.get().get_sequence(haplotype.get_region());
}

void add_ref_to_back(const Variant& variant, Haplotype& haplotype)
{
    haplotype.push_back(variant.get_ref_allele());
}

void add_ref_to_front(const Variant& variant, Haplotype& haplotype)
{
    haplotype.push_front(variant.get_ref_allele());
}

void add_alt_to_back(const Variant& variant, Haplotype& haplotype)
{
    haplotype.push_back(variant.get_alt_allele());
}

void add_alt_to_front(const Variant& variant, Haplotype& haplotype)
{
    haplotype.push_front(variant.get_alt_allele());
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

// debug

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
    if (is_reference(haplotype)) {
        std::cout << "< >";
    } else {
        const auto& contig = contig_name(haplotype);
        std::cout << "< ";
        for (const auto& contig_allele : haplotype.explicit_alleles_) {
            Allele allele {GenomicRegion {contig, contig_allele.get_region()}, contig_allele.get_sequence()};
            if (!is_reference(allele, haplotype.reference_)) std::cout << "{" << allele << "} ";
        }
        std::cout << ">";
    }
}
