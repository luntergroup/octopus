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
#include <stdexcept>
#include <cassert>

#include "reference_genome.hpp"
#include "genomic_region.hpp"
#include "mappable_ranges.hpp"
#include "mappable_algorithms.hpp"
#include "mappable_ranges.hpp"
#include "variant.hpp"

#include <iostream> // TEST

// helper

template <typename T, typename M>
auto haplotype_overlap_range(const T& alleles, const M& mappable)
{
    return bases(overlap_range(alleles, mappable, BidirectionallySortedTag {}));
}

template <typename T, typename M>
auto haplotype_contained_range(const T& alleles, const M& mappable)
{
    return bases(contained_range(alleles, mappable));
}

// public methods

const GenomicRegion& Haplotype::get_region() const
{
    return region_;
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
                }
                
                return false;
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

bool Haplotype::contains(const Allele& allele) const
{
    if (!is_same_contig(allele, region_)) return false;
    return contains(demote(allele));
}

bool Haplotype::contains_exact(const ContigAllele& allele) const
{
    if (!::contains(region_.get_contig_region(), allele)) return false;
    
    if (::contains(explicit_allele_region_, allele)) {
        return has_exact_overlap(explicit_alleles_, allele, BidirectionallySortedTag {});
    }
    
    if (overlaps(explicit_allele_region_, allele) || is_indel(allele)) return false;
    
    const auto it = std::next(std::cbegin(cached_sequence_),
                              begin_distance(allele, contig_region(region_)));
    
    return std::equal(std::cbegin(allele.get_sequence()), std::cend(allele.get_sequence()), it);
}

bool Haplotype::contains_exact(const Allele& allele) const
{
    if (!is_same_contig(allele, region_)) return false;
    return contains_exact(demote(allele));
}

bool is_in_reference_flank(const ContigRegion& region, const ContigRegion& explicit_allele_region_,
                           const std::vector<ContigAllele>& explicit_alleles)
{
    if (overlaps(region, explicit_allele_region_)) {
        return false;
    }
    
    if (!are_adjacent(region, explicit_allele_region_)) {
        return true;
    }
    
    if (begins_before(region, explicit_allele_region_)) {
        return !is_insertion(explicit_alleles.front());
    }
    
    return !is_insertion(explicit_alleles.back());
}

Haplotype::SequenceType Haplotype::get_sequence(const ContigRegion& region) const
{
    using std::cbegin; using std::cend;
    
    if (!::contains(region_.get_contig_region(), region)) {
        throw std::out_of_range {"Haplotype: attempting to get_sequence from region not contained by Haplotype region"};
    }
    
    if (explicit_alleles_.empty()) {
        return cached_sequence_.substr(begin_distance(region, region_.get_contig_region()),
                                       region_size(region));
    }
    
    if (is_in_reference_flank(region, explicit_allele_region_, explicit_alleles_)) {
        return get_reference_sequence(region);
    }
    
    SequenceType result {};
    result.reserve(region_size(region)); // may be more or less depending on indels
    
    if (begins_before(region, explicit_allele_region_)) {
        append_reference(result, left_overhang_region(region, explicit_allele_region_));
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
        
        if (overlapped_explicit_alleles.empty()) {
            append_reference(result, right_overhang_region(region, explicit_allele_region_));
            result.shrink_to_fit();
            return result;
        }
    }
    
    bool region_ends_before_last_overlapped_allele {ends_before(region, overlapped_explicit_alleles.back())};
    
    if (region_ends_before_last_overlapped_allele) {
        overlapped_explicit_alleles.advance_end(-1); // as we don't want all of the last allele
        region_ends_before_last_overlapped_allele = true;
    }
    
    append(result, cbegin(overlapped_explicit_alleles), cend(overlapped_explicit_alleles));
    
    if (region_ends_before_last_overlapped_allele) {
        overlapped_explicit_alleles.advance_end(1); // as we previously removed this allele
        append(result, splice(overlapped_explicit_alleles.back(),
                              overlapped_region(overlapped_explicit_alleles.back(), region)));
    } else if (ends_before(explicit_allele_region_, region)) {
        append_reference(result, right_overhang_region(region, explicit_allele_region_));
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

const Haplotype::SequenceType& Haplotype::get_sequence() const noexcept
{
    return cached_sequence_;
}

Haplotype::SizeType Haplotype::sequence_size(const ContigRegion& region) const
{
    return static_cast<SizeType>(get_sequence(region).size()); // TODO: can be improved
}

Haplotype::SizeType Haplotype::sequence_size(const GenomicRegion& region) const
{
    if (!is_same_contig(region, region_)) return 0;
    return sequence_size(region.get_contig_region());
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

std::size_t Haplotype::get_hash() const noexcept
{
    return cached_hash_;
}

// private methods

void Haplotype::append(SequenceType& result, const ContigAllele& allele) const
{
    result.append(allele.get_sequence());
}

void Haplotype::append(SequenceType& result, AlleleIterator first, AlleleIterator last) const
{
    std::for_each(first, last, [this, &result] (const auto& allele) { append(result, allele); });
}

void Haplotype::append_reference(SequenceType& result, const ContigRegion& region) const
{
    if (is_before(region, explicit_allele_region_)) {
        const auto offset = begin_distance(region, region_.get_contig_region());
        const auto it = std::next(std::cbegin(cached_sequence_), offset);
        result.append(it, std::next(it, region_size(region)));
    } else {
        const auto offset = end_distance(region_.get_contig_region(), region);
        const auto it = std::prev(std::cend(cached_sequence_), offset);
        result.append(std::prev(it, region_size(region)), it);
    }
}

Haplotype::SequenceType Haplotype::get_reference_sequence(const ContigRegion& region) const
{
    SequenceType result {};
    result.reserve(region_size(region));
    append_reference(result, region);
    return result;
}

// Builder

Haplotype::Builder::Builder(const GenomicRegion& region, const ReferenceGenome& reference)
:
region_ {region},
reference_ {reference}
{}

void Haplotype::Builder::push_back(const ContigAllele& allele)
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
}

void Haplotype::Builder::push_back(ContigAllele&& allele)
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
}

void Haplotype::Builder::push_front(const ContigAllele& allele)
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
}

void Haplotype::Builder::push_front(ContigAllele&& allele)
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
}

void Haplotype::Builder::push_back(const Allele& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::logic_error {"Haplotype::push_back called with Allele on different contig"};
    }
    push_back(ContigAllele {contig_region(allele), allele.get_sequence()});
}

void Haplotype::Builder::push_front(const Allele& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::logic_error {"Haplotype::push_front called with Allele on different contig"};
    }
    push_front(ContigAllele {contig_region(allele), allele.get_sequence()});
}

void Haplotype::Builder::push_back(Allele&& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::logic_error {"Haplotype::push_back called with Allele on different contig"};
    }
    push_back(ContigAllele {contig_region(allele), allele.get_sequence()});
}

void Haplotype::Builder::push_front(Allele&& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::logic_error {"Haplotype::push_front called with Allele on different contig"};
    }
    push_front(ContigAllele {contig_region(allele), allele.get_sequence()});
}

Haplotype Haplotype::Builder::build()
{
    return Haplotype {
        std::move(region_),
        std::make_move_iterator(std::begin(explicit_alleles_)),
        std::make_move_iterator(std::end(explicit_alleles_)),
        reference_
    };
}

void Haplotype::Builder::update_region(const ContigAllele& allele) noexcept
{
    const auto new_contig_region = encompassing_region(region_.get_contig_region(), allele);
    region_ = GenomicRegion {region_.get_contig_name(), new_contig_region};
}

void Haplotype::Builder::update_region(const Allele& allele)
{
    region_ = encompassing_region(region_, allele);
}

ContigAllele Haplotype::Builder::get_intervening_reference_allele(const ContigAllele& lhs, const ContigAllele& rhs) const
{
    const auto region = intervening_region(lhs, rhs);
    return ContigAllele {region, reference_.get().get_sequence(GenomicRegion {region_.get_contig_name(), region})};
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
    return contains(mapped_region(lhs), mapped_region(rhs))
            && lhs.get_sequence(rhs.region_) == rhs.get_sequence();
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
    
    Haplotype::Builder result {region, haplotype.reference_};
    
    if (haplotype.explicit_alleles_.empty()) return result.build();
    
    const auto& contig_region = region.get_contig_region();
    
    if (contains(contig_region, haplotype.explicit_allele_region_)) {
        result.explicit_alleles_.insert(end(result.explicit_alleles_),
                                        std::cbegin(haplotype.explicit_alleles_),
                                        std::cend(haplotype.explicit_alleles_));
        return result.build();
    }
    
    if (!overlaps(contig_region, haplotype.explicit_allele_region_)) return result.build();
    
    auto overlapped = haplotype_overlap_range(haplotype.explicit_alleles_, region.get_contig_region());
    
    assert(!overlapped.empty());
    
    if (is_empty(contig_region)) {
        if (!is_empty_region(overlapped.front()) && are_adjacent(contig_region, overlapped.front())) {
            overlapped.advance_begin(1);
            assert(!overlapped.empty());
        }
        
        if (is_empty_region(overlapped.front())) {
            result.push_back(overlapped.front());
        } else {
            result.push_back(ContigAllele {contig_region, ""});
        }
        
        return result.build();
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
    
    return result.build();
}

Allele do_splice(const Haplotype& haplotype, const GenomicRegion& region, std::false_type)
{
    return Allele {region, haplotype.get_sequence(region)};
}
} // namespace detail

ContigAllele splice(const Haplotype& haplotype, const ContigRegion& region)
{
    return ContigAllele {region, haplotype.get_sequence(region)};
}

bool is_reference(const Haplotype& haplotype)
{
    if (haplotype.explicit_alleles_.empty()) return true;
    return haplotype.get_sequence() == haplotype.reference_.get().get_sequence(haplotype.get_region());
}

std::vector<Variant> difference(const Haplotype& lhs, const Haplotype& rhs)
{
    auto result = lhs.difference(rhs);
    auto diffs2 = rhs.difference(lhs);
    
    const auto it = result.insert(std::end(result),
                                  std::make_move_iterator(std::begin(diffs2)),
                                  std::make_move_iterator(std::end(diffs2)));
    
    std::inplace_merge(std::begin(result), it, std::end(result));
    
    return result;
}

// Haplotype::Builder

void add_ref_to_back(const Variant& variant, Haplotype::Builder& haplotype)
{
    haplotype.push_back(variant.get_ref_allele());
}

void add_ref_to_front(const Variant& variant, Haplotype::Builder& haplotype)
{
    haplotype.push_front(variant.get_ref_allele());
}

void add_alt_to_back(const Variant& variant, Haplotype::Builder& haplotype)
{
    haplotype.push_back(variant.get_alt_allele());
}

void add_alt_to_front(const Variant& variant, Haplotype::Builder& haplotype)
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

bool IsLessComplex::operator()(const Haplotype& lhs, const Haplotype& rhs) noexcept
{
    return lhs.explicit_alleles_.size() < rhs.explicit_alleles_.size();
}

unsigned unique_least_complex(std::vector<Haplotype>& haplotypes)
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
    
    const auto it = std::unique(begin(haplotypes), end(haplotypes));
    
    const auto result = std::distance(it, end(haplotypes));
    
    haplotypes.erase(it, last);
    
    return static_cast<unsigned>(result);
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

namespace debug
{
    void print_alleles(const Haplotype& haplotype)
    {
        print_alleles(std::cout, haplotype);
    }
    
    void print_variant_alleles(const Haplotype& haplotype)
    {
        print_variant_alleles(std::cout, haplotype);
    }
    
    Haplotype make_haplotype(const std::string& str, const GenomicRegion& region,
                             const ReferenceGenome& reference)
    {
        if (str.size() < 3 || str.front() != '<' || str.back() != '>') {
            throw std::runtime_error {"make_haplotype: given bad input str"};
        }
        
        Haplotype::Builder hb {region, reference};
        
        if (str.size() == 3) {
            return hb.build(); // reference
        }
        
        std::size_t pos {3};
        
        while (pos < str.size()) {
            auto i = str.find(' ', pos);
            auto j = str.find('}', i + 1);
            hb.push_back(make_allele(str.substr(pos, i - pos), str.substr(i + 1, j - i - 1), reference));
            pos = j + 3;
        }
        
        return hb.build();
    }
    
    Haplotype make_haplotype(const std::string& str, const std::string& region,
                             const ReferenceGenome& reference)
    {
        return make_haplotype(str, parse_region(region, reference), reference);
    }
} // namespace debug
