//
//  variant.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant.hpp"

#include <list>
#include <utility>
#include <algorithm>
#include <iterator>
#include <cassert>

#include <boost/range/algorithm.hpp>

#include "reference_genome.hpp"
#include "candidate_variant_generator.hpp"
#include "mappable_algorithms.hpp"

// public methods

Variant::Variant(const Allele& reference, const Allele& alternative)
:
reference_ {reference},
alternative_ {alternative}
{}

Variant::Variant(Allele&& reference, Allele&& alternative)
:
reference_ {std::move(reference)},
alternative_ {std::move(alternative)}
{}

const GenomicRegion& Variant::get_region() const noexcept
{
    return reference_.get_region();
}

const Allele& Variant::get_ref_allele() const noexcept
{
    return reference_;
}

const Allele& Variant::get_alt_allele() const noexcept
{
    return alternative_;
}

// non-member methods

Variant make_variant(const Allele& alt_allele, const ReferenceGenome& reference)
{
    return Variant {make_reference_allele(mapped_region(alt_allele), reference), alt_allele};
}

Variant make_variant(const std::string& region_str, Variant::SequenceType alt_sequence,
                     const ReferenceGenome& reference)
{
    return make_variant(Allele {parse_region(region_str, reference), std::move(alt_sequence)}, reference);
}

const Variant::SequenceType& ref_sequence(const Variant& variant)
{
    return variant.get_ref_allele().get_sequence();
}

const Variant::SequenceType& alt_sequence(const Variant& variant)
{
    return variant.get_alt_allele().get_sequence();
}

Variant::SizeType ref_sequence_size(const Variant& variant)
{
    return sequence_size(variant.get_ref_allele());
}

Variant::SizeType alt_sequence_size(const Variant& variant)
{
    return sequence_size(variant.get_alt_allele());
}

bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_ref_allele() == rhs.get_ref_allele() && alt_sequence(lhs) == alt_sequence(rhs);
}

bool operator<(const Variant& lhs, const Variant& rhs)
{
    return (lhs.get_ref_allele() < rhs.get_ref_allele() ||
            (lhs.get_ref_allele() == rhs.get_ref_allele() && lhs.get_alt_allele() < rhs.get_alt_allele()));
}

void remove_duplicates(std::vector<Variant>& variants)
{
    variants.erase(std::unique(std::begin(variants), std::end(variants)), std::end(variants));
}

std::vector<Allele> decompose(const std::vector<Variant>& variants)
{
    std::vector<Allele> result {};
    result.reserve(2 * variants.size()); // max num alleles (may be less)
    
    if (variants.empty()) return result;
    
    for (const auto& variant : variants) {
        if (result.empty() || !is_same_region(result.back(), variant)) {
            result.emplace_back(variant.get_ref_allele());
        }
        result.emplace_back(variant.get_alt_allele());
    }
    
    result.shrink_to_fit();
    
    return result;
}

std::vector<Allele> extract_intervening_reference_alleles(const std::vector<Variant>& variants,
                                                          const ReferenceGenome& reference)
{
    const auto all_overlapped_regions     = extract_covered_regions(variants);
    const auto regions_between_candidates = extract_intervening_regions(all_overlapped_regions);
    
    std::vector<Allele> result {};
    result.reserve(regions_between_candidates.size());
    
    boost::transform(regions_between_candidates, std::back_inserter(result),
                     [&reference] (const auto& region) {
                         return make_reference_allele(region, reference);
                     });
    
    return result;
}

auto allele_minmax(const Variant::SequenceType& lhs, const Variant::SequenceType& rhs)
{
    return std::minmax(lhs, rhs, [] (const auto& lhs, const auto& rhs) { return lhs.size() < rhs.size(); });
}

template <typename InputIterator>
size_t get_num_redundant_bases(InputIterator first1, InputIterator last1, InputIterator first2)
{
    return std::distance(first1, std::mismatch(first1, last1, first2).first);
}

template <typename SequenceType>
size_t get_num_redundant_front_bases(const SequenceType& smaller_sequence,
                                     const SequenceType& larger_sequence)
{
    return get_num_redundant_bases(std::cbegin(smaller_sequence), std::cend(smaller_sequence),
                                   std::cbegin(larger_sequence));
}

template <typename SequenceType>
size_t get_num_redundant_back_bases(const SequenceType& smaller_sequence,
                                    const SequenceType& larger_sequence)
{
    return get_num_redundant_bases(std::crbegin(smaller_sequence), std::crend(smaller_sequence),
                                   std::crbegin(larger_sequence));
}

bool is_parsimonious(const Variant& variant) noexcept
{
    using std::cbegin; using std::cend; using std::crbegin; using std::crend;
    
    if (ref_sequence_size(variant) == 0 || alt_sequence_size(variant) == 0) {
        return false;
    }
    
    const auto& alleles = allele_minmax(ref_sequence(variant), alt_sequence(variant));
    
    const auto& small_allele = alleles.first;
    const auto& big_allele   = alleles.second;
    
    if (get_num_redundant_front_bases(small_allele, big_allele) > 1) {
        return false;
    }
    
    if (small_allele.size() > 1 && get_num_redundant_back_bases(small_allele, big_allele) > 0) {
        return false;
    }
    
    return true;
}

Variant make_parsimonious(const Variant& variant, const ReferenceGenome& reference)
{
    using std::move; using std::cbegin; using std::cend; using std::crbegin; using std::crend;
    
    if (is_parsimonious(variant)) return variant;
    
    const auto& old_ref_sequence = ref_sequence(variant);
    const auto& old_alt_sequence = alt_sequence(variant);
    
    const auto& alleles = allele_minmax(old_ref_sequence, old_alt_sequence);
    const auto& the_small_allele = alleles.first;
    const auto& the_big_allele   = alleles.second;
    
    const auto& old_ref_region = mapped_region(variant);
    
    if (the_small_allele.size() == 0) {
        if (old_ref_region.get_begin() > 0) {
            GenomicRegion new_ref_region {old_ref_region.get_contig_name(),
                old_ref_region.get_begin() - 1, old_ref_region.get_end()};
            
            auto new_ref_allele = reference.get_sequence(new_ref_region);
            auto new_alt_allele = new_ref_allele.front() + old_alt_sequence;
            
            return Variant {move(new_ref_region), move(new_ref_allele), move(new_alt_allele)};
        } else {
            // In this very rare care the only option is to pad to the right.
            GenomicRegion new_ref_region {old_ref_region.get_contig_name(),
                old_ref_region.get_begin(), old_ref_region.get_end() + 1};
            
            auto new_ref_allele = reference.get_sequence(new_ref_region);
            auto new_alt_allele = old_alt_sequence + new_ref_allele.back();
            
            return Variant {move(new_ref_region), move(new_ref_allele), move(new_alt_allele)};
        }
    }
    
    auto num_redundant_front_bases = get_num_redundant_front_bases(the_small_allele, the_big_allele);
    auto num_redundant_back_bases  = get_num_redundant_back_bases(the_small_allele, the_big_allele);
    
    // We could avoid this check by removing redundant back bases first, but this way is cheaper.
    bool are_same_redundant_bases = num_redundant_front_bases == num_redundant_back_bases &&
    num_redundant_back_bases == the_small_allele.size();
    if (are_same_redundant_bases) num_redundant_back_bases = 0;
    // Don't forget: a parsimonious variant can't have any empty alleles
    if (num_redundant_front_bases == the_small_allele.size()) --num_redundant_front_bases;
    
    auto new_big_allele_size   = the_big_allele.size() - num_redundant_front_bases -
    num_redundant_back_bases;
    auto new_small_allele_size = the_small_allele.size() - num_redundant_front_bases -
    num_redundant_back_bases;
    
    auto new_big_allele   = the_big_allele.substr(num_redundant_front_bases, new_big_allele_size);
    auto new_small_allele = the_small_allele.substr(num_redundant_front_bases, new_small_allele_size);
    
    using SizeType = GenomicRegion::SizeType;
    
    auto new_ref_region_begin = static_cast<SizeType>(old_ref_region.get_begin() +
                                                      num_redundant_front_bases);
    auto new_ref_region_end   = static_cast<SizeType>(old_ref_region.get_end() -
                                                      num_redundant_back_bases);
    
    GenomicRegion new_ref_region {old_ref_region.get_contig_name(), new_ref_region_begin,
        new_ref_region_end};
    
    if (old_ref_sequence.size() > old_alt_sequence.size()) {
        return Variant {move(new_ref_region), move(new_big_allele), move(new_small_allele)};
    } else {
        return Variant {move(new_ref_region), move(new_small_allele), move(new_big_allele)};
    }
}

bool is_left_alignable(const Variant& variant) noexcept
{
    return is_indel(variant);
}

using LeftAlignmentList = std::list<Variant::SequenceType::value_type>;

auto get_allele_lists(const Variant::SequenceType& allele_a, const Variant::SequenceType& allele_b)
{
    const auto& alleles = allele_minmax(allele_a, allele_b);
    
    LeftAlignmentList big_allele {std::cbegin(alleles.second), std::cend(alleles.second)};
    LeftAlignmentList small_allele {std::cbegin(alleles.first), std::cend(alleles.first)};
    
    return std::make_pair(std::move(big_allele), std::move(small_allele));
}

GenomicRegion extend_allele_lists(LeftAlignmentList& big_allele, LeftAlignmentList& small_allele,
                                  const ReferenceGenome& reference, const GenomicRegion& current_region,
                                  const Variant::SizeType extension_size)
{
    using std::begin; using std::cbegin; using std::cend;
    
    const auto new_region = shift(current_region, -extension_size);
    
    GenomicRegion extension_region {new_region.get_contig_name(), new_region.get_begin(),
        new_region.get_begin() + extension_size};
    
    const auto extension_sequence = reference.get_sequence(extension_region);
    
    big_allele.insert(begin(big_allele), cbegin(extension_sequence), cend(extension_sequence));
    small_allele.insert(begin(small_allele), cbegin(extension_sequence), cend(extension_sequence));
    
    return new_region;
}

Variant left_align(const Variant& variant, const ReferenceGenome& reference,
                   const Variant::SizeType extension_size)
{
    using std::move; using std::cbegin; using std::cend; using std::crbegin; using std::crend;
    using std::tie; using std::next; using std::mismatch;
    
    if (!is_left_alignable(variant)) return variant;
    
    static_assert(sizeof(Variant::SizeType) == sizeof(GenomicRegion::SizeType),
                  "Variant and GenomicRegion have different SizeType");
    
    using SizeType = Variant::SizeType;
    
    const auto& ref_allele_sequence = ref_sequence(variant);
    const auto& alt_allele_sequence = alt_sequence(variant);
    
    LeftAlignmentList big_allele {}, small_allele {};
    
    tie(big_allele, small_allele) = get_allele_lists(ref_allele_sequence, alt_allele_sequence);
    
    auto big_allele_ritr   = crbegin(big_allele);
    auto small_allele_ritr = crbegin(small_allele);
    
    const auto big_allele_size   = static_cast<SizeType>(big_allele.size());
    const auto small_allele_size = static_cast<SizeType>(small_allele.size());
    
    GenomicRegion current_region {variant.get_region()};
    
    do {
        if (current_region.get_begin() >= extension_size) {
            current_region = extend_allele_lists(big_allele, small_allele, reference,
                                                 current_region, extension_size);
        } else if (current_region.get_begin() > 0) {
            current_region = extend_allele_lists(big_allele, small_allele, reference,
                                                 current_region, current_region.get_begin());
        } else {
            break;
        }
        
        // We can continue from previous iterators as list iterators remain valid after list modification
        tie(small_allele_ritr, big_allele_ritr) = mismatch(small_allele_ritr, crend(small_allele), big_allele_ritr);
    } while (small_allele_ritr == crend(small_allele));
    
    // Although the alleles can change, their sizes must remain constant, so we pad to the left
    // of the new alleles to retain the original sizes.
    auto removable_extension = std::distance(small_allele_ritr, crend(small_allele));
    
    if (removable_extension >= small_allele_size) {
        removable_extension -= small_allele_size;
    } else {
        auto required_left_padding = static_cast<SizeType>(small_allele_size - removable_extension);
        // Note this will automatically pad to the right if we've reached the start of the contig
        if (current_region.get_begin() >= required_left_padding) {
            current_region = extend_allele_lists(big_allele, small_allele, reference,
                                                 current_region, required_left_padding);
            removable_extension = 0;
        }
    }
    
    const auto new_big_allele_begin   = next(cbegin(big_allele), removable_extension);
    const auto new_small_allele_begin = next(cbegin(small_allele), removable_extension);
    
    Variant::SequenceType new_big_allele {new_big_allele_begin, next(new_big_allele_begin, big_allele_size)};
    Variant::SequenceType new_small_allele {new_small_allele_begin, next(new_small_allele_begin, small_allele_size)};
    
    const auto new_ref_region_begin = current_region.get_begin() + static_cast<SizeType>(removable_extension);
    const auto new_ref_region_end   = new_ref_region_begin + static_cast<SizeType>(ref_allele_sequence.size());
    
    GenomicRegion new_ref_region {current_region.get_contig_name(), new_ref_region_begin, new_ref_region_end};
    
    if (ref_allele_sequence.size() > alt_allele_sequence.size()) {
        return Variant {move(new_ref_region), move(new_big_allele), move(new_small_allele)};
    } else {
        return Variant {move(new_ref_region), move(new_small_allele), move(new_big_allele)};
    }
}

Variant normalise(const Variant& variant, const ReferenceGenome& reference,
                  const Variant::SizeType extension_size)
{
    return make_parsimonious(left_align(variant, reference, extension_size), reference);
}

Variant pad_left(const Variant& variant, const Variant::SequenceType& sequence)
{
    return Variant {
        expand_lhs(variant, static_cast<GenomicRegion::DifferenceType>(sequence.size())),
        sequence + ref_sequence(variant),
        sequence + alt_sequence(variant)
    };
}

Variant pad_right(const Variant& variant, const Variant::SequenceType& sequence)
{
    return Variant {
        expand_rhs(variant, static_cast<GenomicRegion::DifferenceType>(sequence.size())),
        ref_sequence(variant) + sequence,
        alt_sequence(variant) + sequence
    };
}

Variant pad_left(const Variant& variant, const ReferenceGenome& reference,
                 const Variant::SizeType n)
{
    const auto pad_region = expand_lhs(head_region(variant), static_cast<GenomicRegion::DifferenceType>(n));
    
    const auto pad_sequence = reference.get_sequence(pad_region);
    
    return Variant {
        encompassing_region(pad_region, variant),
        pad_sequence + ref_sequence(variant),
        pad_sequence + alt_sequence(variant)
    };
}

Variant pad_right(const Variant& variant, const ReferenceGenome& reference,
                  const Variant::SizeType n)
{
    const auto pad_region = expand_rhs(tail_region(variant), static_cast<GenomicRegion::DifferenceType>(n));
    
    const auto pad_sequence = reference.get_sequence(pad_region);
    
    return Variant {
        encompassing_region(variant, pad_region),
        ref_sequence(variant) + pad_sequence,
        alt_sequence(variant) + pad_sequence
    };
}

std::vector<Variant> unique_left_align(const std::vector<Variant>& variants,
                                       const ReferenceGenome& reference)
{
    std::vector<Variant> result {};
    result.reserve(variants.size());
    
    boost::transform(variants, std::back_inserter(result),
                     [&reference] (const Variant& variant) {
                         return left_align(variant, reference);
                     });
    
    result.erase(boost::unique<boost::return_found>(boost::sort(result)), std::end(result));
    
    return result;
}

std::vector<Variant> unique_left_align(std::vector<Variant>&& variants,
                                       const ReferenceGenome& reference)
{
    assert(std::is_sorted(std::cbegin(variants), std::cend(variants)));
    assert(std::adjacent_find(std::cbegin(variants), std::cend(variants)) == std::cend(variants));
    
    const auto it = std::stable_partition(std::begin(variants), std::end(variants),
                                          [] (const Variant& variant) {
                                              return !is_left_alignable(variant);
                                          });
    
    std::transform(std::make_move_iterator(it), std::make_move_iterator(end(variants)), it,
                   [&reference] (Variant&& variant) {
                       return left_align(std::move(variant), reference);
                   });
    
    std::sort(it, std::end(variants));
    
    variants.erase(std::unique(it, std::end(variants)), std::end(variants));
    
    std::inplace_merge(std::begin(variants), it, std::end(variants));
    
    return variants;
}

std::vector<Variant> parsimonise_each(const std::vector<Variant>& variants,
                                      const ReferenceGenome& reference)
{
    std::vector<Variant> result {};
    result.reserve(variants.size());
    
    boost::transform(variants, std::back_inserter(result),
                     [&reference] (const auto& variant) {
                         return make_parsimonious(variant, reference);
                     });
    
    return result;
}

std::vector<Variant> parsimonise_each(std::vector<Variant>&& variants,
                                      const ReferenceGenome& reference)
{
    std::vector<Variant> result {};
    result.reserve(variants.size());
    
    std::transform(std::make_move_iterator(std::begin(variants)),
                   std::make_move_iterator(std::end(variants)),
                   std::back_inserter(result),
                   [&reference] (Variant&& variant) {
                       if (!is_parsimonious(variant)) {
                           return make_parsimonious(std::move(variant), reference);
                       }
                       return variant;
                   });
    
    return result;
}

std::vector<Variant> parsimonise_together(const std::vector<Variant>& segment,
                                          const ReferenceGenome& reference)
{
    if (segment.empty()) return {};
    
    auto parsimonised_variants = parsimonise_each(segment, reference);
    
    const auto leftmost  = *leftmost_mappable(parsimonised_variants);
    const auto rightmost = *rightmost_mappable(parsimonised_variants);
    
    std::vector<Variant> result {};
    result.reserve(segment.size());
    
    for (auto&& variant : parsimonised_variants) {
        if (begins_before(leftmost, variant)) {
            variant = pad_left(variant, reference, left_overhang_size(leftmost, variant));
        }
        
        if (ends_before(variant, rightmost)) {
            variant = pad_right(variant, reference, right_overhang_size(rightmost, variant));
        }
        
        result.push_back(std::move(variant));
    }
    
    return result;
}

bool is_snp(const Variant& variant) noexcept
{
    return ref_sequence_size(variant) == 1 && alt_sequence_size(variant) == 1;
}

bool is_insertion(const Variant& variant) noexcept
{
    return ref_sequence_size(variant) < alt_sequence_size(variant);
}

bool is_deletion(const Variant& variant) noexcept
{
    return ref_sequence_size(variant) > alt_sequence_size(variant);
}

bool is_indel(const Variant& variant) noexcept
{
    return is_insertion(variant) || is_deletion(variant);
}

bool is_mnv(const Variant& variant) noexcept
{
    return ref_sequence_size(variant) > 1 && alt_sequence_size(variant) > 1;
}

bool is_transition(const Variant& variant) noexcept
{
    return is_snp(variant)
        && ((ref_sequence(variant) == "A" && alt_sequence(variant) == "G")
         || (ref_sequence(variant) == "G" && alt_sequence(variant) == "A")
         || (ref_sequence(variant) == "C" && alt_sequence(variant) == "T")
         || (ref_sequence(variant) == "T" && alt_sequence(variant) == "C"));
}

bool is_transversion(const Variant& variant) noexcept
{
    return is_snp(variant) && !is_transition(variant);
}

std::vector<Allele::SequenceType> extract_alt_allele_sequences(const std::vector<Variant>& variants)
{
    std::vector<Allele::SequenceType> result {};
    result.reserve(variants.size());
    
    boost::transform(variants, std::back_inserter(result),
                     [] (const auto& variant) { return alt_sequence(variant); });
    
    return result;
}

std::ostream& operator<<(std::ostream& os, const Variant& variant)
{
    os << mapped_region(variant) << " " << ref_sequence(variant) << " " << alt_sequence(variant);
    return os;
}
