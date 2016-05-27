//
//  variant.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant.hpp"

#include <list>
#include <ostream>
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

const GenomicRegion& Variant::mapped_region() const noexcept
{
    return reference_.mapped_region();
}

const Allele& Variant::ref_allele() const noexcept
{
    return reference_;
}

const Allele& Variant::alt_allele() const noexcept
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
    return variant.ref_allele().sequence();
}

const Variant::SequenceType& alt_sequence(const Variant& variant)
{
    return variant.alt_allele().sequence();
}

Variant::SizeType ref_sequence_size(const Variant& variant)
{
    return sequence_size(variant.ref_allele());
}

Variant::SizeType alt_sequence_size(const Variant& variant)
{
    return sequence_size(variant.alt_allele());
}

bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.ref_allele() == rhs.ref_allele() && alt_sequence(lhs) == alt_sequence(rhs);
}

bool operator<(const Variant& lhs, const Variant& rhs)
{
    return (lhs.ref_allele() < rhs.ref_allele() ||
            (lhs.ref_allele() == rhs.ref_allele() && lhs.alt_allele() < rhs.alt_allele()));
}

void remove_duplicates(std::vector<Variant>& variants)
{
    variants.erase(std::unique(std::begin(variants), std::end(variants)), std::end(variants));
}

std::vector<Variant> split_mnv(const Variant& variant)
{
    using std::begin; using std::end; using std::next; using std::prev; using std::distance;
    
    std::vector<Variant> result {};
    result.reserve(ref_sequence_size(variant));
    
    const auto& contig   = contig_name(variant);
    const auto begin_pos = mapped_begin(variant);
    const auto& ref      = ref_sequence(variant);
    const auto& alt      = alt_sequence(variant);
    
    result.emplace_back(contig, begin_pos, ref.front(), alt.front());
    
    const auto last = prev(end(ref));
    
    auto p = std::mismatch(next(begin(ref)), last, next(begin(alt)));
    
    while (p.first != prev(end(ref))) {
        result.emplace_back(contig, begin_pos + distance(begin(ref), p.first), *p.first, *p.second);
        p = std::mismatch(next(p.first), last, next(p.second));
    }
    
    result.emplace_back(contig, mapped_end(variant) - 1, ref.back(), alt.back());
    
    result.shrink_to_fit();
    
    return result;
}

std::vector<Allele> decompose(const std::vector<Variant>& variants)
{
    std::vector<Allele> result {};
    result.reserve(2 * variants.size()); // max num alleles (may be less)
    
    if (variants.empty()) return result;
    
    for (const auto& variant : variants) {
        if (result.empty() || !is_same_region(result.back(), variant)) {
            result.emplace_back(variant.ref_allele());
        }
        result.emplace_back(variant.alt_allele());
    }
    
    result.shrink_to_fit();
    
    return result;
}

std::vector<std::reference_wrapper<const Allele>>
decompose_ref(const std::vector<Variant>& variants)
{
    std::vector<std::reference_wrapper<const Allele>> result {};
    result.reserve(2 * variants.size());
    
    if (variants.empty()) return result;
    
    for (const auto& variant : variants) {
        if (result.empty() || !is_same_region(result.back().get(), variant)) {
            result.emplace_back(variant.ref_allele());
        }
        result.emplace_back(variant.alt_allele());
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
std::size_t get_num_redundant_bases(InputIterator first1, InputIterator last1, InputIterator first2)
{
    return std::distance(first1, std::mismatch(first1, last1, first2).first);
}

template <typename SequenceType>
std::size_t get_num_redundant_front_bases(const SequenceType& smaller_sequence,
                                          const SequenceType& larger_sequence)
{
    return get_num_redundant_bases(std::cbegin(smaller_sequence), std::cend(smaller_sequence),
                                   std::cbegin(larger_sequence));
}

template <typename SequenceType>
std::size_t get_num_redundant_back_bases(const SequenceType& smaller_sequence,
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
    return make_parsimonious(variant, [&reference] (const auto& region) -> char {
        return reference.fetch_sequence(region).front();
    });
}

bool is_left_alignable(const Variant& variant) noexcept
{
    return is_indel(variant);
}

namespace
{
    using AlleleSequence = std::list<Variant::SequenceType::value_type>;

    auto get_alleles(const Variant::SequenceType& allele_a, const Variant::SequenceType& allele_b)
    {
        const auto& alleles = allele_minmax(allele_a, allele_b);
        
        AlleleSequence big_allele {std::cbegin(alleles.second), std::cend(alleles.second)};
        AlleleSequence small_allele {std::cbegin(alleles.first), std::cend(alleles.first)};
        
        return std::make_pair(std::move(big_allele), std::move(small_allele));
    }

    void prepend(const ReferenceGenome::SequenceType& src, AlleleSequence& dst)
    {
        dst.insert(std::begin(dst), std::cbegin(src), std::cend(src));
    }
    
    GenomicRegion extend_alleles(AlleleSequence& big_allele, AlleleSequence& small_allele,
                                 const ReferenceGenome& reference, const GenomicRegion& current_region,
                                 const Variant::SizeType extension_size)
    {
        const auto new_region = shift(current_region, -extension_size);
        
        const auto extension_region = expand_rhs(head_region(new_region), extension_size);
        
        const auto extension_sequence = reference.fetch_sequence(extension_region);
        
        prepend(extension_sequence, big_allele);
        prepend(extension_sequence, small_allele);
        
        return new_region;
    }
} // namespace

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
    
    AlleleSequence big_allele {}, small_allele {};
    
    tie(big_allele, small_allele) = get_alleles(ref_allele_sequence, alt_allele_sequence);
    
    auto big_allele_ritr   = crbegin(big_allele);
    auto small_allele_ritr = crbegin(small_allele);
    
    const auto big_allele_size   = static_cast<SizeType>(big_allele.size());
    const auto small_allele_size = static_cast<SizeType>(small_allele.size());
    
    GenomicRegion current_region {variant.mapped_region()};
    
    do {
        if (current_region.begin() >= extension_size) {
            current_region = extend_alleles(big_allele, small_allele, reference,
                                            current_region, extension_size);
        } else if (current_region.begin() > 0) {
            current_region = extend_alleles(big_allele, small_allele, reference,
                                            current_region, current_region.begin());
        } else {
            break;
        }
        
        // We can continue from previous iterators as list iterators remain valid after modification
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
        if (current_region.begin() >= required_left_padding) {
            current_region = extend_alleles(big_allele, small_allele, reference,
                                            current_region, required_left_padding);
            removable_extension = 0;
        }
    }
    
    const auto new_big_allele_begin   = next(cbegin(big_allele), removable_extension);
    const auto new_small_allele_begin = next(cbegin(small_allele), removable_extension);
    
    Variant::SequenceType new_big_allele {new_big_allele_begin, next(new_big_allele_begin, big_allele_size)};
    Variant::SequenceType new_small_allele {new_small_allele_begin, next(new_small_allele_begin, small_allele_size)};
    
    const auto new_ref_region_begin = current_region.begin() + static_cast<SizeType>(removable_extension);
    const auto new_ref_region_end   = new_ref_region_begin + static_cast<SizeType>(ref_allele_sequence.size());
    
    GenomicRegion new_ref_region {current_region.contig_name(), new_ref_region_begin, new_ref_region_end};
    
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
    
    const auto pad_sequence = reference.fetch_sequence(pad_region);
    
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
    
    const auto pad_sequence = reference.fetch_sequence(pad_region);
    
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
    return ref_sequence_size(variant) == alt_sequence_size(variant) && ref_sequence_size(variant) > 1;
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
    os << mapped_region(variant) << " " << ref_sequence(variant) << " -> " << alt_sequence(variant);
    return os;
}
