// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant.hpp"

#include <list>
#include <ostream>
#include <cassert>

#include <boost/range/algorithm.hpp>

#include "io/reference/reference_genome.hpp"
#include "utils/mappable_algorithms.hpp"

namespace octopus {

// public methods

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

const Variant::NucleotideSequence& ref_sequence(const Variant& variant)
{
    return variant.ref_allele().sequence();
}

const Variant::NucleotideSequence& alt_sequence(const Variant& variant)
{
    return variant.alt_allele().sequence();
}

Variant::NucleotideSequence::size_type  ref_sequence_size(const Variant& variant)
{
    return sequence_size(variant.ref_allele());
}

Variant::NucleotideSequence::size_type  alt_sequence_size(const Variant& variant)
{
    return sequence_size(variant.alt_allele());
}

bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.ref_allele() == rhs.ref_allele() && alt_sequence(lhs) == alt_sequence(rhs);
}

bool operator<(const Variant& lhs, const Variant& rhs)
{
    return (lhs.ref_allele() < rhs.ref_allele()
            || (lhs.ref_allele() == rhs.ref_allele() && lhs.alt_allele() < rhs.alt_allele()));
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
    for (const auto& variant : variants) {
        if (result.empty() || !is_same_region(result.back(), variant)) {
            result.emplace_back(variant.ref_allele());
        }
        result.emplace_back(variant.alt_allele());
    }
    result.shrink_to_fit();
    return result;
}

std::vector<std::reference_wrapper<const Allele>> decompose_ref(const std::vector<Variant>& variants)
{
    std::vector<std::reference_wrapper<const Allele>> result {};
    result.reserve(2 * variants.size());
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

auto allele_minmax(const Variant::NucleotideSequence& lhs, const Variant::NucleotideSequence& rhs)
{
    return std::minmax(lhs, rhs, [] (const auto& lhs, const auto& rhs) { return lhs.size() < rhs.size(); });
}

template <typename ForwardIt>
auto count_redundant_bases(ForwardIt first1, ForwardIt last1, ForwardIt first2)
{
    return std::distance(first1, std::mismatch(first1, last1, first2).first);
}

template <typename Sequence>
auto count_redundant_front_bases(const Sequence& smaller_sequence, const Sequence& larger_sequence)
{
    return count_redundant_bases(std::cbegin(smaller_sequence), std::cend(smaller_sequence),
                                 std::cbegin(larger_sequence));
}

template <typename Sequence>
auto count_redundant_back_bases(const Sequence& smaller_sequence, const Sequence& larger_sequence)
{
    return count_redundant_bases(std::crbegin(smaller_sequence), std::crend(smaller_sequence),
                                 std::crbegin(larger_sequence));
}

bool can_trim(const Variant& v)
{
    if (ref_sequence_size(v) == 0 || alt_sequence_size(v) == 0) return false;
    
    const auto& p = allele_minmax(ref_sequence(v), alt_sequence(v));
    
    return count_redundant_front_bases(p.first, p.second) > 0
            || count_redundant_back_bases(p.first, p.second) > 0;
}

Variant trim(const Variant& v)
{
    using std::cbegin; using std::cend; using std::crbegin; using std::make_reverse_iterator;
    
    const auto& ref = ref_sequence(v);
    const auto& alt = alt_sequence(v);
    const auto m1 = std::mismatch(cbegin(ref), cend(ref), cbegin(alt), cend(alt));
    const auto m2 = std::mismatch(crbegin(ref), make_reverse_iterator(m1.first),
                                  crbegin(alt), make_reverse_iterator(m1.second));
    using NucleotideSequence = Variant::NucleotideSequence;
    NucleotideSequence new_ref {m1.first, m2.first.base()};
    NucleotideSequence new_alt {m1.second, m2.second.base()};
    const auto pad_front = std::distance(m1.first, cbegin(ref));
    const auto pad_back  = std::distance(m2.first, crbegin(ref));
    auto new_region = expand(mapped_region(v), pad_front, pad_back);
    return Variant {std::move(new_region), std::move(new_ref), std::move(new_alt)};
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
    if (count_redundant_front_bases(small_allele, big_allele) > 1) {
        return false;
    }
    if (small_allele.size() > 1 && count_redundant_back_bases(small_allele, big_allele) > 0) {
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

namespace {

using NucleotideList = std::list<char>;

auto get_alleles(const Variant::NucleotideSequence& allele_a, const Variant::NucleotideSequence& allele_b)
{
    const auto& alleles = allele_minmax(allele_a, allele_b);
    NucleotideList big_allele {std::cbegin(alleles.second), std::cend(alleles.second)};
    NucleotideList small_allele {std::cbegin(alleles.first), std::cend(alleles.first)};
    return std::make_pair(std::move(big_allele), std::move(small_allele));
}

void prepend(const ReferenceGenome::GeneticSequence& src, NucleotideList& dst)
{
    dst.insert(std::begin(dst), std::cbegin(src), std::cend(src));
}

GenomicRegion extend_alleles(NucleotideList& big_allele, NucleotideList& small_allele,
                             const ReferenceGenome& reference, const GenomicRegion& current_region,
                             const GenomicRegion::Distance extension_size)
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
                   const GenomicRegion::Size extension_size)
{
    using std::cbegin; using std::cend; using std::crbegin; using std::crend;
    using std::tie; using std::next; using std::mismatch;
    
    if (!is_left_alignable(variant)) return variant;
    
    using SizeType = Variant::MappingDomain::Size;
    const auto& ref_allele_sequence = ref_sequence(variant);
    const auto& alt_allele_sequence = alt_sequence(variant);
    NucleotideList big_allele {}, small_allele {};
    tie(big_allele, small_allele) = get_alleles(ref_allele_sequence, alt_allele_sequence);
    const auto big_allele_size   = static_cast<SizeType>(big_allele.size());
    const auto small_allele_size = static_cast<SizeType>(small_allele.size());
    GenomicRegion current_region {variant.mapped_region()};
    auto big_allele_ritr   = crbegin(big_allele);
    auto small_allele_ritr = crbegin(small_allele);
    
    do {
        if (current_region.begin() > 0) {
            current_region = extend_alleles(big_allele, small_allele, reference, current_region,
                                            std::min(extension_size, current_region.begin()));
        } else {
            break;
        }
        // We can continue from previous iterators as list iterators remain valid after modification
        tie(small_allele_ritr, big_allele_ritr) = mismatch(small_allele_ritr, crend(small_allele), big_allele_ritr);
    } while (small_allele_ritr == crend(small_allele));
    
    // Although the alleles can change, their sizes must remain constant, so we pad to the left
    // of the new alleles to retain the original sizes.
    auto removable_extension = static_cast<SizeType>(std::distance(small_allele_ritr, crend(small_allele)));
    if (removable_extension >= small_allele_size) {
        removable_extension -= small_allele_size;
    } else {
        const auto required_left_padding = small_allele_size - removable_extension;
        // Note this will automatically pad to the right if we've reached the start of the contig
        if (current_region.begin() >= required_left_padding) {
            current_region = extend_alleles(big_allele, small_allele, reference, current_region, required_left_padding);
            removable_extension = 0;
        }
    }
    
    const auto new_big_allele_begin   = next(cbegin(big_allele), removable_extension);
    const auto new_small_allele_begin = next(cbegin(small_allele), removable_extension);
    Variant::NucleotideSequence new_big_allele {new_big_allele_begin, next(new_big_allele_begin, big_allele_size)};
    Variant::NucleotideSequence new_small_allele {new_small_allele_begin, next(new_small_allele_begin, small_allele_size)};
    const auto new_ref_region_begin = current_region.begin() + removable_extension;
    const auto new_ref_region_end   = new_ref_region_begin + static_cast<SizeType>(ref_allele_sequence.size());
    GenomicRegion new_ref_region {current_region.contig_name(), new_ref_region_begin, new_ref_region_end};
    
    if (ref_allele_sequence.size() > alt_allele_sequence.size()) {
        return Variant {std::move(new_ref_region), std::move(new_big_allele), std::move(new_small_allele)};
    } else {
        return Variant {std::move(new_ref_region), std::move(new_small_allele), std::move(new_big_allele)};
    }
}

Variant normalise(const Variant& variant, const ReferenceGenome& reference,
                  const unsigned extension_size)
{
    return make_parsimonious(left_align(variant, reference, extension_size), reference);
}

Variant pad_left(const Variant& variant, const Variant::NucleotideSequence& sequence)
{
    return Variant {
        expand_lhs(mapped_region(variant), static_cast<GenomicRegion::Distance>(sequence.size())),
        sequence + ref_sequence(variant),
        sequence + alt_sequence(variant)
    };
}

Variant pad_right(const Variant& variant, const Variant::NucleotideSequence& sequence)
{
    return Variant {
        expand_rhs(mapped_region(variant), static_cast<GenomicRegion::Distance>(sequence.size())),
        ref_sequence(variant) + sequence,
        alt_sequence(variant) + sequence
    };
}

Variant pad_left(const Variant& variant, const ReferenceGenome& reference, const unsigned n)
{
    const auto pad_region   = expand_lhs(head_region(variant), static_cast<GenomicRegion::Distance>(n));
    const auto pad_sequence = reference.fetch_sequence(pad_region);
    return Variant {
        encompassing_region(pad_region, variant),
        pad_sequence + ref_sequence(variant),
        pad_sequence + alt_sequence(variant)
    };
}

Variant pad_right(const Variant& variant, const ReferenceGenome& reference, const unsigned n)
{
    const auto pad_region = expand_rhs(tail_region(variant), static_cast<GenomicRegion::Distance>(n));
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
                     [&reference] (const Variant& variant) { return left_align(variant, reference); });
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
    return std::move(variants);
}

std::vector<Variant> parsimonise_each(const std::vector<Variant>& variants,
                                      const ReferenceGenome& reference)
{
    std::vector<Variant> result {};
    result.reserve(variants.size());
    boost::transform(variants, std::back_inserter(result),
                     [&reference] (const auto& variant) { return make_parsimonious(variant, reference); });
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
                       return std::move(variant);
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

bool is_snv(const Variant& variant) noexcept
{
    return ref_sequence_size(variant) == 1 && alt_sequence_size(variant) == 1;
}

bool is_mnv(const Variant& variant) noexcept
{
    return ref_sequence_size(variant) == alt_sequence_size(variant) && ref_sequence_size(variant) > 1;
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

bool is_simple_insertion(const Variant& variant) noexcept
{
    return is_insertion(variant) && is_empty_region(variant);
}

bool is_simple_deletion(const Variant& variant) noexcept
{
    return is_deletion(variant) && is_sequence_empty(variant.alt_allele());
}

bool is_simple_indel(const Variant& variant) noexcept
{
    return is_simple_insertion(variant) || is_simple_deletion(variant);
}

bool are_same_type(const Variant& lhs, const Variant& rhs) noexcept
{
    if (ref_sequence_size(lhs) == alt_sequence_size(rhs)) {
        return (ref_sequence_size(lhs) == 1) ? is_snv(rhs) : is_mnv(rhs);
    }
    if (is_insertion(lhs)) return is_insertion(rhs);
    return is_deletion(rhs);
}

bool is_transition(const Variant& variant) noexcept
{
    return is_snv(variant)
        && ((ref_sequence(variant) == "A" && alt_sequence(variant) == "G")
         || (ref_sequence(variant) == "G" && alt_sequence(variant) == "A")
         || (ref_sequence(variant) == "C" && alt_sequence(variant) == "T")
         || (ref_sequence(variant) == "T" && alt_sequence(variant) == "C"));
}

bool is_transversion(const Variant& variant) noexcept
{
    return is_snv(variant) && !is_transition(variant);
}

Variant::NucleotideSequence::size_type indel_size(const Variant& variant) noexcept
{
    const auto p = std::minmax({ref_sequence_size(variant), alt_sequence_size(variant)});
    return p.second - p.first;
}

std::vector<Allele::NucleotideSequence> extract_alt_allele_sequences(const std::vector<Variant>& variants)
{
    std::vector<Allele::NucleotideSequence> result {};
    result.reserve(variants.size());
    boost::transform(variants, std::back_inserter(result), [] (const auto& variant) { return alt_sequence(variant); });
    return result;
}

std::ostream& operator<<(std::ostream& os, const Variant& variant)
{
    os << mapped_region(variant) << " " << ref_sequence(variant) << " -> " << alt_sequence(variant);
    return os;
}

} // namespace octopus
