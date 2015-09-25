//
//  variant_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_utils.hpp"

#include <list>
#include <utility> // std::move

#include "reference_genome.hpp"
#include "candidate_variant_generator.hpp"
#include "mappable_algorithms.hpp"

using std::cbegin;
using std::cend;
using std::crbegin;
using std::crend;

template <typename InputIterator>
std::size_t num_redundant_bases(InputIterator first1, InputIterator last1, InputIterator first2)
{
    return std::distance(first1, std::mismatch(first1, last1, first2).first);
}

void remove_duplicates(std::vector<Variant>& variants)
{
    variants.erase(std::unique(std::begin(variants), std::end(variants)), variants.end());
}

std::vector<Allele> decompose(const std::vector<Variant>& variants)
{
    std::vector<Allele> result {};
    result.reserve(2 * variants.size()); // max num alleles (may be less)
    
    if (variants.empty()) return result;
    
    for (const auto& variant : variants) {
        if (result.empty() || get_region(result.back()) != get_region(variant)) {
            result.emplace_back(variant.get_reference_allele());
        }
        result.emplace_back(variant.get_alternative_allele());
    }
    
    result.shrink_to_fit();
    
    return result;
}

std::vector<Allele> get_intervening_reference_alleles(const std::vector<Variant>& the_variants,
                                                      ReferenceGenome& reference)
{
    auto all_overlapped_regions = get_covered_regions(std::cbegin(the_variants), std::cend(the_variants));
    auto regions_between_candidates = get_all_intervening(std::cbegin(all_overlapped_regions),
                                                          std::cend(all_overlapped_regions));
    
    std::vector<Allele> result {};
    result.reserve(regions_between_candidates.size());
    
    std::transform(std::cbegin(regions_between_candidates), std::cend(regions_between_candidates),
                   std::back_inserter(result), [&reference] (const auto& region) {
                       return get_reference_allele(region, reference);
                   });
    
    return result;
}

auto allele_minmax(const Variant::SequenceType& lhs, const Variant::SequenceType& rhs)
{
    return std::minmax(lhs, rhs, [] (const auto& a1, const auto& a2) { return a1.size() < a2.size(); });
}

bool is_parsimonious(const Variant& variant) noexcept
{
    if (variant.reference_allele_size() == 0 || variant.alternative_allele_size() == 0) {
        return false;
    }
    
    const auto& ref_allele = variant.get_reference_allele_sequence();
    const auto& alt_allele = variant.get_alternative_allele_sequence();
    
    const auto& alleles = allele_minmax(ref_allele, alt_allele);
    const auto& the_small_allele = alleles.first;
    const auto& the_big_allele   = alleles.second;
    
    if (num_redundant_bases(cbegin(the_small_allele), cend(the_small_allele), cbegin(the_big_allele)) > 1) {
        return false;
    }
    
    if (the_small_allele.size() > 1 &&
        num_redundant_bases(crbegin(the_small_allele), crend(the_small_allele), crbegin(the_big_allele)) > 0) {
        return false;
    }
    
    return true;
}

Variant make_parsimonious(const Variant& variant, ReferenceGenome& reference)
{
    if (is_parsimonious(variant)) return variant;
    
    const auto& old_ref_allele = variant.get_reference_allele_sequence();
    const auto& old_alt_allele = variant.get_alternative_allele_sequence();
    
    const auto& alleles = allele_minmax(old_ref_allele, old_alt_allele);
    const auto& the_small_allele = alleles.first;
    const auto& the_big_allele   = alleles.second;
    
    const auto& old_ref_region = variant.get_region();
    
    if (the_small_allele.size() == 0) {
        if (old_ref_region.get_begin() > 0) {
            GenomicRegion new_ref_region {old_ref_region.get_contig_name(),
                old_ref_region.get_begin() - 1, old_ref_region.get_end()};
            
            auto new_ref_allele = reference.get_sequence(new_ref_region);
            auto new_alt_allele = new_ref_allele.front() + variant.get_alternative_allele_sequence();
            
            return Variant {std::move(new_ref_region), std::move(new_ref_allele),
                            std::move(new_alt_allele)};
        } else {
            // In this very rare care the only option is to pad to the right.
            GenomicRegion new_ref_region {old_ref_region.get_contig_name(),
                old_ref_region.get_begin(), old_ref_region.get_end() + 1};
            
            auto new_ref_allele = reference.get_sequence(new_ref_region);
            auto new_alt_allele = variant.get_alternative_allele_sequence() + new_ref_allele.back();
            
            return Variant {std::move(new_ref_region), std::move(new_ref_allele),
                            std::move(new_alt_allele)};
        }
    }
    
    auto num_redundant_back_bases  = num_redundant_bases(crbegin(the_small_allele), crend(the_small_allele),
                                                         crbegin(the_big_allele));
    auto num_redundant_front_bases = num_redundant_bases(cbegin(the_small_allele), cend(the_small_allele),
                                                         cbegin(the_big_allele));
    
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
    
    return (old_ref_allele.size() > old_alt_allele.size()) ?
            Variant {std::move(new_ref_region), std::move(new_big_allele),
                    std::move(new_small_allele)}
            :
            Variant {std::move(new_ref_region), std::move(new_small_allele),
                    std::move(new_big_allele)};
}

bool is_left_alignable(const Variant& variant) noexcept
{
    return is_indel(variant);
}

using LeftAlignmentList = std::list<Variant::SequenceType::value_type>;

auto get_allele_lists(const Variant::SequenceType& allele_a, const Variant::SequenceType& allele_b)
{
    const auto& alleles = allele_minmax(allele_a, allele_b);
    const auto& the_small_allele = alleles.first;
    const auto& the_big_allele   = alleles.second;
    
    LeftAlignmentList big_allele {cbegin(the_big_allele), cend(the_big_allele)};
    LeftAlignmentList small_allele {cbegin(the_small_allele), cend(the_small_allele)};
    
    return std::make_pair(big_allele, small_allele);
}

GenomicRegion extend_allele_lists(LeftAlignmentList& big_allele, LeftAlignmentList& small_allele,
                                  ReferenceGenome& reference, const GenomicRegion& current_region,
                                  Variant::SizeType extension_size)
{
    auto new_region = shift(current_region, -extension_size);
    GenomicRegion extension_region {new_region.get_contig_name(), new_region.get_begin(),
        new_region.get_begin() + extension_size};
    auto the_extension = reference.get_sequence(extension_region);
    
    big_allele.insert(begin(big_allele), cbegin(the_extension), cend(the_extension));
    small_allele.insert(begin(small_allele), cbegin(the_extension), cend(the_extension));
    
    return new_region;
}

Variant left_align(const Variant& variant, ReferenceGenome& reference,
                   Variant::SizeType extension_size)
{
    if (!is_left_alignable(variant)) {
        return variant;
    }
    
    static_assert(sizeof(Variant::SizeType) == sizeof(GenomicRegion::SizeType),
                  "Variant and GenomicRegion have different SizeType");
    using SizeType = Variant::SizeType;
    
    const auto& ref_allele = variant.get_reference_allele_sequence();
    const auto& alt_allele = variant.get_alternative_allele_sequence();
    
    LeftAlignmentList big_allele {}, small_allele {};
    std::tie(big_allele, small_allele) = get_allele_lists(ref_allele, alt_allele);
    auto big_allele_ritr   = crbegin(big_allele);
    auto small_allele_ritr = crbegin(small_allele);
    
    auto big_allele_size   = static_cast<SizeType>(big_allele.size());
    auto small_allele_size = static_cast<SizeType>(small_allele.size());
    
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
        
        // We can continue from previous iterators as list iterators aren't invalidated
        // by modifications to the list.
        std::tie(small_allele_ritr, big_allele_ritr)
                = std::mismatch(small_allele_ritr, crend(small_allele), big_allele_ritr);
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
    
    auto new_big_allele_begin   = std::next(cbegin(big_allele), removable_extension);
    auto new_small_allele_begin = std::next(cbegin(small_allele), removable_extension);
    
    Variant::SequenceType new_big_allele {new_big_allele_begin,
                                            std::next(new_big_allele_begin, big_allele_size)};
    Variant::SequenceType new_small_allele {new_small_allele_begin,
                                            std::next(new_small_allele_begin, small_allele_size)};
    
    auto new_ref_region_begin = current_region.get_begin() + static_cast<SizeType>(removable_extension);
    auto new_ref_region_end   = new_ref_region_begin + static_cast<SizeType>(ref_allele.size());
    GenomicRegion new_ref_region {current_region.get_contig_name(), new_ref_region_begin,
                                    new_ref_region_end};
    
    return (ref_allele.size() > alt_allele.size()) ?
            Variant {std::move(new_ref_region), std::move(new_big_allele),
                    std::move(new_small_allele)}
            :
            Variant {std::move(new_ref_region), std::move(new_small_allele),
                    std::move(new_big_allele)};
}

Variant normalise(const Variant& variant, ReferenceGenome& reference, Variant::SizeType extension_size)
{
    return make_parsimonious(left_align(variant, reference, extension_size), reference);
}

std::vector<Variant> unique_left_align(const std::vector<Variant>& variants, ReferenceGenome& reference)
{
    std::vector<Variant> result {};
    result.reserve(variants.size());
    
    std::transform(variants.cbegin(), variants.cend(), std::back_inserter(result),
                   [&reference] (const Variant& variant) {
                       return left_align(variant, reference);
                   });
    
    if (!std::is_sorted(result.cbegin(), result.cend())) {
        std::sort(result.begin(), result.end());
    }
    
    result.erase(std::unique(result.begin(), result.end()), result.end());
    
    return result;
}

std::vector<Variant> make_parsimonious(const std::vector<Variant>& variants, ReferenceGenome& reference)
{
    std::vector<Variant> result {};
    result.reserve(variants.size());
    std::transform(std::cbegin(variants), std::cend(variants), std::back_inserter(result),
                   [&reference] (const auto& variant) { return make_parsimonious(variant, reference); });
    return result;
}

bool is_snp(const Variant& variant) noexcept
{
    return variant.reference_allele_size() == 1 && variant.alternative_allele_size() == 1;
}

bool is_insertion(const Variant& variant) noexcept
{
    return variant.reference_allele_size() < variant.alternative_allele_size();
}

bool is_deletion(const Variant& variant) noexcept
{
    return variant.reference_allele_size() > variant.alternative_allele_size();
}

bool is_indel(const Variant& variant) noexcept
{
    return is_insertion(variant) || is_deletion(variant);
}

bool is_mnv(const Variant& variant) noexcept
{
    return variant.reference_allele_size() > 1 && variant.alternative_allele_size() > 1;
}

bool is_transition(const Variant& variant) noexcept
{
    return is_snp(variant) && (
           (variant.get_reference_allele_sequence() == "A" && variant.get_alternative_allele_sequence() == "G")
        || (variant.get_reference_allele_sequence() == "G" && variant.get_alternative_allele_sequence() == "A")
        || (variant.get_reference_allele_sequence() == "C" && variant.get_alternative_allele_sequence() == "T")
        || (variant.get_reference_allele_sequence() == "T" && variant.get_alternative_allele_sequence() == "C")
    );
}

bool is_transversion(const Variant& variant) noexcept
{
    return is_snp(variant) && !is_transition(variant);
}
