//
//  variant_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_utils.h"

#include <iterator> // std::tie, std::next, std::distance
#include <list>

#include "reference_genome.h"
#include "genomic_region.h"
#include "variant_factory.h"
#include "variant_candidate_generator.h"

using std::cbegin;
using std::cend;
using std::crbegin;
using std::crend;

auto allele_minmax(const Variant::SequenceType& allele_a, const Variant::SequenceType& allele_b)
{
    static auto is_bigger = [] (const auto& a1, const auto& a2) { return a1.size() < a2.size(); };
    return std::minmax(allele_a, allele_b, is_bigger);
}

bool is_parsimonious(const Variant& a_variant) noexcept
{
    if (reference_allele_size(a_variant) == 0 || alternative_allele_size(a_variant) == 0) {
        return false;
    }
    
    const auto& ref_allele = a_variant.get_reference_allele();
    const auto& alt_allele = a_variant.get_alternative_allele();
    
    const auto& alleles = allele_minmax(ref_allele, alt_allele);
    const auto& the_small_allele = alleles.first;
    const auto& the_big_allele   = alleles.second;
    
    if (num_redundant_bases(cbegin(the_small_allele), cend(the_small_allele),
                            cbegin(the_big_allele)) > 1) {
        return false;
    }
    
    if (the_small_allele.size() > 1 &&
        num_redundant_bases(crbegin(the_small_allele), crend(the_small_allele),
                            crbegin(the_big_allele)) > 0) {
        return false;
    }
    
    return true;
}

Variant make_parsimonious(const Variant& a_variant, ReferenceGenome& the_reference,
                          VariantFactory& a_variant_factory)
{
    if (is_parsimonious(a_variant)) return a_variant;
    
    const auto& old_ref_allele = a_variant.get_reference_allele();
    const auto& old_alt_allele = a_variant.get_alternative_allele();
    
    const auto& alleles = allele_minmax(old_ref_allele, old_alt_allele);
    const auto& the_small_allele = alleles.first;
    const auto& the_big_allele   = alleles.second;
    
    const auto& old_ref_region = a_variant.get_reference_allele_region();
    
    if (the_small_allele.size() == 0) {
        if (old_ref_region.get_begin() > 0) {
            GenomicRegion new_ref_region {old_ref_region.get_contig_name(),
                old_ref_region.get_begin() - 1, old_ref_region.get_end()};
            
            auto new_ref_allele = the_reference.get_sequence(new_ref_region);
            auto new_alt_allele = new_ref_allele.front() + a_variant.get_alternative_allele();
            
            return a_variant_factory.make(std::move(new_ref_region), std::move(new_ref_allele),
                                          std::move(new_alt_allele));
        } else {
            // In this very rare care the only option is to pad to the right.
            GenomicRegion new_ref_region {old_ref_region.get_contig_name(),
                old_ref_region.get_begin(), old_ref_region.get_end() + 1};
            
            auto new_ref_allele = the_reference.get_sequence(new_ref_region);
            auto new_alt_allele = a_variant.get_alternative_allele() + new_ref_allele.back();
            
            return a_variant_factory.make(std::move(new_ref_region), std::move(new_ref_allele),
                                          std::move(new_alt_allele));
        }
    }
    
    auto num_redundant_back_bases  = num_redundant_bases(crbegin(the_small_allele),
                                                         crend(the_small_allele), crbegin(the_big_allele));
    auto num_redundant_front_bases = num_redundant_bases(cbegin(the_small_allele),
                                                         cend(the_small_allele), cbegin(the_big_allele));
    
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
        a_variant_factory.make(std::move(new_ref_region), std::move(new_big_allele),
                               std::move(new_small_allele)) :
        a_variant_factory.make(std::move(new_ref_region), std::move(new_small_allele),
                               std::move(new_big_allele));
}

bool is_left_alignable(const Variant& a_variant) noexcept
{
    return is_indel(a_variant);
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
                                  ReferenceGenome& the_reference, const GenomicRegion& current_region,
                                  Variant::SizeType extension_size)
{
    auto new_region = shift(current_region, -extension_size);
    GenomicRegion extension_region {new_region.get_contig_name(), new_region.get_begin(),
        new_region.get_begin() + extension_size};
    auto the_extension = the_reference.get_sequence(extension_region);
    
    big_allele.insert(begin(big_allele), cbegin(the_extension), cend(the_extension));
    small_allele.insert(begin(small_allele), cbegin(the_extension), cend(the_extension));
    
    return new_region;
}

Variant left_align(const Variant& a_variant, ReferenceGenome& the_reference,
                   VariantFactory& a_variant_factory, Variant::SizeType extension_size)
{
    if (!is_left_alignable(a_variant)) {
        return a_variant;
    }
    
    static_assert(sizeof(Variant::SizeType) == sizeof(GenomicRegion::SizeType),
                  "Variant and GenomicRegion have different SizeType");
    using SizeType = Variant::SizeType;
    
    const auto& ref_allele = a_variant.get_reference_allele();
    const auto& alt_allele = a_variant.get_alternative_allele();
    
    LeftAlignmentList big_allele {}, small_allele {};
    std::tie(big_allele, small_allele) = get_allele_lists(ref_allele, alt_allele);
    auto big_allele_ritr   = crbegin(big_allele);
    auto small_allele_ritr = crbegin(small_allele);
    
    auto big_allele_size   = static_cast<SizeType>(big_allele.size());
    auto small_allele_size = static_cast<SizeType>(small_allele.size());
    
    GenomicRegion current_region {a_variant.get_reference_allele_region()};
    
    do {
        if (current_region.get_begin() >= extension_size) {
            current_region = extend_allele_lists(big_allele, small_allele, the_reference,
                                                 current_region, extension_size);
        } else if (current_region.get_begin() > 0) {
            current_region = extend_allele_lists(big_allele, small_allele, the_reference,
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
            current_region = extend_allele_lists(big_allele, small_allele, the_reference,
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
        a_variant_factory.make(std::move(new_ref_region), std::move(new_big_allele),
                               std::move(new_small_allele)) :
        a_variant_factory.make(std::move(new_ref_region), std::move(new_small_allele),
                               std::move(new_big_allele));
}

Variant normalise(const Variant& a_variant, ReferenceGenome& the_reference,
                  VariantFactory& a_variant_factory, Variant::SizeType extension_size)
{
    return make_parsimonious(left_align(a_variant, the_reference, a_variant_factory, extension_size),
                             the_reference, a_variant_factory);
}

bool is_snp(const Variant& a_variant) noexcept
{
    return reference_allele_size(a_variant) == 1 && alternative_allele_size(a_variant) == 1;
}

bool is_insertion(const Variant& a_variant) noexcept
{
    return reference_allele_size(a_variant) < alternative_allele_size(a_variant);
}

bool is_deletion(const Variant& a_variant) noexcept
{
    return reference_allele_size(a_variant) > alternative_allele_size(a_variant);
}

bool is_indel(const Variant& a_variant) noexcept
{
    return is_insertion(a_variant) || is_deletion(a_variant);
}

bool is_mnv(const Variant& a_variant) noexcept
{
    return reference_allele_size(a_variant) > 1 && alternative_allele_size(a_variant) > 1;
}
