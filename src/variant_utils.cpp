//
//  variant_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_utils.h"

#include <algorithm> // std::mismatch, std::max
#include <iterator> // std::tie, std::next, std::distance
#include <list>

#include "reference_genome.h"
#include "genomic_region.h"
#include "variant_factory.h"

using std::cbegin;
using std::cend;
using std::crbegin;
using std::crend;

bool is_parsimonious(const Variant& a_variant) noexcept
{
    if (reference_allele_size(a_variant) == 0 || alternative_allele_size(a_variant) == 0) {
        return false;
    }
    
    const auto& ref_allele = a_variant.get_reference_allele();
    const auto& alt_allele = a_variant.get_alternative_allele();
    
    if (std::mismatch(cbegin(ref_allele), cend(ref_allele),
                      cbegin(alt_allele)).first != cbegin(ref_allele)) {
        return false;
    }
    
    if (std::mismatch(crbegin(ref_allele), crend(ref_allele),
                      crbegin(alt_allele)).first != crbegin(ref_allele)) {
        return false;
    }
    
    return true;
}

Variant make_parsimonious(const Variant& a_variant)
{
    return a_variant;
}

bool is_left_alignable(const Variant& a_variant) noexcept
{
    return is_indel(a_variant);
}

using LeftAlignmentList = std::list<Variant::StringType::value_type>;

auto get_allele_lists(const Variant::StringType& allele_a, const Variant::StringType& allele_b)
{
    bool is_a_bigger = allele_a.size() > allele_b.size();
    const auto& the_big_allele   = (is_a_bigger) ? allele_a : allele_b;
    const auto& the_small_allele = (is_a_bigger) ? allele_b : allele_a;
    
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
                   Variant::SizeType extension_size)
{
    if (!is_left_alignable(a_variant)) {
        return a_variant;
    }
    
    const auto& ref_allele = a_variant.get_reference_allele();
    const auto& alt_allele = a_variant.get_alternative_allele();
    
    LeftAlignmentList big_allele {}, small_allele {};
    std::tie(big_allele, small_allele) = get_allele_lists(ref_allele, alt_allele);
    auto big_allele_ritr = crbegin(big_allele);
    auto small_allele_ritr = crbegin(small_allele);
    
    auto big_allele_size   = static_cast<Variant::SizeType>(big_allele.size());
    auto small_allele_size = static_cast<Variant::SizeType>(small_allele.size());
    
    GenomicRegion current_region {a_variant.get_reference_allele_region()};
    
    do {
        current_region = extend_allele_lists(big_allele, small_allele, the_reference,
                                             current_region, extension_size);
        
        std::tie(small_allele_ritr, big_allele_ritr)
                = std::mismatch(small_allele_ritr, crend(small_allele), big_allele_ritr);
    } while (small_allele_ritr == crend(small_allele));
    
    auto removable_extension = std::distance(small_allele_ritr, crend(small_allele));
    auto new_big_allele_begin   = std::next(cbegin(big_allele), removable_extension);
    auto new_small_allele_begin = std::next(cbegin(small_allele), removable_extension);
    
    Variant::StringType new_big_allele {new_big_allele_begin,
                                        std::next(new_big_allele_begin, big_allele_size)};
    Variant::StringType new_small_allele {new_small_allele_begin,
                                        std::next(new_small_allele_begin, small_allele_size)};
    
    auto new_ref_region_begin = current_region.get_begin() +
                                static_cast<GenomicRegion::SizeType>(removable_extension);
    auto new_ref_region_end   = new_ref_region_begin +
                                static_cast<GenomicRegion::SizeType>(ref_allele.size());
    GenomicRegion new_ref_region {current_region.get_contig_name(), new_ref_region_begin,
                                    new_ref_region_end};
    
    VariantFactory a_factory {};
    return (ref_allele.size() > alt_allele.size()) ?
        a_factory.make(new_ref_region, new_big_allele, new_small_allele) :
        a_factory.make(new_ref_region, new_small_allele, new_big_allele);
}

bool is_normalised(const Variant& a_variant) noexcept
{
    return is_parsimonious(a_variant) && is_normalised(a_variant);
}

Variant normalise(const Variant& a_variant, ReferenceGenome& the_reference)
{
    return left_align(make_parsimonious(a_variant), the_reference);
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
