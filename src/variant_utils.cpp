//
//  variant_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_utils.h"

#include <algorithm> // std::mismatch, std::max
#include <iterator> // std::tie
#include <type_traits> // std::remove_reference_t
#include <deque>

#include "variant.h"
#include "reference_genome.h"
#include "genomic_region.h"
#include "variant_factory.h"

#include <iostream> // Just for testing

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

bool is_left_aligned(const Variant& a_variant, ReferenceGenome& the_reference) noexcept
{
    return false;
}

Variant left_align(const Variant& a_variant, ReferenceGenome& the_reference)
{
    const auto& ref_allele = a_variant.get_reference_allele();
    const auto& alt_allele = a_variant.get_alternative_allele();
    
    auto max_allele_size = std::max(reference_allele_size(a_variant),
                                    alternative_allele_size(a_variant));
    
    using StringType = std::remove_reference_t<decltype(ref_allele)>;
    
    std::deque<StringType::value_type> ref_allele_q {cbegin(ref_allele), cend(ref_allele)};
    std::deque<StringType::value_type> alt_allele_q {cbegin(alt_allele), cend(alt_allele)};
    
    std::deque<StringType::value_type>::const_reverse_iterator ref_allele_reverse_itr {};
    std::deque<StringType::value_type>::const_reverse_iterator alt_allele_reverse_itr {};
    
    GenomicRegion current_ref_region = a_variant.get_reference_allele_region();
    
    while (true) {
        if (ref_allele_q.size() <= alt_allele_q.size()) {
            std::tie(ref_allele_reverse_itr, alt_allele_reverse_itr)
                    = std::mismatch(crbegin(ref_allele_q), crend(ref_allele_q), crbegin(alt_allele_q));
        } else {
            std::tie(alt_allele_reverse_itr, ref_allele_reverse_itr)
                = std::mismatch(crbegin(alt_allele_q), crend(alt_allele_q), crbegin(ref_allele_q));
        }
        
        if (ref_allele_reverse_itr == crend(ref_allele_q) ||
            alt_allele_reverse_itr == crend(alt_allele_q))
        {
            current_ref_region = shift(current_ref_region, -max_allele_size);
            auto some_reference_sequence = the_reference.get_sequence(current_ref_region);
            std::cout << some_reference_sequence << std::endl;
            ref_allele_q.insert(begin(ref_allele_q), cbegin(some_reference_sequence),
                                cend(some_reference_sequence));
            alt_allele_q.insert(begin(alt_allele_q), cbegin(some_reference_sequence),
                                cend(some_reference_sequence));
        } else {
            break;
        }
    }
    
    auto p = std::mismatch(cbegin(ref_allele_q), cend(ref_allele_q), cbegin(alt_allele_q));
    
    StringType new_ref_allele {p.first, cend(ref_allele_q)};
    StringType new_alt_allele {p.second, cend(alt_allele_q)};
    
    std::cout << new_ref_allele << std::endl;
    std::cout << new_alt_allele << std::endl;
    
    VariantFactory a_factory {};
    return a_factory.make(current_ref_region, new_ref_allele, new_alt_allele);
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
