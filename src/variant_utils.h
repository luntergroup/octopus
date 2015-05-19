//
//  variant_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_utils__
#define __Octopus__variant_utils__

#include <algorithm> // std::mismatch, std::max, std::minmax
#include <iterator>  // std::distance
#include <cstddef>   // std::size_t

#include "allele.h"
#include "variant.h"
#include "genomic_region.h"

class ReferenceGenome;
class VariantFactory;

void merge_equal_variants(std::vector<Variant>& the_variants);

std::vector<Allele> get_reference_alleles_between_variants(const std::vector<Variant>& the_variants,
                                                           ReferenceGenome& the_reference);

/*
 A variant is parsimonious if and only if it is represented in as few nucleotides as possible
 without an allele of length 0.
 */
bool is_parsimonious(const Variant& a_variant) noexcept;

Variant make_parsimonious(const Variant& a_variant, ReferenceGenome& the_reference);

bool is_left_alignable(const Variant& a_variant) noexcept;

// Note there are no 'is_left_aligned' or 'is_normalised' methods because they
// would take as much work as calling 'left_align'.

/*
 A variant is left aligned if and only if it is no longer possible to shift its position
 to the left while keeping the length of all its alleles constant.
 */
Variant left_align(const Variant& a_variant, ReferenceGenome& the_reference,
                   Variant::SizeType extension_size=30);

/*
  A variant is normalised if and only if it is parsimonious and left aligned.
 */
Variant normalise(const Variant& a_variant, ReferenceGenome& the_reference,
                  Variant::SizeType extension_size=30);

bool is_snp(const Variant& a_variant) noexcept;

bool is_insertion(const Variant& a_variant) noexcept;

bool is_deletion(const Variant& a_variant) noexcept;

bool is_indel(const Variant& a_variant) noexcept;

bool is_mnv(const Variant& a_variant) noexcept;

namespace detail
{
    template <typename InputIterator>
    std::size_t num_redundant_bases(InputIterator first1, InputIterator last1, InputIterator first2)
    {
        return std::distance(first1, std::mismatch(first1, last1, first2).first);
    }
} // end namespace detail

#endif /* defined(__Octopus__variant_utils__) */
