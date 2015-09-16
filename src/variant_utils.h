//
//  variant_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_utils__
#define __Octopus__variant_utils__

#include <vector>
#include <algorithm> // std::mismatch, std::max, std::minmax
#include <iterator>  // std::distance
#include <cstddef>   // std::size_t

#include "allele.h"
#include "variant.h"
#include "genomic_region.h"

class ReferenceGenome;
class VariantFactory;

void remove_duplicates(std::vector<Variant>& variants);

std::vector<Allele> decompose(const std::vector<Variant>& variants);

std::vector<Allele> get_intervening_reference_alleles(const std::vector<Variant>& variants,
                                                      ReferenceGenome& reference);

/*
 A variant is parsimonious if and only if it is represented in as few nucleotides as possible
 without an allele of length 0.
 */
bool is_parsimonious(const Variant& variant) noexcept;

Variant make_parsimonious(const Variant& variant, ReferenceGenome& reference);

bool is_left_alignable(const Variant& variant) noexcept;

// Note there are no 'is_left_aligned' or 'is_normalised' methods because they
// would take as much work as calling 'left_align'.

/*
 A variant is left aligned if and only if it is no longer possible to shift its position
 to the left while keeping the length of all its alleles constant.
 */
Variant left_align(const Variant& variant, ReferenceGenome& reference,
                   Variant::SizeType extension_size = 30);

/*
  A variant is normalised if and only if it is parsimonious and left aligned.
 */
Variant normalise(const Variant& variant, ReferenceGenome& reference,
                  Variant::SizeType extension_size = 30);

/*
 Left aligns all input Variants and removes any resulting duplicates. The returned variants are sorted.
 */
std::vector<Variant> unique_left_align(const std::vector<Variant>& variants, ReferenceGenome& reference);

bool is_snp(const Variant& variant) noexcept;

bool is_insertion(const Variant& variant) noexcept;

bool is_deletion(const Variant& variant) noexcept;

bool is_indel(const Variant& variant) noexcept;

bool is_mnv(const Variant& variant) noexcept;

bool is_transition(const Variant& variant) noexcept;

bool is_transversion(const Variant& variant) noexcept;

#endif /* defined(__Octopus__variant_utils__) */
