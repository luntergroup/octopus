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
#include <cstddef>   // std::size_t
#include "variant.h"

class ReferenceGenome;
class VariantFactory;
class IVariantCandidateGenerator;

using VariantGenerators = std::vector<IVariantCandidateGenerator>;

//double get_variant_prior_probability(const Variant& a_variant,
//                                     VariantCandidateGenerator& a_variant_candidate_generator,
//                                     );

/*
 A variant is parsimonious if and only if it is represented in as few nucleotides as possible
 without an allele of length 0.
 */
bool is_parsimonious(const Variant& a_variant) noexcept;

Variant make_parsimonious(const Variant& a_variant, ReferenceGenome& the_reference,
                          VariantFactory& a_variant_factory);

bool is_left_alignable(const Variant& a_variant) noexcept;

// Note there are no 'is_left_aligned' or 'is_normalised' methods because they
// would take as much work as calling 'left_align'.

/*
 A variant is left aligned if and only if it is no longer possible to shift its position
 to the left while keeping the length of all its alleles constant.
 */
Variant left_align(const Variant& a_variant, ReferenceGenome& the_reference,
                   VariantFactory& a_variant_factory, Variant::SizeType extension_size=30);

/*
  A variant is normalised if and only if it is parsimonious and left aligned.
 */
Variant normalise(const Variant& a_variant, ReferenceGenome& the_reference,
                  VariantFactory& a_variant_factory, Variant::SizeType extension_size=30);

bool is_snp(const Variant& a_variant) noexcept;

bool is_insertion(const Variant& a_variant) noexcept;

bool is_deletion(const Variant& a_variant) noexcept;

bool is_indel(const Variant& a_variant) noexcept;

bool is_mnv(const Variant& a_variant) noexcept;

template <typename Iterator>
std::size_t num_redundant_bases(Iterator first1, Iterator last1, Iterator first2)
{
    auto matching_begin_itrs = std::mismatch(first1, last1, first2);
    return std::distance(first1, matching_begin_itrs.first);
}

#endif /* defined(__Octopus__variant_utils__) */
