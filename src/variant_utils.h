//
//  variant_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_utils__
#define __Octopus__variant_utils__

#include "variant.h"

class ReferenceGenome;

bool is_parsimonious(const Variant& a_variant) noexcept;

Variant make_parsimonious(const Variant& a_variant);

bool is_left_alignable(const Variant& a_variant) noexcept;

Variant left_align(const Variant& a_variant, ReferenceGenome& the_reference,
                   Variant::SizeType extension_size=10);

bool is_normalised(const Variant& a_variant) noexcept;

Variant normalise(const Variant& a_variant);

bool is_snp(const Variant& a_variant) noexcept;

bool is_insertion(const Variant& a_variant) noexcept;

bool is_deletion(const Variant& a_variant) noexcept;

bool is_indel(const Variant& a_variant) noexcept;

bool is_mnv(const Variant& a_variant) noexcept;

#endif /* defined(__Octopus__variant_utils__) */
