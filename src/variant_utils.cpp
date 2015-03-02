//
//  variant_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_utils.h"

#include <algorithm>

#include "variant.h"
#include "reference_genome.h"

using std::crbegin;
using std::crend;

Variant left_align(const Variant& a_variant, ReferenceGenome& the_reference)
{
    auto ref_allele = a_variant.get_reference_allele();
    auto alt_allele = a_variant.get_alternative_allele();
    
    auto it = std::mismatch(crbegin(ref_allele), crend(ref_allele), crbegin(alt_allele));
    
    return a_variant;
}
