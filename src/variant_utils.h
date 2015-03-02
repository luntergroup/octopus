//
//  variant_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_utils__
#define __Octopus__variant_utils__

class Variant;
class ReferenceGenome;

Variant left_align(const Variant& a_variant, ReferenceGenome& the_reference);

#endif /* defined(__Octopus__variant_utils__) */
