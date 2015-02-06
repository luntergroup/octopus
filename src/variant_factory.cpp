//
//  variant_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_factory.h"

#include "variant_prior_model.h"

std::unique_ptr<Variant>
VariantFactory::make(GenomeRegion ref_region, std::string sequence_added,
                     std::string sequence_removed) const
{
    if (sequence_added.length() == sequence_removed.length()) {
        if (sequence_added.length() == 1) {
            
        } else {
            return std::make_unsigned<Variant> {};
        }
    }
    return std::make_unsigned<Variant> {ref_region, sequence_added, sequence_removed, };
}
