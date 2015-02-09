//
//  variant_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include <functional>

#include "variant_factory.h"

std::unique_ptr<Variant>
VariantFactory::make(GenomicRegion ref_region, std::string sequence_added,
                     std::string sequence_removed) const
{
    std::function<double()> prior_model {};
    if (sequence_added.length() == sequence_removed.length()) {
        if (sequence_added.length() == 1) {
            prior_model = [] () { return 1e-5; };
        } else {
            prior_model = [] () { return 1e-6; };
        }
    } else {
        if (sequence_added.length() < sequence_removed.length()) {
            prior_model = [] () { return 1e-7; };
        } else {
            prior_model = [] () { return 1e-8; };
        }
    }
    return std::make_unique<Variant> (ref_region, sequence_added, sequence_removed, prior_model);
}
