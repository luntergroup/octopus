//
//  variant_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include <functional>

#include "variant_factory.h"

Variant
VariantFactory::make(std::string contig_name, __uint32_t contig_begin_pos,
                     std::string sequence_added, std::string sequence_removed) const
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
    return Variant(contig_name, contig_begin_pos, sequence_added, sequence_removed, prior_model);
}
