//
//  somatic_call.cpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "somatic_call.hpp"

#include "string_utils.hpp"

namespace Octopus
{
    void SomaticCall::decorate(VcfRecord::Builder& record) const
    {
        record.set_somatic();
        
        record.set_alt_allele(variant_.get_alt_allele().get_sequence());
        
        record.add_format("SCR");
        
        for (const auto& p : credible_regions_) {
            record.add_genotype_field(p.first, "SCR", {
                Octopus::to_string(p.second.somatic_credible_region.first, 2),
                Octopus::to_string(p.second.somatic_credible_region.second, 2)
            });
        }
    }
} // namespace Octopus
