//
//  haplotype_prior_model_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "haplotype_prior_model_factory.hpp"

#include "basic_haplotype_prior_model.hpp"

namespace Octopus
{
    std::unique_ptr<HaplotypePriorModel>
    HaplotypePriorModelFactory::make(const ReferenceGenome& reference) const
    {
        return std::make_unique<BasicHaplotypePriorModel>(reference);
    }
} // namespace Octopus
