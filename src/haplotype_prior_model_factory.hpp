//
//  haplotype_prior_model_factory.hpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef haplotype_prior_model_factory_hpp
#define haplotype_prior_model_factory_hpp

#include <memory>

#include "haplotype_prior_model.hpp"
#include "reference_genome.hpp"

namespace Octopus
{
    class HaplotypePriorModelFactory
    {
    public:
        HaplotypePriorModelFactory() = default;
        ~HaplotypePriorModelFactory() = default;
        
        HaplotypePriorModelFactory(const HaplotypePriorModelFactory&)            = default;
        HaplotypePriorModelFactory& operator=(const HaplotypePriorModelFactory&) = default;
        HaplotypePriorModelFactory(HaplotypePriorModelFactory&&)                 = default;
        HaplotypePriorModelFactory& operator=(HaplotypePriorModelFactory&&)      = default;
        
        std::unique_ptr<HaplotypePriorModel> make(const ReferenceGenome& reference) const;
        
    private:
    };
} // namespace Octopus

#endif /* haplotype_prior_model_factory_hpp */
