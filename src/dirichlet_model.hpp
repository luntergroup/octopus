//
//  dirichlet_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef dirichlet_model_hpp
#define dirichlet_model_hpp

#include <unordered_map>
#include <functional>
#include <iterator>
#include <algorithm>
#include <numeric>

#include "haplotype.hpp"
#include "maths.hpp"

namespace Octopus
{
namespace GenotypeModel
{
    using HaplotypePriorCountMap = std::unordered_map<std::reference_wrapper<const Haplotype>, double>;
    using HaplotypeFrequencyMap  = std::unordered_map<std::reference_wrapper<const Haplotype>, double>;
    
    HaplotypeFrequencyMap init_haplotype_frequencies(const std::vector<Haplotype>& haplotypes);
    
    HaplotypeFrequencyMap init_haplotype_frequencies(const HaplotypePriorCountMap& haplotype_counts);
    
    HaplotypePriorCountMap compute_haplotype_prior_counts(const HaplotypeFrequencyMap& haplotype_priors);
} // namespace GenotypeModel
} // namespace Octopus

#endif /* dirichlet_model_hpp */
