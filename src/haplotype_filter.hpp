//
//  haplotype_filter.hpp
//  Octopus
//
//  Created by Daniel Cooke on 02/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef haplotype_filter_hpp
#define haplotype_filter_hpp

#include <vector>
#include <cstddef>
#include <functional>
#include <unordered_map>

#include "common.hpp"
#include "haplotype.hpp"

namespace octopus {

class HaplotypeLikelihoodCache;

std::vector<Haplotype>
filter_to_n(std::vector<Haplotype>& haplotypes, const std::vector<SampleName>& samples,
            const HaplotypeLikelihoodCache& haplotype_likelihoods, const std::size_t n);

using HaplotypeReference    = std::reference_wrapper<const Haplotype>;
using HaplotypePosteriorMap = std::unordered_map<HaplotypeReference, double>;

std::vector<HaplotypeReference>
extract_removable(const std::vector<Haplotype>& haplotypes,
                  const HaplotypePosteriorMap& haplotype_posteriors,
                  const std::vector<SampleName>& samples,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods,
                  std::size_t max_to_remove, double min_posterior);

} // namespace octopus

#endif /* haplotype_filter_hpp */
