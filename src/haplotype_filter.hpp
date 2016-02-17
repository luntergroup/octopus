//
//  haplotype_filter.hpp
//  Octopus
//
//  Created by Daniel Cooke on 22/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef haplotype_filter_hpp
#define haplotype_filter_hpp

#include <vector>
#include <cstddef>
#include <functional>

#include "common.hpp"
#include "haplotype.hpp"
#include "haplotype_likelihood_cache.hpp"

namespace Octopus
{
    void filter_haplotypes(std::vector<Haplotype>& haplotypes, const ReadMap& reads, size_t n,
                           const HaplotypeLikelihoodCache& haplotype_likelihoods);
    
    std::vector<Haplotype>
    copy_filtered_haplotypes(const std::vector<Haplotype>& haplotypes,
                             const ReadMap& reads, size_t n,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods);
}

#endif /* haplotype_filter_hpp */
