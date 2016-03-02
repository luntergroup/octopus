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

class Haplotype;
class HaplotypeLikelihoodCache;
class ReadMap;

namespace Octopus
{
    std::vector<Haplotype> filter_n_haplotypes(std::vector<Haplotype>& haplotypes,
                                               const ReadMap& reads,
                                               const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                               const std::size_t n);
} // namespace Octopus

#endif /* haplotype_filter_hpp */
