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

#include "common.hpp"
#include "haplotype.hpp"

namespace Octopus {
    std::vector<Haplotype> filter_haplotypes(const std::vector<Haplotype>& haplotypes, const ReadMap& reads, size_t n);
} // namespace Octopus

#endif /* haplotype_filter_hpp */
