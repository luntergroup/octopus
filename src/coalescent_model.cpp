//
//  coalescent_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "coalescent_model.hpp"

#include <stdexcept>

namespace Octopus
{
CoalescentModel::CoalescentModel(Haplotype reference,
                                 double snp_heterozygosity,
                                 double indel_heterozygosity,
                                 unsigned max_haplotypes)
:
reference_ {reference},
snp_heterozygosity_ {snp_heterozygosity},
indel_heterozygosity_ {indel_heterozygosity}
{
    if (snp_heterozygosity <= 0 || indel_heterozygosity) {
        throw std::domain_error {"CoalescentModel: snp and indel heterozygosity must be > 0"};
    }
    
    site_buffer1_.reserve(128);
    site_buffer2_.reserve(128);
    difference_cache_.reserve(max_haplotypes);
    difference_cache_.emplace(std::piecewise_construct, std::forward_as_tuple(reference),
                              std::forward_as_tuple());
    result_cache_.reserve(max_haplotypes);
}

void CoalescentModel::set_reference(Haplotype reference)
{
    reference_ = std::move(reference_);
    difference_cache_.clear();
    difference_cache_.emplace(std::piecewise_construct, std::forward_as_tuple(reference),
                              std::forward_as_tuple());
}
} // namespace Octopus
