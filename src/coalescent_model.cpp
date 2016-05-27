//
//  coalescent_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "coalescent_model.hpp"

namespace Octopus
{
CoalescentModel::CoalescentModel(Haplotype reference,
                                 double snp_heterozygosity,
                                 double indel_heterozygosity)
:
reference_ {reference},
snp_heterozygosity_ {snp_heterozygosity},
indel_heterozygosity_ {indel_heterozygosity}
{
    difference_cache_.reserve(1024);
    difference_cache_.emplace(std::piecewise_construct, std::forward_as_tuple(reference),
                              std::forward_as_tuple());
}

void CoalescentModel::set_reference(Haplotype reference)
{
    reference_ = std::move(reference_);
    difference_cache_.clear();
    difference_cache_.reserve(1024);
    difference_cache_.emplace(std::piecewise_construct, std::forward_as_tuple(reference),
                              std::forward_as_tuple());
}
} // namespace Octopus
