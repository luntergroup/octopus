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
CoalescentModel::CoalescentModel(Haplotype reference_haplotype,
                                 double snp_heterozygosity,
                                 double indel_heterozygosity)
:
reference_haplotypes_ {std::move(reference_haplotype)},
snp_heterozygosity_ {snp_heterozygosity},
indel_heterozygosity_ {indel_heterozygosity}
{}

CoalescentModel::CoalescentModel(std::vector<Haplotype> reference_haplotypes,
                                 double snp_heterozygosity,
                                 double indel_heterozygosity)
:
reference_haplotypes_ {std::move(reference_haplotypes)},
snp_heterozygosity_ {snp_heterozygosity},
indel_heterozygosity_ {indel_heterozygosity}
{
    difference_cache_.reserve(128);
    difference_cache_.emplace(std::cref(reference_haplotypes_.front()), std::vector<Variant> {});
}

void CoalescentModel::set_reference(Haplotype reference)
{
    reference_haplotypes_.clear();
    reference_haplotypes_.emplace_back(std::move(reference));
    difference_cache_.clear();
    difference_cache_.reserve(128);
    difference_cache_.emplace(std::cref(reference_haplotypes_.front()), std::vector<Variant> {});
}
} // namespace Octopus
