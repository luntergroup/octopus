// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_mutation_model.hpp"

#include <utility>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "utils/mappable_algorithms.hpp"

namespace octopus {

SomaticMutationModel::SomaticMutationModel(Parameters params, std::size_t num_haplotypes_hint,
                                           CachingStrategy caching)
: model_ {params, num_haplotypes_hint, caching}
{}

void SomaticMutationModel::prime(MappableBlock<Haplotype> haplotypes)
{
    model_.prime(std::move(haplotypes));
}

void SomaticMutationModel::unprime() noexcept
{
    model_.unprime();
}

bool SomaticMutationModel::is_primed() const noexcept
{
    return model_.is_primed();
}

SomaticMutationModel::LogProbability SomaticMutationModel::evaluate(const Haplotype& somatic, const Haplotype& germline) const
{
    return model_.evaluate(somatic, germline);
}

SomaticMutationModel::LogProbability SomaticMutationModel::evaluate(unsigned somatic, unsigned germline) const
{
    return model_.evaluate(somatic, germline);
}

} // namespace octopus
