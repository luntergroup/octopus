// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_mutation_model.hpp"

#include <utility>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "utils/mappable_algorithms.hpp"

namespace octopus {

SomaticMutationModel::SomaticMutationModel(Parameters params)
: params_ {params}
{
    if (params_.somatic_mutation_rate <= 0) {
        throw std::domain_error {"SomaticMutationModel: somatic mutation rate must be > 0"};
    }
}

double SomaticMutationModel::evaluate(const Haplotype& somatic, const Haplotype& germline) const
{
    // TODO: implement a proper model for this (snv/indel).
    const auto mutations = difference(somatic, germline);
    const auto num_mutation_sites = count_mutually_exclusive_regions(mutations);
    return num_mutation_sites * std::log(params_.somatic_mutation_rate);
}

} // namespace octopus
