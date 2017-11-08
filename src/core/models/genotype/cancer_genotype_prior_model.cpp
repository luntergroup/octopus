// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cancer_genotype_prior_model.hpp"

#include <utility>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "utils/mappable_algorithms.hpp"

namespace octopus {

CancerGenotypePriorModel::CancerGenotypePriorModel(const GenotypePriorModel& germline_model,
                                                   SomaticMutationModel mutation_model)
: germline_model_ {germline_model}
, mutation_model_ {mutation_model}
{}

double CancerGenotypePriorModel::evaluate(const CancerGenotype<Haplotype>& genotype) const
{
    const auto& germline = genotype.germline_genotype();
    const auto& somatic = genotype.somatic_element();
    const auto germline_log_prior = germline_model_.get().evaluate(germline);
    return germline_log_prior + ln_probability_of_somatic(somatic, germline);
}

// p(somatic | germline) = 1 / M sum [k = 1 -> M] p(somatic | germline_k) (M = germline ploidy)
double CancerGenotypePriorModel::ln_probability_of_somatic(const Haplotype& somatic, const Genotype<Haplotype>& germline) const
{
    switch (germline.ploidy()) {
        case 1: return ln_probability_of_somatic(somatic, germline);
        case 2:
        {
            const static double ln2 {std::log(2)};
            const auto a = ln_probability_of_somatic(somatic, germline[0]);
            const auto b = ln_probability_of_somatic(somatic, germline[1]);
            return maths::log_sum_exp(a, b) - ln2;
        }
        case 3:
        {
            const static double ln3 {std::log(2)};
            const auto a = ln_probability_of_somatic(somatic, germline[0]);
            const auto b = ln_probability_of_somatic(somatic, germline[1]);
            const auto c = ln_probability_of_somatic(somatic, germline[3]);
            return maths::log_sum_exp(a, b, c) - ln3;
        }
        default:
        {
            std::vector<double> tmp(germline.ploidy());
            std::transform(std::cbegin(germline), std::cend(germline), std::begin(tmp),
                           [this, &somatic] (const Haplotype& haplotype) {
                               return ln_probability_of_somatic(somatic, haplotype);
                           });
            return maths::log_sum_exp(tmp) - std::log(germline.ploidy());
        }
    }
}

double CancerGenotypePriorModel::ln_probability_of_somatic(const Haplotype& somatic, const Haplotype& germline) const
{
    return mutation_model_.evaluate(somatic, germline);
}
    
} // namespace octopus
