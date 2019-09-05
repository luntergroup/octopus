// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "independent_population_model.hpp"

namespace octopus { namespace model {

IndependentPopulationModel::IndependentPopulationModel(const GenotypePriorModel& genotype_prior_model,
                                                       boost::optional<logging::DebugLogger> debug_log,
                                                       boost::optional<logging::TraceLogger> trace_log)
: individual_model_ {genotype_prior_model, debug_log, trace_log}
{}

IndependentPopulationModel::InferredLatents
IndependentPopulationModel::evaluate(const SampleVector& samples,
                                     const GenotypeVector& genotypes,
                                     const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    InferredLatents result {};
    result.posteriors.genotype_probabilities.reserve(samples.size());
    for (const auto& sample : samples) {
        haplotype_likelihoods.prime(sample);
        auto sample_results = individual_model_.evaluate(genotypes, haplotype_likelihoods);
        result.posteriors.genotype_probabilities.push_back(std::move(sample_results.posteriors.genotype_probabilities));
        result.log_evidence += sample_results.log_evidence;
    }
    return result;
}

IndependentPopulationModel::InferredLatents
IndependentPopulationModel::evaluate(const SampleVector& samples,
                                     const std::vector<unsigned>& sample_ploidies,
                                     const GenotypeVector& genotypes,
                                     const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto max_ploidy = *std::max_element(std::cbegin(sample_ploidies), std::cend(sample_ploidies));
    std::vector<GenotypeVector> genotypes_by_ploidy(max_ploidy + 1);
    for (auto& gs : genotypes_by_ploidy) gs.reserve(genotypes.size());
    for (const auto& genotype : genotypes) {
        if (genotype.ploidy() <= max_ploidy) genotypes_by_ploidy[genotype.ploidy()].push_back(genotype);
    }
    InferredLatents result {};
    for (std::size_t s {0}; s < samples.size(); ++s) {
        haplotype_likelihoods.prime(samples[s]);
        const auto sample_ploidy = sample_ploidies[s];
        auto sample_results = individual_model_.evaluate(genotypes_by_ploidy[sample_ploidy], haplotype_likelihoods);
        Latents::GenotypeProbabilityVector genotype_posteriors(genotypes.size());
        for (std::size_t genotype_idx {0}, sample_genotype_idx {0}; genotype_idx < genotypes.size(); ++genotype_idx) {
            if (genotypes[genotype_idx].ploidy() == sample_ploidy) {
                genotype_posteriors[genotype_idx] = sample_results.posteriors.genotype_probabilities[sample_genotype_idx++];
            }
        }
        result.posteriors.genotype_probabilities.push_back(std::move(genotype_posteriors));
        result.log_evidence += sample_results.log_evidence;
    }
    return result;
}
    
} // namesapce model
} // namespace octopus
