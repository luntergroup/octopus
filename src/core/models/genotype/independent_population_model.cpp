// Copyright (c) 2016 Daniel Cooke
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
                                     const HaplotypeLikelihoodCache& haplotype_likelihoods) const
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
                                     const std::vector<GenotypeVectorReference>& genotypes,
                                     const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(samples.size() == genotypes.size());
    InferredLatents result {};
    for (std::size_t s {0}; s < samples.size(); ++s) {
        haplotype_likelihoods.prime(samples[s]);
        auto sample_results = individual_model_.evaluate(genotypes[s], haplotype_likelihoods);
        result.posteriors.genotype_probabilities.push_back(std::move(sample_results.posteriors.genotype_probabilities));
        result.log_evidence += sample_results.log_evidence;
    }
    return result;
}
    
} // namesapce model
} // namespace octopus
