// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef population_model_hpp
#define population_model_hpp

#include <vector>
#include <functional>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "population_prior_model.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"
#include "containers/probability_matrix.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace model {

class PopulationModel
{
public:
    struct Latents
    {
        std::vector<std::vector<std::size_t>> genotype_combinations;
        using GenotypeProbabilityVector = std::vector<double>;
        GenotypeProbabilityVector joint_genotype_probabilities;
    };
    
    struct InferredLatents
    {
        Latents posteriors;
        double log_evidence;
        bool overflowed = false;
    };
    
    struct Options
    {
        std::size_t min_to_keep = 50, max_to_keep = 500;
        double max_removal_posterior_mass = 1e-20;
    };
    
    using SampleVector            = std::vector<SampleName>;
    using GenotypeVector          = std::vector<Genotype<Haplotype>>;
    using GenotypeVectorReference = std::reference_wrapper<const GenotypeVector>;
    
    PopulationModel() = delete;
    
    PopulationModel(const PopulationPriorModel& prior_model,
                    boost::optional<logging::DebugLogger> debug_log = boost::none);
    
    PopulationModel(const PopulationModel&)            = delete;
    PopulationModel& operator=(const PopulationModel&) = delete;
    PopulationModel(PopulationModel&&)                 = delete;
    PopulationModel& operator=(PopulationModel&&)      = delete;
    
    ~PopulationModel() = default;
    
    // All samples have same ploidy
    InferredLatents evaluate(const SampleVector& samples,
                             const GenotypeVector& genotypes,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    // Samples have different ploidy
    InferredLatents evaluate(const SampleVector& samples,
                             const std::vector<GenotypeVectorReference>& genotypes,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
private:
    const PopulationPriorModel& prior_model_;
    mutable boost::optional<logging::DebugLogger> debug_log_;
};

} // namesapce model
} // namespace octopus

#endif
