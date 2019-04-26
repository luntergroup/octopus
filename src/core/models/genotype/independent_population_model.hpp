// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef independent_population_model_hpp
#define independent_population_model_hpp

#include <vector>
#include <functional>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "genotype_prior_model.hpp"
#include "individual_model.hpp"
#include "core/models/haplotype_likelihood_array.hpp"
#include "containers/probability_matrix.hpp"
#include "containers/mappable_block.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace model {

class IndependentPopulationModel
{
public:
    using SampleVector            = std::vector<SampleName>;
    using GenotypeVector          = MappableBlock<Genotype<Haplotype>>;
    using GenotypeVectorReference = std::reference_wrapper<const GenotypeVector>;
    
    struct Latents
    {
        using GenotypeProbabilityVector = std::vector<double>;
        std::vector<GenotypeProbabilityVector> genotype_probabilities;
    };
    
    struct InferredLatents
    {
        Latents posteriors;
        double log_evidence;
    };
    
    IndependentPopulationModel() = delete;
    
    IndependentPopulationModel(const GenotypePriorModel& genotype_prior_model,
                               boost::optional<logging::DebugLogger> debug_log = boost::none,
                               boost::optional<logging::TraceLogger> trace_log = boost::none);
    
    IndependentPopulationModel(const IndependentPopulationModel&)            = delete;
    IndependentPopulationModel& operator=(const IndependentPopulationModel&) = delete;
    IndependentPopulationModel(IndependentPopulationModel&&)                 = delete;
    IndependentPopulationModel& operator=(IndependentPopulationModel&&)      = delete;
    
    ~IndependentPopulationModel() = default;
    
    // All samples have same ploidy
    InferredLatents evaluate(const SampleVector& samples,
                             const GenotypeVector& genotypes,
                             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    
    // Samples have different ploidy
    InferredLatents evaluate(const SampleVector& samples,
                             const std::vector<GenotypeVectorReference>& genotypes,
                             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;

private:
    IndividualModel individual_model_;
};

} // namesapce model
} // namespace octopus

#endif
