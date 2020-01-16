// Copyright (c) 2015-2019 Daniel Cooke
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
#include "core/models/haplotype_likelihood_array.hpp"
#include "containers/probability_matrix.hpp"
#include "containers/mappable_block.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace model {

class PopulationModel
{
public:
    struct Options
    {
        boost::optional<std::size_t> max_joint_genotypes = boost::none;
        unsigned max_em_iterations = 100;
        double em_epsilon = 0.001;
    };
    struct Latents
    {
        using ProbabilityVector = std::vector<double>;
        std::vector<ProbabilityVector> marginal_genotype_probabilities;
    };
    struct InferredLatents
    {
        Latents posteriors;
        double log_evidence;
    };
    
    using SampleVector            = std::vector<SampleName>;
    using GenotypeVector          = MappableBlock<Genotype<Haplotype>>;
    using GenotypeVectorReference = std::reference_wrapper<const GenotypeVector>;
    
    PopulationModel() = delete;
    
    PopulationModel(const PopulationPriorModel& prior_model,
                    boost::optional<logging::DebugLogger> debug_log = boost::none);
    PopulationModel(const PopulationPriorModel& prior_model,
                    Options options,
                    boost::optional<logging::DebugLogger> debug_log = boost::none);
    
    PopulationModel(const PopulationModel&)            = delete;
    PopulationModel& operator=(const PopulationModel&) = delete;
    PopulationModel(PopulationModel&&)                 = delete;
    PopulationModel& operator=(PopulationModel&&)      = delete;
    
    ~PopulationModel() = default;
    
    const PopulationPriorModel& prior_model() const noexcept;
    
    // All samples have same ploidy
    InferredLatents
    evaluate(const SampleVector& samples,
             const GenotypeVector& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    // All samples have same ploidy
    InferredLatents
    evaluate(const SampleVector& samples,
             const GenotypeVector& genotypes,
             const std::vector<GenotypeIndex>& genotype_indices,
             const std::vector<Haplotype>& haplotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    // Samples have different ploidy
    InferredLatents
    evaluate(const SampleVector& samples,
             const std::vector<unsigned>& sample_ploidies,
             const GenotypeVector& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    InferredLatents
    evaluate(const SampleVector& samples,
             const std::vector<unsigned>& sample_ploidies,
             const GenotypeVector& genotypes,
             const std::vector<GenotypeIndex>& genotype_indices,
             const std::vector<Haplotype>& haplotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    
private:
    Options options_;
    const PopulationPriorModel& prior_model_;
    mutable boost::optional<logging::DebugLogger> debug_log_;
};

} // namesapce model
} // namespace octopus

#endif
