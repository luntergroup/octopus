// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef family_model_hpp
#define family_model_hpp

#include <vector>
#include <functional>
#include <cstddef>

#include <boost/optional.hpp>

#include "basics/pedigree.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/models/haplotype_likelihood_array.hpp"
#include "core/types/genotype.hpp"
#include "logging/logging.hpp"
#include "pedigree_prior_model.hpp"

namespace octopus { namespace model {

class FamilyModel
{
public:
    using SampleVector   = std::vector<SampleName>;
    using GenotypeVector = MappableBlock<Genotype<IndexedHaplotype<>>>;
    
    struct Latents
    {
        using GenotypeCombination = std::vector<unsigned>;
        struct GenotypeCombinationProbability
        {
            GenotypeCombination combination;
            double log_probability;
        };
        using JointProbabilityVector = std::vector<GenotypeCombinationProbability>;
        JointProbabilityVector joint_genotype_probabilities;
    };
    
    struct InferredLatents
    {
        Latents posteriors;
        double log_evidence;
    };
    
    struct Options
    {
        boost::optional<std::size_t> max_genotype_combinations = boost::none;
    };
    
    FamilyModel() = delete;
    
    FamilyModel(const Pedigree& family,
                const PedigreePriorModel& prior_model,
                Options options,
                boost::optional<logging::DebugLogger> debug_log = boost::none);
    
    FamilyModel(const FamilyModel&)            = delete;
    FamilyModel& operator=(const FamilyModel&) = delete;
    FamilyModel(FamilyModel&&)                 = delete;
    FamilyModel& operator=(FamilyModel&&)      = delete;
    
    ~FamilyModel() = default;
    
    static unsigned max_ploidy() noexcept;
    
    const PedigreePriorModel& prior_model() const noexcept;
    
    // All samples have same ploidy
    InferredLatents
    evaluate(const SampleVector& samples,
             const MappableBlock<Haplotype>& haplotypes,
             const GenotypeVector& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    // Samples have different ploidy
    InferredLatents
    evaluate(const SampleVector& samples,
             const std::vector<unsigned>& sample_ploidies,
             const MappableBlock<Haplotype>& haplotypes,
             const GenotypeVector& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    
private:
    const Pedigree& family_;
    const PedigreePriorModel& prior_model_;
    Options options_;
    mutable boost::optional<logging::DebugLogger> debug_log_;
};
    
} // namespace model
} // namespace octopus

#endif
