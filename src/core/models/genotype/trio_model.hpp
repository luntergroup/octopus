// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef trio_model_hpp
#define trio_model_hpp

#include <vector>
#include <functional>
#include <cstddef>

#include <boost/optional.hpp>

#include "core/types/trio.hpp"
#include "core/types/haplotype.hpp"
#include "core/models/mutation/coalescent_model.hpp"
#include "core/models/mutation/denovo_model.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"
#include "core/types/genotype.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace model {

class TrioModel
{
public:
    using GenotypeVector = std::vector<Genotype<Haplotype>>;
    
    struct Latents
    {
        struct JointProbability
        {
            std::reference_wrapper<const Genotype<Haplotype>> maternal, paternal, child;
            double probability;
        };
        std::vector<JointProbability> joint_genotype_probabilities;
    };
    
    struct InferredLatents
    {
        Latents posteriors;
        double log_evidence;
    };
    
    TrioModel() = delete;
    
    TrioModel(const Trio& trio,
              const CoalescentModel& genotype_prior_model,
              const DeNovoModel& mutation_model,
              boost::optional<logging::DebugLogger> debug_log = boost::none);
    
    TrioModel(const TrioModel&)            = delete;
    TrioModel& operator=(const TrioModel&) = delete;
    TrioModel(TrioModel&&)                 = delete;
    TrioModel& operator=(TrioModel&&)      = delete;
    
    ~TrioModel() = default;
    
    InferredLatents evaluate(const GenotypeVector& maternal_genotypes,
                             const GenotypeVector& paternal_genotypes,
                             const GenotypeVector& child_genotypes,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;

private:
    const Trio& trio_;
    const CoalescentModel& genotype_prior_model_;
    const DeNovoModel& mutation_model_;
    
    std::size_t max_search_size_ = 100;
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
};
    
} // namespace model
} // namespace octopus

#endif
