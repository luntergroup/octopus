// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef individual_model_hpp
#define individual_model_hpp

#include <vector>

#include <boost/optional.hpp>

#include "genotype_prior_model.hpp"
#include "core/types/haplotype.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"
#include "core/types/genotype.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace model {

class IndividualModel
{
public:
    struct Latents
    {
        using GenotypeProbabilityVector = std::vector<double>;
        GenotypeProbabilityVector genotype_probabilities;
    };
    
    struct InferredLatents
    {
        Latents posteriors;
        double log_evidence;
    };
    
    IndividualModel() = delete;
    
    IndividualModel(const GenotypePriorModel& genotype_prior_model,
                    boost::optional<logging::DebugLogger> debug_log = boost::none,
                    boost::optional<logging::TraceLogger> trace_log = boost::none);
    
    IndividualModel(const IndividualModel&)            = delete;
    IndividualModel& operator=(const IndividualModel&) = delete;
    IndividualModel(IndividualModel&&)                 = delete;
    IndividualModel& operator=(IndividualModel&&)      = delete;
    
    ~IndividualModel() = default;
    
    InferredLatents evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
private:
    const GenotypePriorModel& genotype_prior_model_;
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
    mutable boost::optional<logging::TraceLogger> trace_log_;
};

} // namesapce model
} // namespace octopus

#endif
