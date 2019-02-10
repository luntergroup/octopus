// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef individual_model_hpp
#define individual_model_hpp

#include <vector>

#include <boost/optional.hpp>

#include "genotype_prior_model.hpp"
#include "core/types/haplotype.hpp"
#include "core/models/haplotype_likelihood_array.hpp"
#include "core/types/genotype.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace model {

class IndividualModel
{
public:
    struct Latents
    {
        using ProbabilityVector = std::vector<double>;
        ProbabilityVector genotype_probabilities;
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
    
    const GenotypePriorModel& prior_model() const noexcept;
    
    void prime(const std::vector<Haplotype>& haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    InferredLatents evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
                             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    
    InferredLatents evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
                             const std::vector<GenotypeIndex>& genotype_indices,
                             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    
private:
    const GenotypePriorModel& genotype_prior_model_;
    const std::vector<Haplotype>* haplotypes_;
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
    mutable boost::optional<logging::TraceLogger> trace_log_;
};

} // namesapce model
} // namespace octopus

#endif
