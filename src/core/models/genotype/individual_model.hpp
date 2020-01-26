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
#include "containers/mappable_block.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace model {

class IndividualModel
{
public:
    using LogProbability = double;
    
    struct Latents
    {
        using ProbabilityVector = std::vector<LogProbability>;
        ProbabilityVector genotype_log_probabilities, genotype_probabilities;
    };
    
    struct InferredLatents
    {
        Latents posteriors;
        LogProbability log_evidence;
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
    
    void prime(const MappableBlock<Haplotype>& haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    InferredLatents
    evaluate(const MappableBlock<Genotype<IndexedHaplotype<>>>& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    
private:
    const GenotypePriorModel& genotype_prior_model_;
    const MappableBlock<Haplotype>* haplotypes_;
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
    mutable boost::optional<logging::TraceLogger> trace_log_;
};

} // namesapce model
} // namespace octopus

#endif
