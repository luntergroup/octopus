// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variational_bayes_mixture_mixture_model_hpp
#define variational_bayes_mixture_mixture_model_hpp

#include <vector>

#include <boost/optional.hpp>

#include "core/models/haplotype_likelihood_array.hpp"
#include "utils/parallel_transform.hpp"
#include "variational_bayes_mixture_model.hpp"

namespace octopus { namespace model {

class VariationalBayesMixtureMixtureModel
{
public:
    struct Options
    {
        double null_log_probability = -10'000;
        double epsilon = 0.05;
        unsigned max_iterations = 1000;
        double save_memory = false;
        bool parallel_execution = false;
    };
    
    using Probability = double;
    using LogProbability = double;
    using ProbabilityVector = std::vector<Probability>;
    using LogProbabilityVector = std::vector<LogProbability>;
    
    using HaplotypeLikelihoodVector = std::vector<VBReadLikelihoodArray>; // One element per haplotype in genotype
    using GenotypeLikelihoodVector = std::vector<HaplotypeLikelihoodVector>; // One element per genotype in combination (i.e. num groups)
    using GenotypeCombinationLikelihoodVector = std::vector<GenotypeLikelihoodVector>; // One element per genotype combination
    using HaplotypeLikelihoodMatrix = std::vector<GenotypeCombinationLikelihoodVector>; // One element per sample
    
    using GroupConcentrationVector = std::vector<double>; // One element per group
    
    using ComponentConcentrationVector = std::vector<double>; // One element per component in group
    using MixtureConcentrationVector = std::vector<ComponentConcentrationVector>; // One element per group
    using MixtureConcentrationArray = std::vector<MixtureConcentrationVector>; // One element per sample
    
    using Sigma = std::vector<double>; // One element per group
    using GroupResponsibilityVector = std::vector<Sigma>; // One element per sample
    
    using GroupOptionalPriorVector = boost::optional<ProbabilityVector>; // One element per group
    using GroupOptionalPriorArray = std::vector<GroupOptionalPriorVector>; // One element per sample
    
    struct Inferences
    {
        ProbabilityVector genotype_posteriors;
        LogProbabilityVector genotype_log_posteriors;
        GroupConcentrationVector group_concentrations;
        MixtureConcentrationArray mixture_concentrations;
        GroupResponsibilityVector group_responsibilities;
        double approx_log_evidence;
    };
    
    VariationalBayesMixtureMixtureModel() = default;
    VariationalBayesMixtureMixtureModel(Options options);
    
    VariationalBayesMixtureMixtureModel(const VariationalBayesMixtureMixtureModel&)            = default;
    VariationalBayesMixtureMixtureModel& operator=(const VariationalBayesMixtureMixtureModel&) = default;
    VariationalBayesMixtureMixtureModel(VariationalBayesMixtureMixtureModel&&)                 = default;
    VariationalBayesMixtureMixtureModel& operator=(VariationalBayesMixtureMixtureModel&&)      = default;
    
    ~VariationalBayesMixtureMixtureModel() = default;
    
    Inferences
    evaluate(const LogProbabilityVector& genotype_log_priors,
             const HaplotypeLikelihoodMatrix& log_likelihoods,
             const GroupOptionalPriorArray& group_priors,
             const GroupConcentrationVector& group_concentrations,
             const MixtureConcentrationArray& mixture_concentrations,
             std::vector<LogProbabilityVector> seeds) const;
    
    Inferences
    evaluate(const LogProbabilityVector& genotype_log_priors,
             const HaplotypeLikelihoodMatrix& log_likelihoods,
             const GroupConcentrationVector& group_concentrations,
             const MixtureConcentrationArray& mixture_concentrations,
             std::vector<LogProbabilityVector> seeds) const;
    
    Inferences
    evaluate(const LogProbabilityVector& genotype_log_priors,
             const HaplotypeLikelihoodMatrix& log_likelihoods,
             const GroupOptionalPriorArray& group_priors,
             double group_concentration,
             const std::vector<double>& mixture_concentrations,
             std::vector<LogProbabilityVector> seeds) const;
    
    Inferences
    evaluate(const LogProbabilityVector& genotype_log_priors,
             const HaplotypeLikelihoodMatrix& log_likelihoods,
             double group_concentration,
             const std::vector<double>& mixture_concentrations,
             std::vector<LogProbabilityVector> seeds) const;
    
    Inferences
    evaluate(const LogProbabilityVector& genotype_log_priors,
             const HaplotypeLikelihoodMatrix& log_likelihoods,
             const GroupOptionalPriorArray& group_priors,
             double group_concentration,
             double mixture_concentration,
             std::vector<LogProbabilityVector> seeds) const;
    
    Inferences
    evaluate(const LogProbabilityVector& genotype_log_priors,
             const HaplotypeLikelihoodMatrix& log_likelihoods,
             double group_concentration,
             double mixture_concentration,
             std::vector<LogProbabilityVector> seeds) const;
    
private:
    using Tau = std::vector<double>; // One element per read
    using ComponentResponsibilityVector = std::vector<Tau>; // One element per haplotype in genotype (max)
    using ComponentResponsibilityMatrix = std::vector<ComponentResponsibilityVector>; // One element per sample
    
    using GroupOptionalLogPriorVector = boost::optional<LogProbabilityVector>; // One element per group
    using GroupOptionalLogPriorArray = std::vector<GroupOptionalLogPriorVector>; // One element per sample
    
    Options options_;
    
    GroupOptionalLogPriorVector to_logs(const GroupOptionalPriorVector& prior) const;
    GroupOptionalLogPriorArray to_logs(const GroupOptionalPriorArray& priors) const;
    Inferences
    evaluate(const LogProbabilityVector& genotype_log_priors,
             const HaplotypeLikelihoodMatrix& log_likelihoods,
             const GroupOptionalLogPriorArray& group_log_priors,
             const GroupConcentrationVector& group_concentrations,
             const MixtureConcentrationArray& mixture_concentrations,
             LogProbabilityVector genotype_log_posteriors) const;
    GroupResponsibilityVector
    init_responsibilities(const GroupOptionalLogPriorArray& group_log_priors,
                          const GroupConcentrationVector& group_concentrations,
                          const MixtureConcentrationArray& mixture_concentrations,
                          const ProbabilityVector& genotype_priors,
                          const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    void
    update_responsibilities(GroupResponsibilityVector& result,
                            const GroupOptionalLogPriorArray& group_log_priors,
                            const GroupConcentrationVector& group_concentrations,
                            const MixtureConcentrationArray& mixture_concentrations,
                            const ProbabilityVector& genotype_posteriors,
                            const ComponentResponsibilityMatrix& component_responsibilities,
                            const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    ComponentResponsibilityMatrix
    init_responsibilities(const GroupConcentrationVector& group_concentrations,
                          const MixtureConcentrationArray& mixture_concentrations,
                          const ProbabilityVector& genotype_priors,
                          const GroupResponsibilityVector& group_responsibilities,
                          const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    void
    update_responsibilities(ComponentResponsibilityMatrix& result,
                            const GroupConcentrationVector& group_concentrations,
                            const MixtureConcentrationArray& mixture_concentrations,
                            const ProbabilityVector& genotype_posteriors,
                            const GroupResponsibilityVector& group_responsibilities,
                            const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    void
    update_genotype_log_posteriors(LogProbabilityVector& result,
                                   const LogProbabilityVector& genotype_log_priors,
                                   const GroupResponsibilityVector& group_responsibilities,
                                   const ComponentResponsibilityMatrix& component_responsibilities,
                                   const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    LogProbability
    marginalise(const GroupResponsibilityVector& group_responsibilities,
                const ComponentResponsibilityMatrix& component_responsibilities,
                const HaplotypeLikelihoodMatrix& log_likelihoods,
                std::size_t g) const noexcept;
    LogProbability
    marginalise(const ComponentResponsibilityVector& responsibilities,
                const HaplotypeLikelihoodVector& log_likelihoods) const noexcept;
    void
    update_group_concentrations(GroupConcentrationVector& result,
                                const GroupConcentrationVector& prior_group_concentrations,
                                const GroupResponsibilityVector& group_responsibilities) const;
    void
    update_mixture_concentrations(MixtureConcentrationArray& result,
                                  const MixtureConcentrationArray& prior_mixture_concentrations,
                                  const GroupResponsibilityVector& group_responsibilities,
                                  const ComponentResponsibilityMatrix& component_responsibilities) const;
    double
    calculate_evidence(const GroupConcentrationVector& prior_group_concentrations,
                       const GroupConcentrationVector& posterior_group_concentrations,
                       const MixtureConcentrationArray& prior_mixture_concentrations,
                       const MixtureConcentrationArray& posterior_mixture_concentrations,
                       const LogProbabilityVector& genotype_log_priors,
                       const LogProbabilityVector& genotype_log_posteriors,
                       const ProbabilityVector& genotype_posteriors,
                       const GroupResponsibilityVector& group_responsibilities,
                       const ComponentResponsibilityMatrix& component_responsibilities,
                       const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    
    // debug
    void print_concentrations(const GroupConcentrationVector& concentrations) const;
    void print_concentrations(const MixtureConcentrationArray& concentrations) const;
    void print(const ProbabilityVector& probabilities) const;
    void print(const GroupResponsibilityVector& responsibilities) const;
    void print(const ComponentResponsibilityMatrix& responsibilities) const;
    void print(const HaplotypeLikelihoodMatrix& log_likelihoods) const;
};

} // namespace model
} // namespace octopus

#endif
