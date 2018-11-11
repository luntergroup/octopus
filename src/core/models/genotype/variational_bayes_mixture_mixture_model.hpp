// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variational_bayes_mixture_mixture_model_hpp
#define variational_bayes_mixture_mixture_model_hpp

#include <vector>

#include "core/models/haplotype_likelihood_array.hpp"
#include "variational_bayes_mixture_model.hpp"

namespace octopus { namespace model {

class VariationalBayesMixtureMixtureModel
{
public:
    struct Options
    {
        double epsilon = 0.05;
        unsigned max_iterations = 1000;
        double save_memory = false;
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
    using GroupResponsabilityVector = std::vector<Sigma>; // One element per sample
    
    struct Inferences
    {
        ProbabilityVector genotype_posteriors;
        LogProbabilityVector genotype_log_posteriors;
        GroupConcentrationVector group_concentrations;
        MixtureConcentrationArray mixture_concentrations;
        GroupResponsabilityVector group_responsabilities;
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
             const GroupConcentrationVector& group_concentrations,
             const MixtureConcentrationArray& mixture_concentrations,
             std::vector<LogProbabilityVector> seeds) const;
    
    Inferences
    evaluate(const LogProbabilityVector& genotype_log_priors,
             const HaplotypeLikelihoodMatrix& log_likelihoods,
             double group_concentration, double mixture_concentration,
             std::vector<LogProbabilityVector> seeds) const;
    
private:
    using Tau = std::vector<double>; // One element per read
    using ComponentResponsabilityVector = std::vector<Tau>; // One element per haplotype in genotype (max)
    using ComponentResponsabilityMatrix = std::vector<ComponentResponsabilityVector>; // One element per sample
    
    Options options_;
    
    Inferences
    evaluate(const LogProbabilityVector& genotype_log_priors,
             const HaplotypeLikelihoodMatrix& log_likelihoods,
             const GroupConcentrationVector& group_concentrations,
             const MixtureConcentrationArray& mixture_concentrations,
             LogProbabilityVector& genotype_log_posteriors) const;
    GroupResponsabilityVector
    init_responsabilities(const GroupConcentrationVector& group_concentrations,
                          const MixtureConcentrationArray& mixture_concentrations,
                          const ProbabilityVector& genotype_priors,
                          const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    void
    update_responsabilities(GroupResponsabilityVector& result,
                            const GroupConcentrationVector& group_concentrations,
                            const MixtureConcentrationArray& mixture_concentrations,
                            const ProbabilityVector& genotype_posteriors,
                            const ComponentResponsabilityMatrix& component_responsabilities,
                            const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    ComponentResponsabilityMatrix
    init_responsabilities(const GroupConcentrationVector& group_concentrations,
                          const MixtureConcentrationArray& mixture_concentrations,
                          const ProbabilityVector& genotype_priors,
                          const GroupResponsabilityVector& group_responsabilities,
                          const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    void
    update_responsabilities(ComponentResponsabilityMatrix& result,
                            const GroupConcentrationVector& group_concentrations,
                            const MixtureConcentrationArray& mixture_concentrations,
                            const ProbabilityVector& genotype_posteriors,
                            const GroupResponsabilityVector& group_responsabilities,
                            const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    void
    update_genotype_log_posteriors(LogProbabilityVector& result,
                                   const LogProbabilityVector& genotype_log_priors,
                                   const GroupResponsabilityVector& group_responsabilities,
                                   const ComponentResponsabilityMatrix& component_responsabilities,
                                   const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    LogProbability
    marginalise(const GroupResponsabilityVector& group_responsabilities,
                const ComponentResponsabilityMatrix& component_responsabilities,
                const HaplotypeLikelihoodMatrix& log_likelihoods,
                std::size_t g) const noexcept;
    LogProbability
    marginalise(const ComponentResponsabilityVector& responsabilities,
                const HaplotypeLikelihoodVector& log_likelihoods) const noexcept;
    void
    update_group_concentrations(GroupConcentrationVector& result,
                                const GroupConcentrationVector& prior_group_concentrations,
                                const GroupResponsabilityVector& group_responsabilities) const;
    void
    update_mixture_concentrations(MixtureConcentrationArray& result,
                                  const MixtureConcentrationArray& prior_mixture_concentrations,
                                  const GroupResponsabilityVector& group_responsabilities,
                                  const ComponentResponsabilityMatrix& component_responsabilities) const;
    double
    calculate_evidence(const GroupConcentrationVector& prior_group_concentrations,
                       const GroupConcentrationVector& posterior_group_concentrations,
                       const MixtureConcentrationArray& prior_mixture_concentrations,
                       const MixtureConcentrationArray& posterior_mixture_concentrations,
                       const LogProbabilityVector& genotype_log_priors,
                       const LogProbabilityVector& genotype_log_posteriors,
                       const ProbabilityVector& genotype_posteriors,
                       const GroupResponsabilityVector& group_responsabilities,
                       const ComponentResponsabilityMatrix& component_responsabilities,
                       const HaplotypeLikelihoodMatrix& log_likelihoods) const;
    
    // debug
    void print_concentrations(const GroupConcentrationVector& concentrations) const;
    void print_concentrations(const MixtureConcentrationArray& concentrations) const;
    void print(const ProbabilityVector& probabilities) const;
    void print(const GroupResponsabilityVector& responsabilities) const;
    void print(const ComponentResponsabilityMatrix& responsabilities) const;
    void print(const HaplotypeLikelihoodMatrix& log_likelihoods) const;
};

} // namespace model
} // namespace octopus

#endif
