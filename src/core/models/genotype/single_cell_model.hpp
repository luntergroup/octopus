// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef single_cell_model_hpp
#define single_cell_model_hpp

#include <vector>
#include <cstddef>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/phylogeny.hpp"
#include "core/models/haplotype_likelihood_array.hpp"
#include "single_cell_prior_model.hpp"
#include "variational_bayes_mixture_mixture_model.hpp"

namespace octopus { namespace model {

class SingleCellModel
{
public:
    struct Parameters
    {
        double dropout_concentration = 1.5;
        double group_concentration = 1.0;
    };
    
    struct Inferences
    {
        struct GroupInferences
        {
            std::vector<double> genotype_posteriors;
            std::vector<double> sample_attachment_posteriors;
        };
        using InferedPhylogeny = Phylogeny<std::size_t, GroupInferences>;
        InferedPhylogeny phylogeny;
        double log_evidence;
    };
    
    struct AlgorithmParameters
    {
        std::size_t max_genotype_combinations;
        unsigned max_seeds = 5;
    };
    
    SingleCellModel() = delete;
    
    SingleCellModel(std::vector<SampleName> samples,
                    SingleCellPriorModel prior_model,
                    Parameters parameters,
                    AlgorithmParameters config);
    
    SingleCellModel(const SingleCellModel&)            = default;
    SingleCellModel& operator=(const SingleCellModel&) = default;
    SingleCellModel(SingleCellModel&&)                 = default;
    SingleCellModel& operator=(SingleCellModel&&)      = default;
    
    ~SingleCellModel() = default;
    
    Inferences
    evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    Inferences
    evaluate(const std::vector<GenotypeIndex>& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;

private:
    std::vector<SampleName> samples_;
    SingleCellPriorModel prior_model_;
    VariationalBayesMixtureMixtureModel posterior_model_;
    Parameters parameters_;
    AlgorithmParameters config_;
    
    using GenotypeCombination = std::vector<std::size_t>;
    using GenotypeCombinationVector = std::vector<GenotypeCombination>;
    using VBLikelihoodMatrix = VariationalBayesMixtureMixtureModel::HaplotypeLikelihoodMatrix;
    using VBSeedVector = std::vector<VariationalBayesMixtureMixtureModel::LogProbabilityVector>;
    
    GenotypeCombinationVector
    propose_genotype_combinations(const std::vector<Genotype<Haplotype>>& genotypes,
                                  const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    GenotypeCombinationVector
    propose_all_genotype_combinations(const std::vector<Genotype<Haplotype>>& genotypes) const;
    VariationalBayesMixtureMixtureModel::LogProbabilityVector
    calculate_genotype_priors(const GenotypeCombinationVector& genotype_combinations,
                              const std::vector<Genotype<Haplotype>>& genotypes) const;
    VBLikelihoodMatrix
    make_likelihood_matrix(const GenotypeCombinationVector& genotype_combinations,
                           const std::vector<Genotype<Haplotype>>& genotypes,
                           const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    VBSeedVector
    propose_seeds(const GenotypeCombinationVector& genotype_combinations,
                  const std::vector<Genotype<Haplotype>>& genotypes,
                  const VariationalBayesMixtureMixtureModel::LogProbabilityVector& genotype_combination_priors,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
};

} // namespace model
} // namespace octopus

#endif
