// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef single_cell_model_hpp
#define single_cell_model_hpp

#include <vector>
#include <unordered_map>
#include <cstddef>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/phylogeny.hpp"
#include "core/models/haplotype_likelihood_array.hpp"
#include "exceptions/program_error.hpp"
#include "single_cell_prior_model.hpp"
#include "population_prior_model.hpp"
#include "uniform_population_prior_model.hpp"
#include "variational_bayes_mixture_mixture_model.hpp"

namespace octopus { namespace model {

class SingleCellModel
{
public:
    struct Parameters
    {
        using GroupOptionalPriorArray = VariationalBayesMixtureMixtureModel::GroupOptionalPriorArray;
        boost::optional<GroupOptionalPriorArray> group_priors = boost::none;
        std::vector<double> sample_dropout_concentrations = {};
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
        boost::optional<std::size_t> max_genotype_combinations;
        unsigned max_seeds = 20;
        ExecutionPolicy execution_policy = ExecutionPolicy::seq;
    };
    
    struct NoViableGenotypeCombinationsError : public ProgramError
    {
        std::string do_where() const override { return "SingleCellModel"; }
        std::string do_why() const override { return "No viable genotype combinations found"; };
    };
    
    using GenotypeVector = std::vector<Genotype<IndexedHaplotype<>>>;
    using PhylogenyNodePloidyMap = std::unordered_map<SingleCellPriorModel::CellPhylogeny::LabelType, unsigned>;
    
    SingleCellModel() = delete;
    
    SingleCellModel(std::vector<SampleName> samples,
                    SingleCellPriorModel prior_model,
                    Parameters parameters,
                    AlgorithmParameters config,
                    boost::optional<const PopulationPriorModel&> population_prior_model = boost::none);
    
    SingleCellModel(const SingleCellModel&)            = default;
    SingleCellModel& operator=(const SingleCellModel&) = default;
    SingleCellModel(SingleCellModel&&)                 = default;
    SingleCellModel& operator=(SingleCellModel&&)      = default;
    
    ~SingleCellModel() = default;
    
    const SingleCellPriorModel& prior_model() const;
    
    Inferences
    evaluate(const GenotypeVector& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    
    Inferences
    evaluate(const PhylogenyNodePloidyMap& phylogeny_ploidies,
             const GenotypeVector& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;

private:
    const static UniformPopulationPriorModel default_population_prior_model_;
    
    std::vector<SampleName> samples_;
    SingleCellPriorModel prior_model_;
    VariationalBayesMixtureMixtureModel posterior_model_;
    Parameters parameters_;
    AlgorithmParameters config_;
    const PopulationPriorModel* population_prior_model_;
    
    using GenotypeCombination = std::vector<std::size_t>;
    using GenotypeCombinationVector = std::vector<GenotypeCombination>;
    using VBLikelihoodMatrix = VariationalBayesMixtureMixtureModel::HaplotypeLikelihoodMatrix;
    using VBSeedVector = std::vector<VariationalBayesMixtureMixtureModel::LogProbabilityVector>;
    
    std::vector<std::size_t>
    propose_genotypes(const GenotypeVector& genotypes,
                      const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    GenotypeCombinationVector
    propose_genotype_combinations(const GenotypeVector& genotypes,
                                  const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    GenotypeCombinationVector
    propose_all_genotype_combinations(const GenotypeVector& genotypes) const;
    GenotypeCombinationVector
    propose_genotype_combinations(const PhylogenyNodePloidyMap& phylogeny_ploidies,
                                  const GenotypeVector& genotypes,
                                  const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    void
    evaluate(Inferences& result,
             const GenotypeVector& genotypes,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    void
    evaluate(Inferences& result,
             const GenotypeVector& genotypes,
             const GenotypeCombinationVector& genotype_combinations,
             const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    VariationalBayesMixtureMixtureModel::LogProbabilityVector
    calculate_genotype_priors(const GenotypeCombinationVector& genotype_combinations,
                              const GenotypeVector& genotypes) const;
    VBLikelihoodMatrix
    make_likelihood_matrix(const GenotypeCombinationVector& genotype_combinations,
                           const GenotypeVector& genotypes,
                           const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    VBSeedVector
    propose_seeds(const GenotypeCombinationVector& genotype_combinations,
                  const GenotypeVector& genotypes,
                  const VariationalBayesMixtureMixtureModel::LogProbabilityVector& genotype_combination_priors,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods) const;
    VariationalBayesMixtureMixtureModel::Inferences
    evaluate_model(const VariationalBayesMixtureMixtureModel::LogProbabilityVector& genotype_combination_priors,
                   const VBLikelihoodMatrix& haplotype_likelihoods,
                   VBSeedVector seeds) const;
};

} // namespace model
} // namespace octopus

#endif
