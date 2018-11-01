// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef allelic_dropout_genotype_likelihood_model_hpp
#define allelic_dropout_genotype_likelihood_model_hpp

#include <vector>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/cancer_genotype.hpp"
#include "core/models/haplotype_likelihood_array.hpp"
#include "constant_mixture_genotype_likelihood_model.hpp"
#include "variable_mixture_genotype_likelihood_model.hpp"

namespace octopus { namespace model {

class AllelicDropoutGenotypeLikelihoodModel
{
public:
    using LogProbability = double;
    
    struct Parameters
    {
        double dropout_rate;
    };
    
    AllelicDropoutGenotypeLikelihoodModel() = delete;
    
    AllelicDropoutGenotypeLikelihoodModel(Parameters params, const HaplotypeLikelihoodArray& likelihoods);
    AllelicDropoutGenotypeLikelihoodModel(Parameters params, const HaplotypeLikelihoodArray& likelihoods,
                                          const std::vector<Haplotype>& haplotypes);
    
    AllelicDropoutGenotypeLikelihoodModel(const AllelicDropoutGenotypeLikelihoodModel&)            = default;
    AllelicDropoutGenotypeLikelihoodModel& operator=(const AllelicDropoutGenotypeLikelihoodModel&) = default;
    AllelicDropoutGenotypeLikelihoodModel(AllelicDropoutGenotypeLikelihoodModel&&)                 = default;
    AllelicDropoutGenotypeLikelihoodModel& operator=(AllelicDropoutGenotypeLikelihoodModel&&)      = default;
    
    ~AllelicDropoutGenotypeLikelihoodModel() = default;
    
    Parameters parameters() const;
    const HaplotypeLikelihoodArray& cache() const noexcept;
    
    void prime(const std::vector<Haplotype>& haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    LogProbability evaluate(const Genotype<Haplotype>& genotype) const;
    LogProbability evaluate(const GenotypeIndex& genotype) const;
    LogProbability evaluate(const CancerGenotype<Haplotype>& genotype) const;
    LogProbability evaluate(const CancerGenotypeIndex& genotype) const;
    
private:
    Parameters params_;
    ConstantMixtureGenotypeLikelihoodModel constant_mixture_model_;
    mutable VariableMixtureGenotypeLikelihoodModel variable_mixture_model_;
    
    LogProbability dropout_probability() const noexcept;
};

} // namespace model
} // namespace octopus

#endif
