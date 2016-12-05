// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cancer_genotype_prior_model_hpp
#define cancer_genotype_prior_model_hpp

#include <type_traits>
#include <functional>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/cancer_genotype.hpp"
#include "genotype_prior_model.hpp"
#include "../mutation/somatic_mutation_model.hpp"
#include "utils/maths.hpp"

namespace octopus {

class CancerGenotypePriorModel
{
public:
    CancerGenotypePriorModel() = delete;
    
    CancerGenotypePriorModel(const GenotypePriorModel& germline_model,
                             SomaticMutationModel mutation_model);
    
    CancerGenotypePriorModel(const CancerGenotypePriorModel&)            = default;
    CancerGenotypePriorModel& operator=(const CancerGenotypePriorModel&) = default;
    CancerGenotypePriorModel(CancerGenotypePriorModel&&)                 = default;
    CancerGenotypePriorModel& operator=(CancerGenotypePriorModel&&)      = default;
    
    ~CancerGenotypePriorModel() = default;
    
    double evaluate(const CancerGenotype<Haplotype>& genotype) const;

private:
    std::reference_wrapper<const GenotypePriorModel> germline_model_;
    SomaticMutationModel mutation_model_;
    
    // p(somatic | germline)
    double ln_probability_of_somatic(const Haplotype& somatic, const Genotype<Haplotype>& germline) const;
    double ln_probability_of_somatic(const Haplotype& somatic, const Haplotype& germline) const;
};

template <typename Container>
std::vector<double> calculate_log_priors(const Container& genotypes,
                                         const CancerGenotypePriorModel& model)
{
    static_assert(std::is_same<typename Container::value_type, CancerGenotype<Haplotype>>::value,
                  "genotypes must contain CancerGenotype<Haplotype>'s");
    std::vector<double> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&model](const auto& genotype) {
                       return model.evaluate(genotype);
                   });
    maths::normalise_logs(result);
    return result;
}
    
} // namespace octopus

#endif
