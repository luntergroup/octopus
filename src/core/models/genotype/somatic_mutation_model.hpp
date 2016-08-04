// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef somatic_mutation_model_hpp
#define somatic_mutation_model_hpp

#include <type_traits>
#include <functional>

#include <core/models/genotype/coalescent_model.hpp>
#include <core/types/cancer_genotype.hpp>
#include <core/types/haplotype.hpp>
#include <utils/maths.hpp>

namespace octopus {

class SomaticMutationModel
{
public:
    SomaticMutationModel() = delete;
    
    SomaticMutationModel(const CoalescentModel& germline_model,
                         double somatic_mutation_rate = 0.00001);
    
    SomaticMutationModel(const SomaticMutationModel&)            = default;
    SomaticMutationModel& operator=(const SomaticMutationModel&) = default;
    SomaticMutationModel(SomaticMutationModel&&)                 = default;
    SomaticMutationModel& operator=(SomaticMutationModel&&)      = default;
    
    ~SomaticMutationModel() = default;
    
    double evaluate(const CancerGenotype<Haplotype>& genotype) const;
    
private:
    std::reference_wrapper<const CoalescentModel> germline_model_;
    double somatic_mutation_rate_;
};

template <typename Container>
std::vector<double> calculate_log_priors(const Container& genotypes,
                                         const SomaticMutationModel& model)
{
    static_assert(std::is_same<typename Container::value_type, CancerGenotype<Haplotype>>::value,
                  "genotypes must contain CancerGenotype<Haplotype>'s");
    
    std::vector<double> result(genotypes.size());
    
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&model] (const auto& genotype) {
                       return model.evaluate(genotype);
                   });
    
    maths::normalise_logs(result);
    
    return result;
}

} // namespace octopus

#endif /* somatic_mutation_model_hpp */
