//
//  somatic_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef somatic_model_hpp
#define somatic_model_hpp

#include <type_traits>
#include <functional>

#include "coalescent_model.hpp"
#include "cancer_genotype.hpp"
#include "haplotype.hpp"
#include "maths.hpp"

namespace Octopus
{
    class SomaticModel
    {
    public:
        SomaticModel() = default;
        explicit SomaticModel(const CoalescentModel& germline_model, double somatic_mutation_rate = 0.00001);
        ~SomaticModel() = default;
        
        SomaticModel(const SomaticModel&)            = default;
        SomaticModel& operator=(const SomaticModel&) = default;
        SomaticModel(SomaticModel&&)                 = default;
        SomaticModel& operator=(SomaticModel&&)      = default;
        
        const CoalescentModel& get_germline_model() const noexcept;
        
        double evaluate(const CancerGenotype<Haplotype>& genotype) const;
    private:
        std::reference_wrapper<const CoalescentModel> germline_model_;
        double somatic_mutation_rate_;
    };
    
    template <typename Container>
    std::vector<double> calculate_log_priors(const Container& genotypes, const SomaticModel& model)
    {
        static_assert(std::is_same<typename Container::value_type, CancerGenotype<Haplotype>>::value,
                      "genotypes must contain CancerGenotype<Haplotype>'s");
        
        std::vector<double> result(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [&model] (const auto& genotype) {
                           return model.evaluate(genotype);
                       });
        
        Maths::normalise_logs(result);
        
        return result;
    }
} // namespace Octopus

#endif /* somatic_model_hpp */
