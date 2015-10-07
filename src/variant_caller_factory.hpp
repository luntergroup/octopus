//
//  variant_caller_factory.hpp
//  Octopus
//
//  Created by Daniel Cooke on 03/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef variant_caller_factory_hpp
#define variant_caller_factory_hpp

#include <unordered_map>
#include <string>
#include <memory> // std::unique_ptr
#include <stdexcept>

#include "variant_caller.hpp"
#include "population_caller.hpp"
#include "cancer_caller.hpp"

namespace Octopus {

inline std::unique_ptr<VariantCaller>
make_variant_caller(const std::string& model, ReferenceGenome& reference,
                    CandidateVariantGenerator& candidate_generator,
                    VariantCaller::RefCallType refcall_type,
                    double min_variant_posterior, double min_refcall_posterior,
                    unsigned ploidy, const SampleIdType& normal_sample, double min_somatic_posterior)
{
    std::unordered_map<std::string, std::function<std::unique_ptr<VariantCaller>()>> model_map {
        {"population", [&] () {
            return std::make_unique<PopulationVariantCaller>(reference, candidate_generator,
                                                             refcall_type,
                                                             min_variant_posterior, min_refcall_posterior,
                                                             ploidy);
        }},
        {"cancer",     [&] () {
            return std::make_unique<CancerVariantCaller>(reference, candidate_generator,
                                                         refcall_type, min_variant_posterior, min_somatic_posterior, 
                                                         min_refcall_posterior, normal_sample);
        }}
        };
    
    if (model_map.count(model) == 0) {
        throw std::runtime_error {"unknown model " + model};
    }
    
    return model_map[model]();
}

} // namespace Octopus

#endif /* variant_caller_factory_hpp */
