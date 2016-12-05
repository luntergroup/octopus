// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef somatic_mutation_model_hpp
#define somatic_mutation_model_hpp

#include <type_traits>
#include <functional>

#include "coalescent_model.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/cancer_genotype.hpp"
#include "utils/maths.hpp"

namespace octopus {

class SomaticMutationModel
{
public:
    struct Parameters
    {
        double somatic_mutation_rate = 0.00001;
    };
    
    SomaticMutationModel() = delete;
    
    SomaticMutationModel(Parameters params);
    
    SomaticMutationModel(const SomaticMutationModel&)            = default;
    SomaticMutationModel& operator=(const SomaticMutationModel&) = default;
    SomaticMutationModel(SomaticMutationModel&&)                 = default;
    SomaticMutationModel& operator=(SomaticMutationModel&&)      = default;
    
    ~SomaticMutationModel() = default;
    
    // ln p(somatic | germline)
    double evaluate(const Haplotype& somatic, const Haplotype& germline) const;

private:
    Parameters params_;
};

} // namespace octopus

#endif
