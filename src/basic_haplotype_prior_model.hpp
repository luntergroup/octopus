//
//  basic_haplotype_prior_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef basic_haplotype_prior_model_hpp
#define basic_haplotype_prior_model_hpp

#include "haplotype_prior_model.hpp"

#include <vector>
#include <unordered_map>
#include <functional>

#include "reference_genome.hpp"

class Haplotype;

namespace Octopus
{
    class BasicHaplotypePriorModel : public HaplotypePriorModel
    {
    public:
        BasicHaplotypePriorModel() = delete;
        
        explicit BasicHaplotypePriorModel(const ReferenceGenome& reference);
        
        explicit BasicHaplotypePriorModel(const ReferenceGenome& reference,
                                          double transition_rate, double transversion_rate);
        
        ~BasicHaplotypePriorModel() = default;
        
        BasicHaplotypePriorModel(const BasicHaplotypePriorModel&)            = default;
        BasicHaplotypePriorModel& operator=(const BasicHaplotypePriorModel&) = default;
        BasicHaplotypePriorModel(BasicHaplotypePriorModel&&)                 = default;
        BasicHaplotypePriorModel& operator=(BasicHaplotypePriorModel&&)      = default;
        
    private:
        // ln p(to | from)
        double do_evaluate(const Haplotype& to, const Haplotype& from) const override;
        
        HaplotypePriorMap
        do_compute_maximum_entropy_haplotype_set(std::vector<Haplotype>& haplotypes) const override;
        
    private:
        std::reference_wrapper<const ReferenceGenome> reference_;
        
        const double transition_rate_   = 0.000222;
        const double transversion_rate_ = 0.000111;
    };
    
} // namespace Octopus

#endif /* basic_haplotype_prior_model_hpp */
