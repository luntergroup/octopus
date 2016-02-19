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

class Haplotype;

namespace Octopus
{
    class BasicHaplotypePriorModel : public HaplotypePriorModel
    {
    public:
        BasicHaplotypePriorModel()  = default;
        explicit BasicHaplotypePriorModel(double transition_rate, double transversion_rate);
        ~BasicHaplotypePriorModel() = default;
        
        BasicHaplotypePriorModel(const BasicHaplotypePriorModel&)            = default;
        BasicHaplotypePriorModel& operator=(const BasicHaplotypePriorModel&) = default;
        BasicHaplotypePriorModel(BasicHaplotypePriorModel&&)                 = default;
        BasicHaplotypePriorModel& operator=(BasicHaplotypePriorModel&&)      = default;
        
    private:
        using HaplotypePriorModel::HaplotypePriorMap;
        
        // ln p(to | from)
        double do_evaluate(const Haplotype& to, const Haplotype& from) const override;
        
        HaplotypePriorMap do_evaluate(std::vector<Haplotype>::const_iterator first,
                                      std::vector<Haplotype>::const_iterator last,
                                      std::vector<Haplotype>::const_iterator reference) const override;
        
    private:
        const double transition_rate_   = 0.000222;
        const double transversion_rate_ = 0.000111;
    };
    
} // namespace Octopus

#endif /* basic_haplotype_prior_model_hpp */
