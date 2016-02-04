//
//  haplotype_prior_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_haplotype_prior_model_hpp
#define Octopus_haplotype_prior_model_hpp

#include <vector>
#include <unordered_map>

class Haplotype;

namespace Octopus
{

class HaplotypePriorModel
{
public:
    using HaplotypePriorMap = std::unordered_map<Haplotype, double>;
    
    HaplotypePriorModel()  = default;
    explicit HaplotypePriorModel(double transition_rate, double transversion_rate);
    ~HaplotypePriorModel() = default;
    
    HaplotypePriorModel(const HaplotypePriorModel&)            = default;
    HaplotypePriorModel& operator=(const HaplotypePriorModel&) = default;
    HaplotypePriorModel(HaplotypePriorModel&&)                 = default;
    HaplotypePriorModel& operator=(HaplotypePriorModel&&)      = default;
    
    // ln p(to | from)
    double evaluate(const Haplotype& to, const Haplotype& from) const;
    
    HaplotypePriorMap evaluate(const std::vector<Haplotype>& haplotypes, const Haplotype& reference) const;
    
private:
    const double transition_rate_   = 0.000222;
    const double transversion_rate_ = 0.000111;
};

void make_unique(std::vector<Haplotype>& haplotypes, const HaplotypePriorModel& prior_model);

} // namespace Octopus

#endif
