//
//  coalescent_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef coalescent_model_hpp
#define coalescent_model_hpp

#include <vector>

#include "haplotype.hpp"

namespace Octopus
{
    class CoalescentModel
    {
    public:
        CoalescentModel() = default;
        
        explicit CoalescentModel(Haplotype reference_haplotype,
                                 double snp_heterozygosity = 0.001,
                                 double indel_heterozygosity = 0.001);
        explicit CoalescentModel(std::vector<Haplotype> reference_haplotypes,
                                 double snp_heterozygosity = 0.001,
                                 double indel_heterozygosity = 0.001);
        
        ~CoalescentModel() = default;
        
        CoalescentModel(const CoalescentModel&)            = default;
        CoalescentModel& operator=(const CoalescentModel&) = default;
        CoalescentModel(CoalescentModel&&)                 = default;
        CoalescentModel& operator=(CoalescentModel&&)      = default;
        
        double evaluate(const std::vector<Haplotype>& haplotypes) const;
        
    private:
        std::vector<Haplotype> reference_haplotypes_;
        
        double snp_heterozygosity_, indel_heterozygosity_;
    };
} // namespace Octopus

#endif /* coalescent_model_hpp */
