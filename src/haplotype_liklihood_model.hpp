//
//  haplotypeliklihood_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef haplotype_liklihood_model_hpp
#define haplotype_liklihood_model_hpp

#include <unordered_map>
#include <functional>
#include <cstddef>
#include <unordered_map>

#include "common.hpp"
#include "haplotype.hpp"
#include "read_indel_error_model.hpp"

class AlignedRead;

namespace Octopus
{
    class HaplotypeLikelihoodModel
    {
    public:
        HaplotypeLikelihoodModel()  = default;
        ~HaplotypeLikelihoodModel() = default;
        
        HaplotypeLikelihoodModel(const HaplotypeLikelihoodModel&)            = default;
        HaplotypeLikelihoodModel& operator=(const HaplotypeLikelihoodModel&) = default;
        HaplotypeLikelihoodModel(HaplotypeLikelihoodModel&&)                 = default;
        HaplotypeLikelihoodModel& operator=(HaplotypeLikelihoodModel&&)      = default;
        
        // ln p(read | haplotype)
        double log_probability(const AlignedRead& read, const Haplotype& haplotype);
        
    private:
        ReadIndelErrorModel indel_error_model_;
        
        using CacheKeyType       = Haplotype::SequenceType;
        using GapOpenPenalities  = std::vector<short>;
        
        mutable std::unordered_map<CacheKeyType, GapOpenPenalities> cached_gap_penalities_;
    };
} // namespace Octopus

#endif /* haplotype_liklihood_model_hpp */
