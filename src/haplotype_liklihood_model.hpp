//
//  haplotypeliklihood_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef haplotype_liklihood_model_hpp
#define haplotype_liklihood_model_hpp

#include "common.hpp"
#include "haplotype.hpp"
#include "read_indel_error_model.hpp"
#include "kmer_mapper.hpp"

class AlignedRead;

namespace Octopus
{
    class HaplotypeLikelihoodModel
    {
    public:
        enum class InactiveRegionState { Clear, Unclear };
        
        HaplotypeLikelihoodModel()  = delete;
        HaplotypeLikelihoodModel(KmerMapper mapper);
        ~HaplotypeLikelihoodModel() = default;
        
        HaplotypeLikelihoodModel(const HaplotypeLikelihoodModel&)            = default;
        HaplotypeLikelihoodModel& operator=(const HaplotypeLikelihoodModel&) = default;
        HaplotypeLikelihoodModel(HaplotypeLikelihoodModel&&)                 = default;
        HaplotypeLikelihoodModel& operator=(HaplotypeLikelihoodModel&&)      = default;
        
        // ln p(read | haplotype)
        double log_probability(const AlignedRead& read, const Haplotype& haplotype,
                               InactiveRegionState flank_state = InactiveRegionState::Clear) const;
        
    private:
        ReadIndelErrorModel indel_error_model_;
        KmerMapper mapper_;
    };
} // namespace Octopus

#endif /* haplotype_liklihood_model_hpp */
