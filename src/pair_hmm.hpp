//
//  pair_hmm.hpp
//  pair_hmm
//
//  Created by Daniel Cooke on 14/12/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef pair_hmm_hpp
#define pair_hmm_hpp

#include <string>
#include <vector>
#include <cstddef>
#include <functional>

namespace PairHMM
{
    static constexpr unsigned AlignmenetPad {15};
    
    struct Model
    {
        using PenaltyType = std::int8_t;
        
        std::reference_wrapper<const std::vector<PenaltyType>> snv_priors;
        std::reference_wrapper<const std::vector<PenaltyType>> gap_open_penalties;
        
        short gap_extend;
        
        short nuc_prior = 2;
        
        std::size_t lhs_flank_size = 0;
        std::size_t rhs_flank_size = 0;
    };
    
    // p(target | truth, target_qualities, model)
    double align(const std::string& truth, const std::string& target,
                 const std::vector<std::uint8_t>& target_qualities,
                 std::size_t target_offset, const Model& model);
} // namespace PairHMM

#endif /* pair_hmm_hpp */
