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

namespace PairHMM
{
    struct Model
    {
        int nucprior;
        int gapextend;
        bool do_lhs_backtrace = false;
        bool do_rhs_backtrace = false;
    };
    
    // p(target | truth, target_qualities, target_gap_open_penalties, model)
    double align_around_offset(const std::string& truth, const std::string& target,
                               const std::vector<std::uint8_t>& target_qualities,
                               const std::vector<std::int8_t>& truth_gap_open_penalties,
                               std::size_t target_offset, const Model& model);
} // namespace PairHMM

#endif /* pair_hmm_hpp */
