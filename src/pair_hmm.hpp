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

struct Model
{
    int nucprior, gapextend;
};

// p(target | truth, target_qualities, target_gap_open_penalties, model)
double compute_log_conditional_probability(const std::string& truth, const std::string& target,
                                           const std::vector<std::uint8_t>& target_qualities,
                                           const std::vector<std::uint8_t>& target_gap_open_penalties,
                                           std::size_t target_offset_into_truth_hint, const Model& model);

#endif /* pair_hmm_hpp */
