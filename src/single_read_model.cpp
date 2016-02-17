//
//  single_read_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "single_read_model.hpp"

#include <cmath>
#include <algorithm>

#include "mappable.hpp"
#include "pair_hmm.hpp"

#include <iostream> // TEST
#include <chrono>   // TEST

namespace Octopus
{
    // public methods
    
    double SingleReadModel::log_probability(const AlignedRead& read, const Haplotype& haplotype)
    {
        Model model {2, 3};
        
        std::vector<std::uint8_t> gap_open_penalities(sequence_size(haplotype), 20);
        
        const auto offset_hint = begins_before(haplotype, read) ? begin_distance(read, haplotype) : 0;
        
        return compute_log_conditional_probability(haplotype.get_sequence(), read.get_sequence(),
                                                   read.get_qualities(), gap_open_penalities,
                                                   offset_hint, model);
    }
} // namespace Octopus
