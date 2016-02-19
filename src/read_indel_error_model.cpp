//
//  read_indel_error_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "read_indel_error_model.hpp"

#include "tandem.hpp"

namespace Octopus
{
namespace
{
    auto extract_repeats(const Haplotype& haplotype)
    {
        return Tandem::find_maximal_repetitions(haplotype.get_sequence(), 1, 2);
    }
} // namespace

std::vector<std::uint8_t> ReadIndelErrorModel::calculate_gap_open_penalties(const Haplotype& haplotype) const
{
//    const auto repeats = extract_repeats(haplotype);
//    
//    if (repeats.empty()) {
//        
//    }
    
    std::vector<std::uint8_t> result(sequence_size(haplotype), 30);
    
    return result;
}
} // namespace Octopus
