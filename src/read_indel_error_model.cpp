//
//  read_indel_error_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "read_indel_error_model.hpp"

#include <array>
#include <algorithm>
#include <iterator>

#include "tandem.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
namespace
{
    auto extract_repeats(const Haplotype& haplotype)
    {
        return Tandem::find_maximal_repetitions(haplotype.get_sequence(), 1, 1);
    }
}

std::vector<std::int8_t> ReadIndelErrorModel::calculate_gap_open_penalties(const Haplotype& haplotype) const
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
    using std::fill_n; using std::equal;
    
    static constexpr std::array<std::int8_t, 50> homopolymer_errors
    {
        45,42,41,39,37,32,28,23,20,19,17,16,15,14,13,12,11,11,10,
        9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
    };
    
    const auto repeats = extract_repeats(haplotype);
    
    std::vector<std::int8_t> result(sequence_size(haplotype), homopolymer_errors.front());
    
    for (const auto& repeat : repeats) {
        if (repeat.period == 1) {
            const auto e = (repeat.length < homopolymer_errors.size()) ? homopolymer_errors[repeat.length - 1] : homopolymer_errors.back();
            fill_n(next(begin(result), repeat.pos), repeat.length, e);
        } else if (repeat.period == 3) {
            static constexpr std::array<char, 3> GGC {'G', 'G', 'C'};
            static constexpr std::array<char, 3> GCC {'G', 'C', 'C'};
            
            const auto it = next(cbegin(haplotype.get_sequence()), repeat.pos);
            
            std::int8_t e {38};
            
            if (equal(cbegin(GGC), cend(GGC), it)) {
                e = 15;
            } else if (equal(cbegin(GCC), cend(GCC), it)) {
                e = 20;
            }
            
            fill_n(next(begin(result), repeat.pos), repeat.length, e);
        }
    }
    
    return result;
}
} // namespace Octopus
