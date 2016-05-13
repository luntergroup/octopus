//
//  read_indel_error_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "read_indel_error_model.hpp"

#include <algorithm>
#include <iterator>

#include "tandem.hpp"

#include <iostream> // DEBUG

namespace Octopus
{

constexpr decltype(ReadIndelErrorModel::Homopolymer_errors_) ReadIndelErrorModel::Homopolymer_errors_;
constexpr decltype(ReadIndelErrorModel::gap_extension_) ReadIndelErrorModel::gap_extension_;

namespace
{
    auto extract_repeats(const Haplotype& haplotype)
    {
        return Tandem::find_maximal_repetitions(haplotype.sequence(), 1, 3);
    }
}

ReadIndelErrorModel::PenaltyType
ReadIndelErrorModel::calculate_gap_extension_penalty(const Haplotype& haplotype) const noexcept
{
    return gap_extension_;
}

template <typename C, typename T>
static auto get_penalty(const C& penalties, const T length)
{
    return (length < penalties.size()) ? penalties[length - 1] : penalties.back();
}

void ReadIndelErrorModel::fill_gap_open_penalties(const Haplotype& haplotype,
                                                  std::vector<PenaltyType>& result) const
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
    
    const auto repeats = extract_repeats(haplotype);
    
    result.resize(sequence_size(haplotype), Homopolymer_errors_.front());
    
    for (const auto& repeat : repeats) {
        if (repeat.period == 1) {
            std::fill_n(next(begin(result), repeat.pos), repeat.length,
                        get_penalty(Homopolymer_errors_, repeat.length));
        } else if (repeat.period == 2) {
            static constexpr std::array<char, 2> AC {'A', 'C'};
            
            const auto it = next(cbegin(haplotype.sequence()), repeat.pos);
            
            std::int8_t e {30};
            
            if (std::equal(cbegin(AC), cend(AC), it)) {
                e = 20;
            }
            
            std::fill_n(next(begin(result), repeat.pos), repeat.length, e);
        }else if (repeat.period == 3) {
            static constexpr std::array<char, 3> GGC {'G', 'G', 'C'};
            static constexpr std::array<char, 3> GCC {'G', 'C', 'C'};
            
            const auto it = next(cbegin(haplotype.sequence()), repeat.pos);
            
            std::int8_t e {38};
            
            if (std::equal(cbegin(GGC), cend(GGC), it)) {
                e = 15;
            } else if (std::equal(cbegin(GCC), cend(GCC), it)) {
                e = 20;
            }
            
            std::fill_n(next(begin(result), repeat.pos), repeat.length, e);
        }
    }
}
} // namespace Octopus
