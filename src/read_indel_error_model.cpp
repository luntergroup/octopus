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
        return Tandem::find_maximal_repetitions(haplotype.get_sequence(), 1, 1);
    }
}

ReadIndelErrorModel::PenaltyType
ReadIndelErrorModel::calculate_gap_extension_penalty(const Haplotype& haplotype) const noexcept
{
    return gap_extension_;
}

void ReadIndelErrorModel::fill_gap_open_penalties(const Haplotype& haplotype,
                                                  std::vector<PenaltyType>& result) const
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
    using std::fill_n; using std::equal;
    
    const auto repeats = extract_repeats(haplotype);
    
    result.resize(sequence_size(haplotype), Homopolymer_errors_.front());
    
    for (const auto& repeat : repeats) {
        if (repeat.period == 1) {
            const auto e = (repeat.length < Homopolymer_errors_.size())
                ? Homopolymer_errors_[repeat.length - 1] : Homopolymer_errors_.back();
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
}
} // namespace Octopus
