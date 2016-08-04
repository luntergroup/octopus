//
//  indel_error_model.cpp
//  octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "indel_error_model.hpp"

#include <algorithm>
#include <iterator>

#include <tandem/tandem.hpp>

#include <core/types/haplotype.hpp>

#include <iostream> // DEBUG

namespace octopus {

constexpr decltype(IndelErrorModel::Homopolymer_errors_) IndelErrorModel::Homopolymer_errors_;
constexpr decltype(IndelErrorModel::Homopolymer_errors_) IndelErrorModel::Di_nucleotide_tandem_repeat_errors_;
constexpr decltype(IndelErrorModel::Homopolymer_errors_) IndelErrorModel::Tri_nucleotide_tandem_repeat_errors_;
constexpr decltype(IndelErrorModel::Homopolymer_errors_) IndelErrorModel::Poly_nucleotide_tandem_repeat_errors_;
constexpr decltype(IndelErrorModel::default_gap_extension_) IndelErrorModel::default_gap_extension_;

namespace {
    auto extract_repeats(const Haplotype& haplotype)
    {
        return tandem::find_maximal_repetitions(haplotype.sequence(), 1, 3);
    }
}

template <typename C, typename T>
static auto get_penalty(const C& penalties, const T length)
{
    return (length < penalties.size()) ? penalties[length] : penalties.back();
}

IndelErrorModel::PenaltyType
IndelErrorModel::evaluate(const Haplotype& haplotype, PenaltyVector& gap_open_penalities) const
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
    
    const auto repeats = extract_repeats(haplotype);
    
    gap_open_penalities.assign(sequence_size(haplotype), Homopolymer_errors_.front());
    
    tandem::StringRun max_repeat {};
    
    for (const auto& repeat : repeats) {
        std::int8_t e;
        
        switch (repeat.period) {
            case 1:
            {
                e = get_penalty(Homopolymer_errors_, repeat.length);
                break;
            }
            case 2:
            {
                static constexpr std::array<char, 2> AC {'A', 'C'};
                
                e = get_penalty(Di_nucleotide_tandem_repeat_errors_, repeat.length / 2);
                
                const auto it = next(cbegin(haplotype.sequence()), repeat.pos);
                
                if (e > 10 && std::equal(cbegin(AC), cend(AC), it)) {
                    e -= 2;
                }
                break;
            }
            case 3:
            {
                static constexpr std::array<char, 3> GGC {'G', 'G', 'C'};
                static constexpr std::array<char, 3> GCC {'G', 'C', 'C'};
                
                e = get_penalty(Tri_nucleotide_tandem_repeat_errors_, repeat.length / 3);
                
                const auto it = next(cbegin(haplotype.sequence()), repeat.pos);
                
                if (e > 10 && std::equal(cbegin(GGC), cend(GGC), it)) {
                    e -= 2;
                } else if (e > 12 && std::equal(cbegin(GCC), cend(GCC), it)) {
                    e -= 1;
                }
                break;
            }
            default:
                e = get_penalty(Poly_nucleotide_tandem_repeat_errors_, repeat.length / repeat.period);
        }
        
        std::fill_n(next(begin(gap_open_penalities), repeat.pos), repeat.length, e);
        
        if (repeat.length > max_repeat.length) {
            max_repeat = repeat;
        }
    }
    
    switch (max_repeat.period) {
        case 1: return 3;
        case 2: return 5;
        case 3: return 5;
        default: return default_gap_extension_;
    }
}

} // namespace octopus
