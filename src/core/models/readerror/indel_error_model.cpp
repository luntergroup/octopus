// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "indel_error_model.hpp"

#include <algorithm>
#include <iterator>

#include <tandem/tandem.hpp>

#include <core/types/haplotype.hpp>

#include <iostream> // DEBUG

namespace octopus {

constexpr decltype(IndelErrorModel::homopolymerErrors_) IndelErrorModel::homopolymerErrors_;
constexpr decltype(IndelErrorModel::homopolymerErrors_) IndelErrorModel::diNucleotideTandemRepeatErrors_;
constexpr decltype(IndelErrorModel::homopolymerErrors_) IndelErrorModel::triNucleotideTandemRepeatErrors_;
constexpr decltype(IndelErrorModel::homopolymerErrors_) IndelErrorModel::polyNucleotideTandemRepeatErrors_;
constexpr decltype(IndelErrorModel::defaultGapExtension_) IndelErrorModel::defaultGapExtension_;

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
    
    gap_open_penalities.assign(sequence_size(haplotype), homopolymerErrors_.front());
    
    tandem::StringRun max_repeat {};
    
    for (const auto& repeat : repeats) {
        std::int8_t e;
        
        switch (repeat.period) {
            case 1:
            {
                e = get_penalty(homopolymerErrors_, repeat.length);
                break;
            }
            case 2:
            {
                static constexpr std::array<char, 2> AC {'A', 'C'};
                
                e = get_penalty(diNucleotideTandemRepeatErrors_, repeat.length / 2);
                
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
                
                e = get_penalty(triNucleotideTandemRepeatErrors_, repeat.length / 3);
                
                const auto it = next(cbegin(haplotype.sequence()), repeat.pos);
                
                if (e > 10 && std::equal(cbegin(GGC), cend(GGC), it)) {
                    e -= 2;
                } else if (e > 12 && std::equal(cbegin(GCC), cend(GCC), it)) {
                    e -= 1;
                }
                break;
            }
            default:
                e = get_penalty(polyNucleotideTandemRepeatErrors_, repeat.length / repeat.period);
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
        default: return defaultGapExtension_;
    }
}

} // namespace octopus
