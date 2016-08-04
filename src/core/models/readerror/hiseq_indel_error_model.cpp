// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "hiseq_indel_error_model.hpp"

#include <algorithm>
#include <iterator>

#include <tandem/tandem.hpp>

#include <core/types/haplotype.hpp>

#include <iostream> // DEBUG

namespace octopus {
    
    constexpr decltype(HiSeqIndelErrorModel::Homopolymer_errors_) HiSeqIndelErrorModel::Homopolymer_errors_;
    constexpr decltype(HiSeqIndelErrorModel::Homopolymer_errors_) HiSeqIndelErrorModel::Di_nucleotide_tandem_repeat_errors_;
    constexpr decltype(HiSeqIndelErrorModel::Homopolymer_errors_) HiSeqIndelErrorModel::Tri_nucleotide_tandem_repeat_errors_;
    constexpr decltype(HiSeqIndelErrorModel::Homopolymer_errors_) HiSeqIndelErrorModel::Poly_nucleotide_tandem_repeat_errors_;
    constexpr decltype(HiSeqIndelErrorModel::default_gap_extension_) HiSeqIndelErrorModel::default_gap_extension_;
    
    namespace
    {
        auto extract_repeats(const Haplotype& haplotype)
        {
            return tandem::find_maximal_repetitions(haplotype.sequence(), 1, 3);
        }
    }
    
    template <typename C, typename T>
    static auto get_penalty(const C& penalties, const T length)
    {
        return (length < penalties.size()) ? penalties[length - 1] : penalties.back();
    }
    
    HiSeqIndelErrorModel::PenaltyType
    HiSeqIndelErrorModel::evaluate(const Haplotype& haplotype, PenaltyVector& gap_open_penalities) const
    {
        using std::begin; using std::end; using std::cbegin; using std::cend; using std::next;
        
        const auto repeats = extract_repeats(haplotype);
        
        gap_open_penalities.assign(sequence_size(haplotype), 50);
        
        tandem::StringRun max_repeat {};
        
        for (const auto& repeat : repeats) {
            std::int8_t e;
            
            if (repeat.period == 1) {
                e = get_penalty(Homopolymer_errors_, repeat.length);
            } else if (repeat.period == 2) {
                static constexpr std::array<char, 2> AC {'A', 'C'};
                
                const auto it = next(cbegin(haplotype.sequence()), repeat.pos);
                
                e = get_penalty(Di_nucleotide_tandem_repeat_errors_, repeat.length);
                
                if (e > 10 && std::equal(cbegin(AC), cend(AC), it)) {
                    e -= 2;
                }
            } else if (repeat.period == 3) {
                static constexpr std::array<char, 3> GGC {'G', 'G', 'C'};
                static constexpr std::array<char, 3> GCC {'G', 'C', 'C'};
                
                const auto it = next(cbegin(haplotype.sequence()), repeat.pos);
                
                e = get_penalty(Tri_nucleotide_tandem_repeat_errors_, repeat.length);
                
                if (e > 10 && std::equal(cbegin(GGC), cend(GGC), it)) {
                    e -= 3;
                } else if (e > 12 && std::equal(cbegin(GCC), cend(GCC), it)) {
                    e -= 2;
                }
            } else {
                e = get_penalty(Poly_nucleotide_tandem_repeat_errors_, repeat.length);
                ++e;
            }
            
            std::fill_n(next(begin(gap_open_penalities), repeat.pos), repeat.length, e);
            
            if (repeat.length > max_repeat.length) {
                max_repeat = repeat;
            }
        }
        
        switch (max_repeat.period) {
            case 1: return 3;
            case 2: return 2;
            case 3: return 1;
            default: return default_gap_extension_;
        }
    }
} // namespace octopus
