//
//  sequence_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 14/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__sequence_utils__
#define __Octopus__sequence_utils__

#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <cstddef>
#include <iterator>
#include <algorithm> // std::count_if

#include <iostream> // TEST

#include "contig_region.h"
#include "tandem.h"

template <typename SequenceType>
bool is_n_free(const SequenceType& sequence)
{
    return std::none_of(std::cbegin(sequence), std::cend(sequence), 'N');
}

template <typename SequenceType>
bool is_dna(const SequenceType& sequence)
{
    return sequence.find_first_not_of("ACGTN") == SequenceType::npos;
}

template <typename SequenceType>
bool is_rna(const SequenceType& sequence)
{
    return sequence.find_first_not_of("ACGUN") == SequenceType::npos;
}

template <typename SequenceType>
bool is_dna_rna_ambiguous(const SequenceType& sequence)
{
    return sequence.find_first_of("TU") != SequenceType::npos;
}

namespace detail {
    static constexpr std::array<char, 128> rc_table {
        4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4,  4,
        4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4,  4,
        4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4,  4,
        4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4,  4,
        4, 84, 4, 71, 4,  4,  4, 67, 4, 4, 4, 4,  4, 4, 78, 4,
        4, 4,  4, 4,  65, 65, 4, 4,  4, 4, 4, 4,  4, 4, 4,  4,
        4, 84, 4, 71, 4,  4,  4, 67, 4, 4, 4, 4,  4, 4, 4,  4,
        4, 4,  4, 4,  65, 65, 4, 4,  4, 4, 4, 4,  4, 4, 4,  4
    };
}

constexpr char complement(char base)
{
    return detail::rc_table[base];
}

template <typename SequenceType>
SequenceType reverse_complement(SequenceType sequence)
{
    auto f_itr = std::begin(sequence);
    auto r_itr = std::prev(std::end(sequence));
    
    for (; f_itr < r_itr; ++f_itr, --r_itr) {
        char c {complement(*f_itr)};
        *f_itr = complement(*r_itr);
        *r_itr = c;
    }
    
    if (f_itr == r_itr) {
        *f_itr = complement(*f_itr); // complement middle base if sequence is odd length
    }
    
    return sequence;
}

template <typename SequenceType>
bool is_palindromic(const SequenceType& sequence)
{
    if (sequence.size() % 2 != 0) return false;
    
    auto f_itr = std::cbegin(sequence);
    auto r_itr = std::prev(std::cend(sequence));
    
    for (; f_itr < r_itr; ++f_itr, --r_itr) {
        if (*f_itr != complement(*r_itr)) return false;
    }
    
    return true;
}

template <typename SequenceType>
std::unordered_map<char, std::size_t> count_bases(const SequenceType& sequence)
{
    std::unordered_map<char, std::size_t> result {};
    result.reserve(5); // 4 bases + N
    
    for (char base : sequence) {
        ++result[base];
    }
    
    return result;
}

// Replaces all contiguous sub-sequences of N's with a single N, inplace, and returns a map of
// each N position in the new sequence, and how many N's have been removed up to the first non-N
// base past the position
template <typename SequenceType>
std::map<std::size_t, std::size_t> collapse_ns(SequenceType& sequence)
{
    std::map<std::size_t, std::size_t> result {};
    
    auto first = std::begin(sequence);
    auto last  = std::end(sequence);
    
    std::size_t position {}, num_removed {};
    
    while (first != last) {
        auto it1 = std::adjacent_find(first, last, [] (char lhs, char rhs) { return lhs == 'N' && lhs == rhs; });
        
        if (it1 == last) break;
        
        auto it2 = std::find_if_not(it1, last, [] (char base) { return base == 'N'; });
        
        position    += std::distance(first, it1);
        num_removed += std::distance(it1, it2) - 1;
        
        result.emplace(position, num_removed);
        
        first = it2;
    }
    
    sequence.erase(std::unique(std::begin(sequence) + result.cbegin()->first, last,
                               [] (char lhs, char rhs) { return lhs == 'N' && lhs == rhs; }), last);
    
    return result;
}

namespace detail
{
    void rebase(std::vector<StringRun>& runs, const std::map<std::size_t, std::size_t>& shift_map)
    {
        auto shift_map_it = std::cbegin(shift_map);
        for (auto& run : runs) {
            while (std::next(shift_map_it)->first <= run.pos) ++shift_map_it;
            run.pos += static_cast<decltype(run.pos)>(shift_map_it->second);
        }
    }
}

struct TandemRepeat
{
    using SizeType = GenomicRegion::SizeType;
    TandemRepeat() = delete;
    template <typename T>
    TandemRepeat(T region, GenomicRegion::SizeType period) : region {std::forward<T>(region)}, period {period} {}
    
    GenomicRegion region;
    GenomicRegion::SizeType period;
};

template <typename SequenceType>
std::vector<TandemRepeat> find_exact_tandem_repeats(SequenceType sequence, const GenomicRegion& region,
                                                    GenomicRegion::SizeType min_repeat_size = 2,
                                                    GenomicRegion::SizeType max_repeat_size = 10000)
{
    auto n_shift_map = collapse_ns(sequence);
    
    auto maximal_repetitions = find_maximal_repetitions(sequence , min_repeat_size, max_repeat_size);
    
    detail::rebase(maximal_repetitions, n_shift_map);
    
    std::vector<TandemRepeat> result {};
    result.reserve(maximal_repetitions.size());
    
    auto offset = region.get_begin();
    
    for (const auto& run : maximal_repetitions) {
        result.emplace_back(GenomicRegion {region.get_contig_name(),
            static_cast<GenomicRegion::SizeType>(run.pos + offset),
            static_cast<GenomicRegion::SizeType>(run.pos + run.length + offset)
        }, run.period);
    }
    
    return result;
}

template <typename SequenceType>
double gc_bias(const SequenceType& sequence)
{
    auto gc_count = std::count_if(std::cbegin(sequence), std::cend(sequence),
                                  [] (char base) { return base == 'G' || base == 'C'; });
    return static_cast<double>(gc_count) / sequence.length();
}

inline constexpr unsigned base_hash(char base)
{
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'U':
            return 3;
        default:
            return 0; // should never do this
    }
}

template <typename SequenceType>
std::vector<ContigRegion> find_cpg_islands(const SequenceType& sequence)
{
    std::vector<ContigRegion> result {};
    
//    // transitions probabilities taken from Durbin et. al (1998)
//    constexpr std::array<std::array<std::array<std::array<double, 4>, 4>, 2>, 2> transitions
//    {{
//        {{{{0.180, 0.274, 0.426, 0.120},
//           {0.171, 0.368, 0.274, 0.188},
//           {0.161, 0.355, 0.384, 0.125},
//           {0.079, 0.355, 0.384, 0.182}}},
//          {{0.180, 0.274, 0.426, 0.120},
//           {0.171, 0.368, 0.274, 0.188},
//           {0.161, 0.355, 0.384, 0.125},
//           {0.079, 0.355, 0.384, 0.182}}},
//        {{{{0.300, 0.205, 0.285, 0.210},
//           {0.322, 0.298, 0.078, 0.302},
//           {0.248, 0.246, 0.298, 0.208},
//           {0.177, 0.239, 0.292, 0.292}}},
//          {{0.300, 0.205, 0.285, 0.210},
//           {0.322, 0.298, 0.078, 0.302},
//           {0.248, 0.246, 0.298, 0.208},
//           {0.177, 0.239, 0.292, 0.292}}}
//    }}}};
//    
//    using State = std::array<std::array<unsigned, 4>, 2>;
//    
//    std::vector<State> states(sequence.length());
//    auto state = states.begin();
//    
//    for (std::size_t i {}; i < sequence.length(); ++i) {
//        for (unsigned a {}; a < 4; ++a) {
//            states[i][0][a] = (base_hash(sequence[i]) == a) * std::max({
//                transitions[0][0][a] * states[i - 1][0][0],
//                transitions[1][0][a] * states[i - 1][1][0],
//                transitions[0][1][a] * states[i - 1][0][1],
//                transitions[1][1][a] * states[i - 1][1][1],
//                transitions[0][2][a] * states[i - 1][0][2],
//                transitions[1][2][a] * states[i - 1][1][2],
//                transitions[0][3][a] * states[i - 1][0][3],
//                transitions[1][3][a] * states[i - 1][1][3],
//            });
//            
//            states[i][0][a] = (base_hash(sequence[i]) == a) * std::max({
//                transitions[0][0][a] * states[i - 1][0][0],
//                transitions[1][0][a] * states[i - 1][1][0],
//                transitions[0][1][a] * states[i - 1][0][1],
//                transitions[1][1][a] * states[i - 1][1][1],
//                transitions[0][2][a] * states[i - 1][0][2],
//                transitions[1][2][a] * states[i - 1][1][2],
//                transitions[0][3][a] * states[i - 1][0][3],
//                transitions[1][3][a] * states[i - 1][1][3],
//            });
//        }
//    }
    
    return result;
}

template <typename SequenceType>
double sequence_complexity(const SequenceType& sequence)
{
    return 0;
}

template <typename SequenceType>
std::vector<ContigRegion> find_low_complex_regions(const SequenceType& sequence)
{
    std::vector<ContigRegion> result {};
    
    return result;
}

#endif /* defined(__Octopus__sequence_utils__) */
