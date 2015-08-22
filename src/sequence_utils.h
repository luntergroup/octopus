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

std::vector<GenomicRegion> find_exact_tandem_repeats(ReferenceGenome& reference, const GenomicRegion& region,
                                                     GenomicRegion::SizeType min_repeat_size = 2)
{
    auto maximal_repetitions = find_maximal_repetitions(reference.get_sequence(region), min_repeat_size);
    
    std::vector<GenomicRegion> result {};
    result.reserve(maximal_repetitions.size());
    
    auto offset = region.get_begin();
    
    for (auto& run : maximal_repetitions) {
        result.emplace_back(region.get_contig_name(), run.pos + offset, run.pos + run.length + offset);
    }
    
    return result;
}

template <typename SequenceType>
std::map<ContigRegion, SequenceType> find_tandem_repeats(const SequenceType& sequence, unsigned min_repeat_size = 3)
{
    return false;
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
