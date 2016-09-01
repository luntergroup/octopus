// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <vector>
#include <deque>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <tuple>

#include "tandem/tandem.hpp"

inline
std::vector<unsigned>
calculate_positional_repeat_coverage(const std::vector<Tandem::StringRun>& repeats,
                                     const std::size_t sequence_length,
                                     const std::size_t min_repeat_length = 4)
{
    std::vector<unsigned> result(sequence_length, 0);
    
    for (const auto& repeat : repeats) {
        if (repeat.length > min_repeat_length) {
            const auto it = std::next(std::begin(result), repeat.pos);
            std::transform(it, std::next(it, repeat.length), it, [] (const auto c) { return c + 1; });
        }
    }
    
    return result;
}

template <typename T>
std::vector<unsigned>
calculate_positional_repeat_coverage(const T& sequence,
                                     const std::size_t min_repeat_length = 4)
{
    return calculate_positional_repeat_coverage(Tandem::find_maximal_repetitions(sequence), sequence.size(), min_repeat_length);
}

template <typename T>
std::size_t num_repeat_positions(const T& sequence, const std::size_t min_repeat_length = 4)
{
    const auto positions = calculate_positional_repeat_coverage(sequence, min_repeat_length);
    
    return std::count_if(std::cbegin(positions), std::cend(positions),
                         [] (const auto c) { return c > 0; });
}

template <typename T>
std::vector<unsigned char>
mask_sequence(const T& sequence,
              const std::size_t min_mask_size = 20,
              const std::size_t max_join_size = 50,
              const std::size_t window_size = 10,
              const std::size_t min_repeat_length = 4)
{
    std::vector<unsigned char> result(sequence.size(), 0);
    
    if (sequence.size() < window_size) {
        return result;
    }
    
    auto positions = calculate_positional_repeat_coverage(sequence, min_repeat_length);
    
    auto it  = std::cbegin(positions);
    auto it2 = std::next(it, window_size);
    
    auto result_it = std::begin(result);
    
    auto z = static_cast<double>(std::accumulate(it, std::prev(it2), 0));
    
    while (it2 != std::cend(positions)) {
        z += *it2;
        if (z / window_size >= 0.75) {
            std::fill_n(result_it, window_size, 1);
        }
        z -= *it;
        ++result_it;
        ++it;
        ++it2;
    }
    
    result_it = std::find(std::begin(result), std::end(result), 1);
    
    while (result_it != std::end(result)) {
        auto it = std::find(std::next(result_it), std::end(result), 0);
        
        if (std::distance(result_it, it) <= min_mask_size) {
            std::fill(result_it, it, 0);
        }
        
        result_it = it;
    }
    
    result_it = std::find(std::begin(result), std::end(result), 1);
    
    while (result_it != std::end(result)) {
        auto it  = std::find(std::next(result_it), std::end(result), 0);
        
        if (it == std::end(result)) break;
        
        auto it2 = std::find(std::next(result_it), std::end(result), 1);
        
        if (std::distance(it, it2) <= max_join_size) {
            std::fill(it, it2, 1);
        }
        
        result_it = it2;
    }
    
    return result;
}

struct ComplexBlock
{
    ComplexBlock() = default;
    ComplexBlock(std::size_t pos, std::size_t length) : pos {pos}, length {length} {}
    std::size_t pos;
    std::size_t length;
};

template <typename T>
std::deque<ComplexBlock>
find_high_complexity_subsequences(const T& sequence,
                                  const std::size_t min_mask_size = 50,
                                  const std::size_t max_join_size = 60,
                                  const std::size_t window_size = 10,
                                  const std::size_t min_repeat_length = 4)
{
    const auto mask = mask_sequence(sequence, min_mask_size, max_join_size, window_size, min_repeat_length);
    
    auto it = std::find(std::cbegin(mask), std::cend(mask), 1);
    
    std::deque<ComplexBlock> result {};
    
    while (it != std::cend(mask)) {
        const auto it2 = std::find(std::next(it), std::cend(mask), 0);
        
        result.emplace_back(std::distance(std::cbegin(mask), it), std::distance(it, it2));
        
        it = std::find(it2, std::cend(mask), 1);
    }
    
    return result;
}
