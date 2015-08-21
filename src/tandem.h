//
//  tandem.h
//  tandem
//
//  Created by Daniel Cooke on 17/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __tandem__tandem__
#define __tandem__tandem__

#include <vector>
#include <set>
#include <deque>
#include <cstddef>
#include <algorithm> // std::sort, std::mismatch, std::max, std::min, std::transform, std::lower_bound
#include <iterator>

#include "divsufsort.h"

using std::size_t;

struct StringRun
{
    StringRun() = default;
    StringRun(size_t pos, size_t length, size_t period) : pos {pos}, length {length}, period {period} {}
    size_t pos, length, period;
};

inline bool operator==(const StringRun& lhs, const StringRun& rhs)
{
    return lhs.pos == rhs.pos && lhs.length == rhs.length && lhs.period == rhs.period;
}
inline bool operator!=(const StringRun& lhs, const StringRun& rhs) { return !operator==(lhs, rhs); }

// Wrapper for divsufsort
template <typename T>
std::vector<size_t> make_suffix_array(const T& str)
{
    std::vector<int32_t> sa(str.size());
    
    divsufsort(str.c_str(), &sa[0], static_cast<int>(str.size()));
    
    return std::vector<size_t> {sa.cbegin(), sa.cend()};
}

// Naive implementation - O(nlogn)
//template <typename T>
//std::vector<size_t> make_suffix_array(const T& str)
//{
//    std::vector<std::pair<T, size_t>> suffixes {};
//    suffixes.reserve(str.size());
//    
//    for (size_t i {}; i < str.size(); ++i) {
//        suffixes.emplace_back(str.substr(i), i);
//    }
//    
//    std::sort(suffixes.begin(), suffixes.end());
//    
//    std::vector<size_t> result {};
//    result.reserve(str.size());
//    
//    for (auto& suffix : suffixes) {
//        result.push_back(suffix.second);
//    }
//    
//    return result;
//}

// Naive implementation, but can be fast in practice according to Fischer (2011)
template <typename T>
std::vector<size_t> make_longest_common_prefix_array(const T& str, const std::vector<size_t>& suffix_array)
{
    std::vector<size_t> result(str.size());
    
    result[0] = 0;
    
    for (size_t i {1}; i < suffix_array.size(); ++i) {
        auto it = str.cbegin() + suffix_array[i - 1];
        result[i] = std::mismatch(it, str.end(), str.cbegin() + suffix_array[i]).first - it;
    }
    
    return result;
}

std::vector<size_t> make_longest_common_factor_array(std::vector<size_t> sa, std::vector<size_t> lcp);
std::pair<std::vector<size_t>, std::vector<size_t>> make_lpf_and_prev_occ(std::vector<size_t> sa, std::vector<size_t> lcp);

struct LZBlock
{
    LZBlock() = default;
    LZBlock(size_t pos, size_t length) : pos {pos}, length {length} {}
    size_t pos, length;
};

// Implementation of algorithm found in Crochemore et al. (2008)
template <typename T>
std::vector<LZBlock> lempel_ziv_factorisation(const T& str)
{
    auto sa  = make_suffix_array(str);
    auto lcp = make_longest_common_prefix_array(str, sa);
    auto lcf = make_longest_common_factor_array(std::move(sa), std::move(lcp));
    
    std::vector<LZBlock> result {};
    result.reserve(str.size() - 1); // max possible blocks
    
    size_t end {1};
    
    result.emplace_back(0, end);
    
    while (end < str.size()) {
        auto m = std::max(size_t {1}, lcf[end]);
        result.emplace_back(end, m);
        end += m;
    }
    
    result.shrink_to_fit();
    
    return result;
}

// Implementation of algorithm found in Crochemore et al. (2008)
template <typename T>
std::pair<std::vector<LZBlock>, std::vector<size_t>> lempel_ziv_factorisation_with_backpointers(const T& str)
{
    auto sa  = make_suffix_array(str);
    auto lcp = make_longest_common_prefix_array(str, sa);
    
    std::vector<size_t> lpf, prev_occ;
    std::tie(lpf, prev_occ) = make_lpf_and_prev_occ(std::move(sa), std::move(lcp));
    
    std::vector<LZBlock> lz_blocks {};
    lz_blocks.reserve(str.size() - 1); // max possible blocks
    std::vector<size_t> backpointers {};
    backpointers.reserve(str.size() - 1);
    
    size_t end {1};
    
    lz_blocks.emplace_back(0, end);
    backpointers.emplace_back(-1);
    
    while (end < str.size()) {
        auto m = std::max(size_t {1}, lpf[end]);
        lz_blocks.emplace_back(end, m);
        backpointers.emplace_back(prev_occ[end]);
        end += m;
    }
    
    lz_blocks.shrink_to_fit();
    backpointers.shrink_to_fit();
    
    return {lz_blocks, backpointers};
}

// Llie et al (2010) show that this direct computation of LCE actually outperforms other O(n) and O(log n)
// methods in practice

template <typename T>
size_t forward_lce(const T& str, size_t i, size_t j, size_t n)
{
    return std::mismatch(str.cbegin() + i, str.cbegin() + n, str.cbegin() + j).first - (str.cbegin() + i);
}

template <typename T>
size_t forward_lce(const T& str, size_t i, size_t j)
{
    return forward_lce(str, i, j, str.size());
}

template <typename T>
size_t backward_lce(const T& str, size_t i, size_t j, size_t n)
{
    return std::mismatch(str.crend() - i - 1, str.crend() - n, str.crend() - j - 1).first - (str.crend() - i - 1);
}

template <typename T>
size_t backward_lce(const T& str, size_t i, size_t j)
{
    return backward_lce(str, i, j, 0);
}

// Implements Mains algorithm found in Main (1989)
template <typename T>
std::vector<StringRun> find_leftmost_maximal_repetitions(const T& str, const std::vector<LZBlock>& lz_blocks,
                                                         size_t min_period = 1)
{
    std::vector<StringRun> result {};
    
    for (size_t h {1}; h < lz_blocks.size(); ++h) {
        auto u = lz_blocks[h].pos;
        auto n = lz_blocks[h].length;
        auto m = std::min(u, 2 * lz_blocks[h - 1].length + n);
        auto t = u - m;
        
        // rightmax periodicities
        for (size_t j {min_period}; j <= n; ++j) {
            auto ls = backward_lce(str, u - 1, u + j - 1, t);
            auto lp = forward_lce(str, u + j, u, u + n);
            
            if (ls + lp >= j) {
                result.emplace_back(u - ls, j + lp + ls, j);
            }
        }
        
        // leftmax periodicities
        for (size_t j {min_period}; j < m; ++j) {
            auto ls = backward_lce(str, u - j - 1, u - 1, t);
            auto lp = forward_lce(str, u, u - j, u + n);
            
            if (ls + lp >= j) {
                result.emplace_back(u - (ls + j), j + lp + ls, j);
            }
        }
    }
    
    return result;
}

template <typename T>
std::vector<StringRun> find_leftmost_maximal_repetitions(const T& str, size_t min_period = 1)
{
    return find_leftmost_maximal_repetitions(str, lempel_ziv_factorisation(str), min_period);
}

// Implements the algorithm described in Kolpakov & Kucherov (1999)
template <typename T>
std::vector<StringRun> find_maximal_repetitions(const T& str, size_t min_period = 1)
{
    std::vector<LZBlock> lz_blocks;
    std::vector<size_t> prev_lz_block_occurrence;
    std::tie(lz_blocks, prev_lz_block_occurrence) = lempel_ziv_factorisation_with_backpointers(str);
    
    // first bucket sort by end position
    
    std::vector<std::vector<StringRun>> end_buckets(str.size(), std::vector<StringRun> {});
    
    for (const auto& run : find_leftmost_maximal_repetitions(str, lz_blocks, min_period)) {
        auto end = run.pos + run.length - 1;
        if (std::find(end_buckets[end].cbegin(), end_buckets[end].cend(), run) == end_buckets[end].cend()) {
            end_buckets[end].push_back(run);
        }
    }
    
    // second bucket sort by start position
    
    std::vector<std::deque<StringRun>> sorted_buckets(str.size(), std::deque<StringRun> {});
    
    size_t num_runs {};
    
    for (auto& bucket : end_buckets) {
        for (const auto& run : bucket) {
            sorted_buckets[run.pos].push_back(run);
            ++num_runs;
        }
        bucket.clear();
        bucket.shrink_to_fit();
    }
    
    end_buckets.clear();
    end_buckets.shrink_to_fit();
    
    for (size_t k {}; k < lz_blocks.size(); ++k) {
        const auto& block = lz_blocks[k];
        auto block_end    = block.pos + block.length;
        
        auto delta = block.pos - ((prev_lz_block_occurrence[k] == -1) ? 0 : prev_lz_block_occurrence[k]);
        auto v_end = block_end - delta;
        
        for (auto j = block.pos; j < block_end; ++j) {
            auto v = j - delta;
            
            auto last_v = std::lower_bound(sorted_buckets[v].cbegin(), sorted_buckets[v].cend(), v_end,
                                           [] (const auto& run, auto val) {
                                               return run.pos + run.length < val;
                                           });
            
            std::transform(sorted_buckets[v].cbegin(), last_v, std::front_inserter(sorted_buckets[j]),
                           [delta] (const auto& run) {
                               return StringRun {run.pos + delta, run.length, run.period};
                           });
            
            num_runs += std::distance(sorted_buckets[v].cbegin(), last_v);
        }
    }
    
    lz_blocks.clear();
    lz_blocks.shrink_to_fit();
    prev_lz_block_occurrence.clear();
    prev_lz_block_occurrence.shrink_to_fit();
    
    std::vector<StringRun> result {};
    result.reserve(num_runs);
    
    for (auto& bucket : sorted_buckets) {
        result.insert(result.end(), bucket.cbegin(), bucket.cend());
        bucket.clear();
        bucket.shrink_to_fit();
    }
    
    return result;
}

#endif /* defined(__tandem__tandem__) */
