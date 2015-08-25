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
#include <deque>
#include <cstdint>   // uint32_t
#include <cstddef>   // size_t
#include <algorithm> // std::mismatch, std::max, std::min, std::find, std::transform, std::lower_bound
#include <iterator>  // std::distance, std::cbegin, std::cend, std::begin, std::end, std::crend
#include <tuple>     // std::pair, std::tie
#include <utility>   // std::move
#include <numeric>   // std::accumulate

#include "divsufsort.h"

/**
 This is a small library for finding exact tandem repeats in sequences. It also provides a
 simple API for making common string processing data structures, that are used in the repeat
 finding algorithm. The main method is:
 
 template <typename T>
 std::vector<StringRun> find_maximal_repetitions(const T& str, uint32_t min_period, uint32_t max_period)
 
 where T can be any contiguous sequence container with value_type char (e.g. std::string, std::vector<char>)
 that has a size() member function (so no const char*).
 
 The library only works with uint32_t at the moment. This is mostly because divsufsort does not use templates,
 so I'd either have to use uint32_t or uint64_t (or size_t).
 */

namespace Tandem
{
    
    struct StringRun
    {
        StringRun() = default;
        explicit StringRun(uint32_t pos, uint32_t length, uint32_t period) : pos {pos}, length {length}, period {period} {}
        
        uint32_t pos, length, period;
    };
    
    inline bool operator==(const StringRun& lhs, const StringRun& rhs)
    {
        return lhs.pos == rhs.pos && lhs.length == rhs.length;
    }
    inline bool operator!=(const StringRun& lhs, const StringRun& rhs) { return !operator==(lhs, rhs); }
    
    // Wrapper for divsufsort
    template <typename T>
    std::vector<uint32_t> make_suffix_array(const T& str)
    {
        std::vector<saidx_t> sa(str.size()); // divsufsort requires signed integers so need to copy
        
        divsufsort(str.data(), &sa[0], static_cast<int>(str.size()));
        
        return std::vector<uint32_t> {std::cbegin(sa), std::cend(sa)};
    }
    
    // rank array is inverse suffix array
    inline std::vector<uint32_t> make_rank_array(const std::vector<uint32_t>& suffix_array)
    {
        std::vector<uint32_t> result(suffix_array.size());
        
        for (uint32_t i {}; i < suffix_array.size(); ++i) {
            result[suffix_array[i]] = i;
        }
        
        return result;
    }
    
    namespace detail
    {
        // Llie et al (2010) show that this direct computation of LCE actually outperforms other O(n) and O(log n)
        // methods in practice.
        
        template <typename T>
        uint32_t forward_lce(const T& str, uint32_t i, uint32_t j, uint32_t n)
        {
            using std::cbegin;
            return static_cast<uint32_t>(std::distance(cbegin(str) + i, std::mismatch(cbegin(str) + i,
                                                                                      cbegin(str) + n, cbegin(str) + j).first));
        }
        
        template <typename T>
        uint32_t forward_lce(const T& str, uint32_t i, uint32_t j)
        {
            return forward_lce(str, i, j, static_cast<uint32_t>(str.size()));
        }
        
        template <typename T>
        uint32_t backward_lce(const T& str, uint32_t i, uint32_t j, uint32_t n)
        {
            using std::crend;
            return static_cast<uint32_t>(std::distance(crend(str) - i - 1, std::mismatch(crend(str) - i - 1,
                                                                                         crend(str) - n, crend(str) - j - 1).first));
        }
        
        template <typename T>
        uint32_t backward_lce(const T& str, uint32_t i, uint32_t j)
        {
            return backward_lce(str, i, j, uint32_t {});
        }
    } // end namespace detail
    
    // LCP = Longest Common Prefix. O(n) implementation given in Kasai et al (2001).
    template <typename T>
    std::vector<uint32_t> make_lcp_array(const T& str, const std::vector<uint32_t>& suffix_array)
    {
        auto rank = make_rank_array(suffix_array);
        
        std::vector<uint32_t> result(suffix_array.size());
        
        for (uint32_t i {}, h {}; i < suffix_array.size(); ++i) {
            if (rank[i] > 0) {
                h += detail::forward_lce(str, i + h, suffix_array[rank[i] - 1] + h);
                result[rank[i]] = h;
                if (h > 0) --h;
            }
        }
        
        return result;
    }
    
    // LPF = Longest Previous Factor
    std::vector<uint32_t> make_lpf_array(std::vector<uint32_t> sa, std::vector<uint32_t> lcp);
    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> make_lpf_and_prev_occ_arrays(std::vector<uint32_t> sa, std::vector<uint32_t> lcp);
    
    struct LZBlock
    {
        LZBlock() = default;
        explicit LZBlock(uint32_t pos, uint32_t length) : pos {pos}, length {length} {}
        uint32_t pos, length;
    };
    
    // Implementation of algorithm found in Crochemore et al. (2008)
    template <typename T>
    std::vector<LZBlock> lempel_ziv_factorisation(const T& str)
    {
        if (str.empty()) return {};
        
        auto sa  = make_suffix_array(str);
        auto lcp = make_lcp_array(str, sa);
        auto lpf = make_lpf_array(std::move(sa), std::move(lcp));
        
        std::vector<LZBlock> result {};
        result.reserve(str.size()); // max possible blocks
        
        uint32_t end {1};
        result.emplace_back(0, end);
        
        while (end < str.size()) {
            auto m = std::max(uint32_t {1}, lpf[end]);
            result.emplace_back(end, m);
            end += m;
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    // Implementation of algorithm found in Crochemore et al. (2008)
    template <typename T>
    std::pair<std::vector<LZBlock>, std::vector<uint32_t>> lempel_ziv_factorisation_with_prev_block_occurences(const T& str)
    {
        if (str.empty()) return {{}, {}};
        
        auto sa  = make_suffix_array(str);
        auto lcp = make_lcp_array(str, sa);
        
        std::vector<uint32_t> lpf, prev_occ;
        std::tie(lpf, prev_occ) = make_lpf_and_prev_occ_arrays(std::move(sa), std::move(lcp));
        
        std::vector<LZBlock> lz_blocks {};
        lz_blocks.reserve(str.size()); // max possible blocks
        std::vector<uint32_t> prev_lz_block_occurrence {};
        prev_lz_block_occurrence.reserve(str.size());
        
        uint32_t end {1};
        
        lz_blocks.emplace_back(0, end);
        prev_lz_block_occurrence.emplace_back(-1);
        
        while (end < str.size()) {
            auto m = std::max(uint32_t {1}, lpf[end]);
            lz_blocks.emplace_back(end, m);
            prev_lz_block_occurrence.emplace_back(prev_occ[end]);
            end += m;
        }
        
        lz_blocks.shrink_to_fit();
        prev_lz_block_occurrence.shrink_to_fit();
        
        return {lz_blocks, prev_lz_block_occurrence};
    }
    
    namespace detail
    {
        // Implements Mains algorithm found in Main (1989).
        template <typename T>
        std::deque<StringRun> find_leftmost_maximal_repetitions(const T& str, const std::vector<LZBlock>& lz_blocks,
                                                                const uint32_t min_period = 1, const uint32_t max_period = -1)
        {
            std::deque<StringRun> result {}; // using std::deque rather than std::vector for nice growth properties
            
            for (uint32_t h {1}; h < lz_blocks.size(); ++h) {
                const auto u   = lz_blocks[h].pos;
                const auto n   = lz_blocks[h].length;
                const auto m   = std::min(u, 2 * lz_blocks[h - 1].length + n);
                const auto t   = u - m;
                const auto end = u + n;
                
                // rightmax periodicities
                for (uint32_t j {min_period}; j <= std::min(n, max_period); ++j) {
                    auto ls = detail::backward_lce(str, u - 1, u + j - 1, t);
                    auto lp = detail::forward_lce(str, u + j, u, end);
                    
                    if (ls > 0 && ls + lp >= j && j + lp < n) {
                        result.emplace_back(u - ls, j + lp + ls, j);
                    }
                }
                
                // leftmax periodicities
                for (uint32_t j {min_period}; j < std::min(m, max_period); ++j) {
                    auto ls = detail::backward_lce(str, u - j - 1, u - 1, t);
                    auto lp = detail::forward_lce(str, u, u - j, end);
                    
                    if (ls + lp >= j) {
                        result.emplace_back(u - (ls + j), j + lp + ls, j);
                    }
                }
            }
            
            return result;
        }
        
        // just reserves enough space to avoid reallocations
        std::vector<std::vector<StringRun>> get_init_buckets(const size_t n, const std::deque<StringRun>& lmrs);
        
        template <typename T>
        std::vector<std::vector<StringRun>>
        get_end_buckets(const T& str, const std::vector<LZBlock>& lz_blocks, const uint32_t min_period, const uint32_t max_period)
        {
            auto lmrs = find_leftmost_maximal_repetitions(str, lz_blocks, min_period, max_period);
            
            auto result = get_init_buckets(str.size(), lmrs);
            
            for (const auto& run : lmrs) {
                auto& curr = result[run.pos + run.length - 1];
                if (std::find(std::cbegin(curr), std::cend(curr), run) == std::cend(curr)) {
                    curr.push_back(run);
                }
            }
            
            return result;
        }
        
        // just reserves enough space to avoid reallocations
        std::vector<std::vector<StringRun>> get_init_buckets(const size_t n, const std::vector<std::vector<StringRun>>& end_buckets);
        
        template <typename T>
        std::vector<std::vector<StringRun>>
        get_sorted_buckets(const T& str, const std::vector<LZBlock>& lz_blocks, const uint32_t min_period, const uint32_t max_period)
        {
            auto end_buckets = get_end_buckets(str, lz_blocks, min_period, max_period);
            
            auto result = get_init_buckets(str.size(), end_buckets);
            
            for (auto& bucket : end_buckets) {
                for (const auto& run : bucket) {
                    result[run.pos].push_back(run);
                }
                bucket.clear();
                bucket.shrink_to_fit();
            }
            
            return result;
        }
        
        // Implements the algorithm described in Kolpakov & Kucherov (1999)
        template <typename T>
        std::vector<std::vector<StringRun>>
        find_maximal_repetitions(const T& str, const uint32_t min_period, const uint32_t max_period)
        {
            std::vector<LZBlock> lz_blocks;
            std::vector<uint32_t> prev_lz_block_occurrence;
            std::tie(lz_blocks, prev_lz_block_occurrence) = lempel_ziv_factorisation_with_prev_block_occurences(str);
            
            auto sorted_buckets = get_sorted_buckets(str, lz_blocks, min_period, max_period);
            
            for (uint32_t k {}; k < lz_blocks.size(); ++k) {
                const auto& block = lz_blocks[k];
                
                auto block_end = block.pos + block.length;
                auto delta     = block.pos - ((prev_lz_block_occurrence[k] == -1) ? 0 : prev_lz_block_occurrence[k]);
                auto v         = block_end - delta;
                
                for (auto j = block.pos; j < block_end; ++j) {
                    const auto& target = sorted_buckets[j - delta];
                    
                    auto last_target = std::lower_bound(std::cbegin(target), std::cend(target), v,
                                                        [] (const auto& run, auto val) {
                                                            return run.pos + run.length < val;
                                                        });
                    
                    auto num_targets = std::distance(std::cbegin(target), last_target);
                    
                    if (num_targets > 0) {
                        std::vector<StringRun> shifted_targets(num_targets);
                        
                        std::transform(std::cbegin(target), last_target, std::begin(shifted_targets),
                                       [delta] (const auto& run) {
                                           return StringRun {run.pos + delta, run.length, run.period};
                                       });
                        
                        sorted_buckets[j].insert(std::begin(sorted_buckets[j]), std::cbegin(shifted_targets),
                                                 std::cend(shifted_targets));
                    }
                }
            }
            
            return sorted_buckets;
        }
        
        template <typename T>
        size_t get_num_runs(const std::vector<T>& buckets)
        {
            return std::accumulate(std::cbegin(buckets), std::cend(buckets), size_t {},
                                   [] (auto lhs, const auto& rhs) { return lhs + rhs.size(); });
        }
        
    } // end namespace detail
    
    template <typename T>
    std::vector<StringRun>
    find_maximal_repetitions(const T& str, const uint32_t min_period = 1, const uint32_t max_period = -1)
    {
        auto sorted_buckets = detail::find_maximal_repetitions(str, min_period, max_period);
        
        std::vector<StringRun> result {};
        result.reserve(detail::get_num_runs(sorted_buckets));
        
        for (auto& bucket : sorted_buckets) {
            result.insert(std::end(result), std::cbegin(bucket), std::cend(bucket));
            bucket.clear();
            bucket.shrink_to_fit();
        }
        
        return result;
    }
    
} // end namespace Tandem

#endif /* defined(__tandem__tandem__) */
