/*  tandem.hpp -- Finding tadem repeats in sequences
 
 Copyright (C) 2015 University of Oxford.
 
 Author: Daniel Cooke <dcooke@well.ox.ac.uk>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.  */

#ifndef __tandem__tandem__
#define __tandem__tandem__

#include <vector>
#include <deque>
#include <map>
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <utility>
#include <numeric>
#include <stdexcept>
#include <cassert>

#include "libdivsufsort/divsufsort.h"

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

namespace tandem
{
    struct StringRun
    {
        StringRun() = default;
        
        explicit StringRun(uint32_t pos, uint32_t length, uint32_t period)
        : pos {pos}, length {length}, period {period} {}
        
        uint32_t pos, length, period;
    };
    
    inline bool operator==(const StringRun& lhs, const StringRun& rhs)
    {
        return lhs.pos == rhs.pos && lhs.length == rhs.length;
    }
    inline bool operator!=(const StringRun& lhs, const StringRun& rhs)
    {
        return !operator==(lhs, rhs);
    }
    
    // Wrapper for divsufsort
    template <typename T>
    std::vector<uint32_t> make_suffix_array(const T& str)
    {
        std::vector<saidx_t> sa(str.size()); // divsufsort requires signed integers so need to copy
        
        divsufsort(str.data(), sa.data(), static_cast<int>(str.size()));
        
        return std::vector<uint32_t> {std::cbegin(sa), std::cend(sa)};
    }
    
    template <typename T>
    std::vector<uint32_t> make_suffix_array(const T& str, const size_t extra_capacity)
    {
        std::vector<saidx_t> sa(str.size());
        
        divsufsort(str.data(), sa.data(), static_cast<int>(str.size()));
        
        std::vector<uint32_t> result(str.size() + extra_capacity);
        
        const auto it = std::copy(std::cbegin(sa), std::cend(sa), std::begin(result));
        
        std::fill_n(it, extra_capacity, 0);
        
        return result;
    }
    
    // rank array is inverse suffix array
    std::vector<uint32_t> make_rank_array(const std::vector<uint32_t>& suffix_array);
    std::vector<uint32_t> make_rank_array(const std::vector<uint32_t>& suffix_array,
                                          const std::size_t extra_capacity);
    
    namespace detail
    {
        // Llie et al (2010) show that this direct computation of LCE actually outperforms other O(n) and O(log n)
        // methods in practice.
        
        template <typename T>
        uint32_t forward_lce(const T& str, uint32_t i, uint32_t j, uint32_t n)
        {
            using std::cbegin;
            return static_cast<uint32_t>(std::distance(cbegin(str) + i, std::mismatch(cbegin(str) + i,
                                                                                      cbegin(str) + n,
                                                                                      cbegin(str) + j).first));
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
                                                                                         crend(str) - n,
                                                                                         crend(str) - j - 1).first));
        }
        
        template <typename T>
        uint32_t backward_lce(const T& str, uint32_t i, uint32_t j)
        {
            return backward_lce(str, i, j, uint32_t {});
        }
    } // namespace detail
    
    // LCP = Longest Common Prefix. O(n) implementation given in Kasai et al (2001).
    template <typename T>
    std::vector<uint32_t>
    make_lcp_array(const T& str, const std::vector<uint32_t>& suffix_array,
                   const size_t extra_capacity = 0)
    {
        const auto rank = make_rank_array(suffix_array, extra_capacity);
        
        std::vector<uint32_t> result(suffix_array.size() + extra_capacity);
        
        for (uint32_t i {0}, h {0}; i < (suffix_array.size() - extra_capacity); ++i) {
            if (rank[i] > 0) {
                h += detail::forward_lce(str, i + h, suffix_array[rank[i] - 1] + h);
                result[rank[i]] = h;
                if (h > 0) --h;
            }
        }
        
        return result;
    }
    
    // LPF = Longest Previous Factor
    std::vector<uint32_t>
    make_lpf_array(std::vector<uint32_t> sa, std::vector<uint32_t> lcp);
    std::pair<std::vector<uint32_t>, std::vector<uint32_t>>
    make_lpf_and_prev_occ_arrays(std::vector<uint32_t> sa, std::vector<uint32_t> lcp);
    
    template <typename T>
    auto make_lpf_array(const T& str)
    {
        auto sa  = make_suffix_array(str, 1);
        auto lcp = make_lcp_array(str, sa, 1);
        return make_lpf_array(std::move(sa), std::move(lcp));
    }
    
    template <typename T>
    auto make_lpf_and_prev_occ_arrays(const T& str)
    {
        auto sa  = make_suffix_array(str, 1);
        auto lcp = make_lcp_array(str, sa, 1);
        return make_lpf_and_prev_occ_arrays(std::move(sa), std::move(lcp));
    }
    
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
        
        const auto lpf = make_lpf_array(str);
        
        std::vector<LZBlock> result {};
        result.reserve(str.size()); // max possible blocks
        
        uint32_t end {1};  // start at 1 because the first element of lpf is sentinel
        
        result.emplace_back(0, end);
        
        while (end < str.size()) {
            const auto m = std::max(uint32_t {1}, lpf[end]);
            result.emplace_back(end, m);
            end += m;
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    // Implementation of algorithm found in Crochemore et al. (2008)
    template <typename T>
    std::pair<std::vector<LZBlock>, std::vector<uint32_t>>
    lempel_ziv_factorisation_with_prev_block_occurences(const T& str)
    {
        if (str.empty()) return {{}, {}};
        
        std::vector<uint32_t> lpf, prev_occ;
        std::tie(lpf, prev_occ) = make_lpf_and_prev_occ_arrays(str);
        
        std::vector<LZBlock> lz_blocks {};
        lz_blocks.reserve(str.size()); // max possible blocks
        
        std::vector<uint32_t> prev_lz_block_occurrence {};
        prev_lz_block_occurrence.reserve(str.size());
        
        uint32_t end {1}; // start at 1 because the first element of lpf is sentinel
        
        lz_blocks.emplace_back(0, end);
        prev_lz_block_occurrence.emplace_back(-1);
        
        while (end < str.size()) {
            const auto m = std::max(uint32_t {1}, lpf[end]);
            lz_blocks.emplace_back(end, m);
            prev_lz_block_occurrence.emplace_back(prev_occ[end]);
            end += m;
        }
        
        lz_blocks.shrink_to_fit();
        prev_lz_block_occurrence.shrink_to_fit();
        
        return std::make_pair(std::move(lz_blocks), std::move(prev_lz_block_occurrence));
    }
    
    namespace detail
    {
        // Implements Mains algorithm found in Main (1989). Obscure notation as in paper.
        template <typename T>
        std::deque<StringRun>
        find_leftmost_maximal_repetitions(const T& str, const std::vector<LZBlock>& lz_blocks,
                                          const uint32_t min_period = 1, const uint32_t max_period = -1)
        {
            std::deque<StringRun> result {};
            
            for (uint32_t h {1}; h < lz_blocks.size(); ++h) {
                const auto u   = lz_blocks[h].pos;
                const auto n   = lz_blocks[h].length;
                const auto m   = std::min(u, 2 * lz_blocks[h - 1].length + n);
                const auto t   = u - m;
                const auto end = u + n;
                
                // rightmax periodicities
                for (uint32_t j {min_period}; j <= std::min(n, max_period); ++j) {
                    const auto ls = detail::backward_lce(str, u - 1, u + j - 1, t);
                    const auto lp = detail::forward_lce(str, u + j, u, end);
                    
                    if (ls > 0 && ls + lp >= j && j + lp < n) {
                        result.emplace_back(u - ls, j + lp + ls, j);
                    }
                }
                
                // leftmax periodicities
                for (uint32_t j {min_period}; j < std::min(m, max_period); ++j) {
                    const auto ls = detail::backward_lce(str, u - j - 1, u - 1, t);
                    const auto lp = detail::forward_lce(str, u, u - j, end);
                    
                    if (ls + lp >= j) {
                        result.emplace_back(u - (ls + j), j + lp + ls, j);
                    }
                }
            }
            
            return result;
        }
        
        // just reserves enough space to avoid reallocations
        std::vector<std::vector<StringRun>> get_init_buckets(size_t n, const std::deque<StringRun>& lmrs);
        
        template <typename T>
        std::vector<std::vector<StringRun>>
        get_end_buckets(const T& str, const std::vector<LZBlock>& lz_blocks,
                        const uint32_t min_period, const uint32_t max_period)
        {
            const auto lmrs = find_leftmost_maximal_repetitions(str, lz_blocks, min_period, max_period);
            
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
        std::vector<std::vector<StringRun>>
        get_init_buckets(size_t n, const std::vector<std::vector<StringRun>>& end_buckets);
        
        template <typename T>
        std::vector<std::vector<StringRun>>
        get_sorted_buckets(const T& str, const std::vector<LZBlock>& lz_blocks,
                           const uint32_t min_period, const uint32_t max_period)
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
            
            for (uint32_t k {0}; k < lz_blocks.size(); ++k) {
                const auto& block = lz_blocks[k];
                
                const auto block_end = block.pos + block.length;
                const auto delta     = block.pos - ((prev_lz_block_occurrence[k] == -1) ? 0 : prev_lz_block_occurrence[k]);
                const auto v         = block_end - delta;
                
                for (auto j = block.pos; j < block_end; ++j) {
                    const auto& target = sorted_buckets[j - delta];
                    
                    const auto last_target_itr = std::lower_bound(std::cbegin(target), std::cend(target), v,
                                                                  [] (const auto& run, const auto val) {
                                                                      return run.pos + run.length < val;
                                                                  });
                    
                    const auto num_targets = std::distance(std::cbegin(target), last_target_itr);
                    
                    if (num_targets > 0) {
                        std::vector<StringRun> shifted_targets(num_targets);
                        
                        std::transform(std::cbegin(target), last_target_itr, std::begin(shifted_targets),
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
        size_t count_runs(const std::vector<T>& buckets)
        {
            return std::accumulate(std::cbegin(buckets), std::cend(buckets), size_t {},
                                   [] (const auto curr, const auto& bucket) { return curr + bucket.size(); });
        }
        
        template <typename ForwardIt>
        std::vector<StringRun>
        find_homopolymers(const ForwardIt first, const ForwardIt last)
        {
            std::vector<StringRun> result {};
            
            for (auto curr = first; curr != last; ) {
                const auto it = std::adjacent_find(curr, last);
                
                if (it == last) break;
                
                const auto base = *it;
                
                const auto it2 = std::find_if_not(std::next(it), last,
                                                  [base] (const auto b) { return b == base; });
                
                result.emplace_back(static_cast<std::uint32_t>(std::distance(first, it)),
                                    static_cast<std::uint32_t>(std::distance(it, it2)),
                                    std::uint32_t {1});
                
                curr = it2;
            }
            
            result.shrink_to_fit();
            
            return result;
        }
        
        template <std::size_t N, typename ForwardIt>
        std::vector<StringRun>
        find_exact_tandem_repeats(const ForwardIt first, const ForwardIt last)
        {
            std::vector<StringRun> result {};
            
            if (std::distance(first, last) < 2 * N) return result;
            
            auto it1 = std::adjacent_find(first, last, std::not_equal_to<void> {});
            
            if (it1 == last) return result;
            
            result.reserve(std::min(std::distance(first, last) / N, std::size_t {1024}));
            
            for (auto it2 = std::next(it1, N); it2 < last; ) {
                const auto p = std::mismatch(it2, last, it1);
                
                if (p.second >= it2) {
                    result.emplace_back(static_cast<uint32_t>(std::distance(first, it1)),
                                        static_cast<uint32_t>(std::distance(it1, p.first)),
                                        N);
                    it1 = p.second;
                } else {
                    ++it1;
                }
                
                it1 = std::adjacent_find(it1, last, std::not_equal_to<void> {});
                
                if (it1 == last) break;
                
                it2 = std::next(it1, N);
            }
            
            result.shrink_to_fit();
            
            return result;
        }
        
        template <typename T>
        auto find_homopolymers(const T& str)
        {
            return find_homopolymers(std::cbegin(str), std::cend(str));
        }
        
        template <typename T>
        auto find_exact_dinucleotide_tandem_repeats(const T& str)
        {
            return find_exact_tandem_repeats<2>(std::cbegin(str), std::cend(str));
        }
        
        template <typename T>
        auto find_exact_trinucleotide_tandem_repeats(const T& str)
        {
            return find_exact_tandem_repeats<3>(std::cbegin(str), std::cend(str));
        }
        
        template <typename Container1, typename Container2>
        void append(Container2&& src, Container1& dst)
        {
            const auto it = dst.insert(std::end(dst),
                                       std::make_move_iterator(std::begin(src)),
                                       std::make_move_iterator(std::end(src)));
            
            std::inplace_merge(std::begin(dst), it, std::end(dst),
                               [] (const StringRun& lhs, const StringRun& rhs) {
                                   return lhs.pos < rhs.pos;
                               });
        }
    } // namespace detail
    
    /**
     This is the main function for finding exact repeats in a sequence.
     */
    template <typename T>
    std::vector<StringRun>
    find_maximal_repetitions(const T& str, uint32_t min_period = 1, const uint32_t max_period = -1)
    {
        if (min_period == 0) ++min_period;
        
        if (str.empty() || str.size() < min_period) return {};
        
        if (min_period > max_period) {
            throw std::domain_error {"find_maximal_repetitions: given unsatisfiable condition min_period > max_period"};
        }
        
        if (max_period <= 3) { // The naive algorithm is faster in these cases
            if (min_period == max_period) {
                switch(min_period) {
                    case 1: return detail::find_homopolymers(str);
                    case 2: return detail::find_exact_dinucleotide_tandem_repeats(str);
                    case 3: return detail::find_exact_trinucleotide_tandem_repeats(str);
                }
            }
            
            if (min_period == 1) { // known max_period >= 2
                auto result = detail::find_homopolymers(str);
                
                detail::append(detail::find_exact_dinucleotide_tandem_repeats(str), result);
                
                if (max_period == 3) {
                    detail::append(detail::find_exact_trinucleotide_tandem_repeats(str), result);
                }
                
                return result;
            } else { // min_period == 2 && max_period == 3
                auto result = detail::find_exact_dinucleotide_tandem_repeats(str);
                
                detail::append(detail::find_exact_trinucleotide_tandem_repeats(str), result);
                
                return result;
            }
        }
        
        auto sorted_buckets = detail::find_maximal_repetitions(str, min_period, max_period);
        
        std::vector<StringRun> result {};
        result.reserve(detail::count_runs(sorted_buckets));
        
        for (auto& bucket : sorted_buckets) {
            result.insert(std::end(result), std::cbegin(bucket), std::cend(bucket));
            bucket.clear();
            bucket.shrink_to_fit();
        }
        
        return result;
    }
    
    /**
     Replaces all contiguous sub-sequences of c with a single c, inplace, and returns a map of
     each c position in the new sequence, and how many c's have been removed up to the first non-c
     base past the position. This is just a helper that can speed up repetition finding if the sequence
     contains long runs of characters that are not of interest (e.g. unknwown base 'N' in DNA/RNA sequence).
     
     If this function is used, the output StringRun's will need to be rebased to get the correct positions
     in the original sequence. The function rebase does this transformation.
     
     Example:
     std::string str {"NNNACGTNNTGCNANNNN"};
     auto n_shift_map = colapse(str, 'N'); // str is now "NACGTNTGCNAN", n_shift_map contains (0, 2), (4, 3), (9, 6)
     */
    template <typename SequenceType>
    std::map<size_t, size_t> collapse(SequenceType& sequence, const char c)
    {
        std::map<size_t, size_t> result {};
        
        const auto last = std::end(sequence);
        
        size_t position {}, num_removed {};
        
        for (auto first = std::begin(sequence); first != last;) {
            const auto it1 = std::adjacent_find(first, last,
                                                [c] (const char lhs, const char rhs) {
                                                    return lhs == c && lhs == rhs;
                                                });
            
            if (it1 == last) break;
            
            const auto it2 = std::find_if_not(it1, last, [c] (const char b) { return b == c; });
            
            position    += std::distance(first, it1);
            num_removed += std::distance(it1, it2) - 1;
            
            result.emplace(position, num_removed);
            
            first = it2;
        }
        
        if (!result.empty()) {
            sequence.erase(std::unique(std::next(std::begin(sequence), std::cbegin(result)->first), last,
                                       [c] (const char lhs, const char rhs) {
                                           return lhs == c && lhs == rhs;
                                       }),
                           last);
        }
        
        return result;
    }
    
    void rebase(std::vector<tandem::StringRun>& runs, const std::map<size_t, size_t>& shift_map);
    
} // namespace tandem

#endif /* defined(__tandem__tandem__) */
