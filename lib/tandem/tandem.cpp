/*  tandem.cpp
 
 Copyright (C) 2017-2018 University of Oxford.
 
 Author: Daniel Cooke <dcooke@well.ox.ac.uk>
 
 Use of this source code is governed by the MIT license that can be found in the LICENSE file. */

#include "tandem.hpp"

#include <stack>

namespace tandem  {

std::vector<std::uint32_t> make_rank_array(const std::vector<std::uint32_t>& suffix_array)
{
    std::vector<std::uint32_t> result(suffix_array.size());
    for (std::uint32_t i {0}; i < suffix_array.size(); ++i) {
        result[suffix_array[i]] = i;
    }
    return result;
}

std::vector<std::uint32_t> make_rank_array(const std::vector<std::uint32_t>& suffix_array,
                                           const std::size_t extra_capacity)
{
    std::vector<std::uint32_t> result(suffix_array.size() - extra_capacity);
    for (std::uint32_t i {0}; i < (suffix_array.size() - extra_capacity); ++i) {
        result[suffix_array[i]] = i;
    }
    return result;
}

// Implementation of algorithm found in Crochemore et al. (2008)
std::vector<std::uint32_t> make_lpf_array(std::vector<std::uint32_t> sa, std::vector<std::uint32_t> lcp)
{
    const auto n = sa.size();
    
    std::vector<std::uint32_t> result(n, 0);
    
    static constexpr auto sentinal = std::numeric_limits<std::uint32_t>::max();
    sa.push_back(sentinal);
    lcp.push_back(0);
    
    std::stack<std::uint32_t> stack {};
    stack.push(0);
    
    for (std::uint32_t i {1}; i <= n; ++i) {
        while (!stack.empty()
               && (sa[i] == sentinal || sa[i] < sa[stack.top()]
                   || (sa[i] > sa[stack.top()] && lcp[i] <= lcp[stack.top()]))) {
                   if (sa[i] < sa[stack.top()]) {
                       std::tie(lcp[i], result[sa[stack.top()]]) = std::minmax(lcp[stack.top()], std::uint32_t {lcp[i]});
                   } else {
                       result[sa[stack.top()]] = lcp[stack.top()];
                   }
                   stack.pop();
               }
        if (i < n) {
            stack.push(i);
        }
    }
    
    return result;
}

// Implementation of algorithm found in Crochemore et al. (2008)
std::pair<std::vector<std::uint32_t>, std::vector<std::uint32_t>>
make_lpf_and_prev_occ_arrays(std::vector<std::uint32_t> sa, std::vector<std::uint32_t> lcp)
{
    const auto n = sa.size();
    
    std::vector<std::uint32_t> lpf(n, 0), prev_occ(n, 0);
    
    static constexpr auto sentinal = std::numeric_limits<std::uint32_t>::max();
    sa.push_back(sentinal);
    lcp.push_back(0);
    
    std::stack<std::pair<std::uint32_t, std::uint32_t>> stack {}; // default std::stack uses std::deque
    stack.emplace(0, sa[0]);
    
    for (std::uint32_t i {1}; i <= n; ++i) {
        auto u = lcp[i];
        while (!stack.empty() && (sa[i] == sentinal || sa[i] < stack.top().second)) {
            std::tie(u, lpf[stack.top().second]) = std::minmax(stack.top().first, std::uint32_t {u});
            const auto v = stack.top();
            stack.pop();
            if (lpf[v.second] == 0) {
                prev_occ[v.second] = sentinal;
            } else if (v.first > u) {
                prev_occ[v.second] = stack.top().second;
            } else {
                prev_occ[v.second] = sa[i];
            }
        }
        if (i < n) {
            stack.emplace(u, sa[i]);
        }
    }
    
    return std::make_pair(std::move(lpf), std::move(prev_occ));
}

namespace detail {

std::vector<std::vector<Repeat>>
get_init_buckets(const std::size_t n, const LMRVector& lmrs)
{
    std::vector<std::uint32_t> counts(n, 0);
    for (const auto& run : lmrs) {
        ++counts[run.pos + run.length - 1];
    }
    std::vector<std::vector<Repeat>> result(n, std::vector<Repeat> {});
    for (std::size_t i {0}; i < n; ++i) {
        result[i].reserve(counts[i]);
    }
    return result;
}

std::vector<std::vector<Repeat>>
get_init_buckets(const std::size_t n, const std::vector<std::vector<Repeat>>& end_buckets)
{
    std::vector<std::uint32_t> counts(n, 0);
    for (const auto& bucket : end_buckets) {
        for (const auto& run : bucket) {
            ++counts[run.pos];
        }
    }
    std::vector<std::vector<Repeat>> result(n, std::vector<Repeat> {});
    for (std::size_t i {0}; i < n; ++i) {
        result[i].reserve(counts[i]);
    }
    return result;
}

} // namespace detail

void rebase(std::vector<tandem::Repeat>& runs, const std::map<std::size_t, std::size_t>& shift_map)
{
    if (shift_map.empty()) return;
    auto shift_map_it = std::cbegin(shift_map);
    for (auto& run : runs) {
        while (std::next(shift_map_it) != std::cend(shift_map) && std::next(shift_map_it)->first <= run.pos) {
            ++shift_map_it;
        }
        run.pos += static_cast<decltype(run.pos)>(shift_map_it->second);
    }
}

} // namespace tandem
