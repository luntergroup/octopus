//
//  tandem.cpp
//  tandem
//
//  Created by Daniel Cooke on 17/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "tandem.h"

#include <stack>
#include <numeric> // std::accumulate

// Implementation of algorithm found in Crochemore et al. (2008)
std::vector<uint32_t> make_lcf_array(std::vector<uint32_t> sa, std::vector<uint32_t> lcp)
{
    auto n = sa.size();
    
    std::vector<uint32_t> result(n, 0);
    
    sa.push_back(-1);
    lcp.push_back(0);
    
    std::stack<uint32_t> stack {};
    
    stack.push(0);
    
    for (uint32_t i {1}; i <= n; ++i) {
        while (!stack.empty() && (sa[i] == -1 || sa[i] < sa[stack.top()] || (sa[i] > sa[stack.top()] && lcp[i] <= lcp[stack.top()]))) {
            if (sa[i] < sa[stack.top()]) {
                std::tie(lcp[i], result[sa[stack.top()]]) = std::minmax(lcp[stack.top()], uint32_t {lcp[i]});
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
std::pair<std::vector<uint32_t>, std::vector<uint32_t>> make_lpf_and_prev_occ_arrays(std::vector<uint32_t> sa, std::vector<uint32_t> lcp)
{
    auto n = sa.size();
    
    std::vector<uint32_t> lpf(n, 0), prev_occ(n, 0);
    
    sa.push_back(-1);
    lcp.push_back(0);
    
    std::stack<std::pair<uint32_t, uint32_t>> stack {}; // default std::stack uses std::deque
    
    stack.emplace(0, sa[0]);
    
    for (uint32_t i {1}; i <= n; ++i) {
        auto u = lcp[i];
        
        while (!stack.empty() && (sa[i] == -1 || sa[i] < stack.top().second)) {
            std::tie(u, lpf[stack.top().second]) = std::minmax(stack.top().first, uint32_t {u});
            
            auto v = stack.top();
            stack.pop();
            
            if (lpf[v.second] == 0) {
                prev_occ[v.second] = -1;
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
    
    return {lpf, prev_occ};
}

namespace detail
{
    std::size_t get_num_runs(const std::vector<std::deque<StringRun>>& buckets)
    {
        return std::accumulate(std::cbegin(buckets), std::cend(buckets), 0,
                               [] (auto lhs, const auto& rhs) {
                                   return lhs + rhs.size();
                               });
    }
    
    std::size_t get_num_runs(const std::unordered_map<uint32_t, std::deque<StringRun>>& buckets)
    {
        return std::accumulate(std::cbegin(buckets), std::cend(buckets), 0,
                               [] (auto lhs, const auto& rhs) {
                                   return lhs + rhs.second.size();
                               });
    }
} // end namespace detail

void remove_non_primitives(std::vector<StringRun>& repetitions)
{
    auto first = std::begin(repetitions);
    auto last  = std::end(repetitions);
    
    while (first != last) {
        
    }
}
