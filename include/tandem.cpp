//
//  tandem.cpp
//  tandem
//
//  Created by Daniel Cooke on 17/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "tandem.h"

#include <stack>

// Implementation of algorithm found in Crochemore et al. (2008)
std::vector<size_t> make_longest_common_factor_array(std::vector<size_t> sa, std::vector<size_t> lcp)
{
    auto n = sa.size();
    
    std::vector<size_t> result(n, 0);
    
    sa.push_back(-1);
    lcp.push_back(0);
    
    std::stack<size_t> stack {};
    
    stack.push(0);
    
    for (size_t i {1}; i <= n; ++i) {
        while (!stack.empty() && (sa[i] == -1 || sa[i] < sa[stack.top()] || (sa[i] > sa[stack.top()] && lcp[i] <= lcp[stack.top()]))) {
            if (sa[i] < sa[stack.top()]) {
                auto tmp = lcp[i];
                std::tie(lcp[i], result[sa[stack.top()]]) = std::minmax(lcp[stack.top()], tmp);
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
std::pair<std::vector<size_t>, std::vector<size_t>> make_lpf_and_prev_occ(std::vector<size_t> sa, std::vector<size_t> lcp)
{
    auto n = sa.size();
    
    std::vector<size_t> lpf(n, 0), prev_occ(n, 0);
    
    sa.push_back(-1);
    lcp.push_back(0);
    
    std::stack<std::pair<size_t, size_t>> stack {};
    
    stack.emplace(0, sa[0]);
    
    for (size_t i {1}; i <= n; ++i) {
        auto u = lcp[i];
        
        while (!stack.empty() && (sa[i] == -1 || sa[i] < stack.top().second)) {
            auto tmp = u;
            std::tie(u, lpf[stack.top().second]) = std::minmax(stack.top().first, tmp);
            
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
