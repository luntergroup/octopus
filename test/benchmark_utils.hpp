//
//  benchmark_utils.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_benchmark_utils_hpp
#define Octopus_benchmark_utils_hpp

#include <chrono>
#include <vector>
#include <numeric>

template <typename D=std::chrono::nanoseconds, typename F>
D benchmark(F f, unsigned num_tests)
{
    D total {0};
    for (unsigned i = 0; i < num_tests; ++i) {
        auto start = std::chrono::system_clock::now();
        f();
        auto end = std::chrono::system_clock::now();
        total += std::chrono::duration_cast<D>(end - start);
    }
    return D {total / num_tests};
}

#endif
