//
//  benchmark_utils.hpp
//  octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_benchmark_utils_hpp
#define Octopus_benchmark_utils_hpp

#include <chrono>

template <typename D = std::chrono::nanoseconds, typename F>
D benchmark(F f, unsigned num_tests)
{
    D total {0};
    
    for (; num_tests > 0; --num_tests) {
        const auto start = std::chrono::system_clock::now();
        f();
        const auto end = std::chrono::system_clock::now();
        total += std::chrono::duration_cast<D>(end - start);
    }
    
    return D {total / num_tests};
}

#endif
