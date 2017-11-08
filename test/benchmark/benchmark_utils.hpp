// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

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
