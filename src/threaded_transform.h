//
//  threaded_transform.h
//  trilobyte
//
//  Created by Daniel Cooke on 29/10/2014.
//  Copyright (c) 2014 Oxford University. All rights reserved.
//

#ifndef trilobyte_threaded_transform_h
#define trilobyte_threaded_transform_h

#include <vector>
#include <thread>
#include <algorithm>
#include <cstddef>

template <typename InputIterator, typename OutputIterator, typename UnaryOperation>
OutputIterator threaded_transform(InputIterator first, InputIterator last, OutputIterator result,
                                  UnaryOperation op, unsigned num_threads)
{
    std::size_t num_values_per_threads = std::distance(first, last) / num_threads;
    
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    
    for (int i = 1; i <= num_threads; ++i) {
        if (i == num_threads) {
            // The last thread processes the remaining values.
            threads.push_back(std::thread(std::transform<InputIterator, OutputIterator, UnaryOperation>,
                                          first, last, result, op));
        } else {
            threads.push_back(std::thread(std::transform<InputIterator, OutputIterator, UnaryOperation>,
                                          first, first + num_values_per_threads, result, op));
        }
        
        first  += num_values_per_threads;
        result += num_values_per_threads;
    }
    
    for (auto& thread : threads)
        thread.join();
    
    return result;
}

template <typename InputIterator, typename OutputIterator, typename UnaryOperation>
OutputIterator threaded_transform(InputIterator first, InputIterator last, OutputIterator result,
                                  UnaryOperation op)
{
    return threaded_transform<InputIterator, OutputIterator, UnaryOperation>
                                (first, last, result, op, std::thread::hardware_concurrency());
}

#endif
