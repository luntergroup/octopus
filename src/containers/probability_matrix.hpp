// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef probability_matrix_hpp
#define probability_matrix_hpp

#include <iterator>
#include <utility>
#include <vector>
#include <functional>

#include <config/common.hpp>

#include "matrix_map.hpp"

namespace octopus {

template <typename T>
using ProbabilityMatrix = MatrixMap<SampleName, T, double>;

template <typename T>
using SampleProbabilities = typename ProbabilityMatrix<T>::ZipSlice;

template <typename T>
auto num_samples(const ProbabilityMatrix<T>& matrix)
{
    return matrix.size1();
}

template <typename T>
auto extract_keys(const ProbabilityMatrix<T>& matrix)
{
    std::vector<T> result {};
    
    if (matrix.empty1()) return result;
    
    result.reserve(matrix.size2());
    
    const auto& test_sample = matrix.begin()->first;
    
    for (const auto& p : matrix[test_sample]) {
        result.emplace_back(p.first);
    }
    
    return result;
}

template <typename T>
auto extract_key_refs(const ProbabilityMatrix<T>& matrix)
{
    std::vector<std::reference_wrapper<const T>> result {};
    
    if (matrix.empty1()) return result;
    
    result.reserve(matrix.size2());
    
    const auto& test_sample = matrix.begin()->first;
    
    for (const auto& p : matrix[test_sample]) {
        result.emplace_back(p.first);
    }
    
    return result;
}

template <typename Container, typename T>
void assign_keys(const Container& keys, ProbabilityMatrix<T>& matrix)
{
    matrix.assign_keys(std::cbegin(keys), std::cend(keys));
}

template <typename S, typename Container, typename T>
void insert_sample(S&& sample, const Container& probabilities, ProbabilityMatrix<T>& matrix)
{
    matrix.insert_at(std::forward<S>(sample), std::cbegin(probabilities), std::cend(probabilities));
}

template <typename Map, typename T>
void insert_samples(const Map& probabilities, ProbabilityMatrix<T>& matrix)
{
    matrix.reserve1(probabilities.size());
    for (const auto& s : probabilities) {
        insert_sample(s.first, s.second, matrix);
    }
}

} // namespace octopus

#endif
