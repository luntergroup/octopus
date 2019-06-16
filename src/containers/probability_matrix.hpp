// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef probability_matrix_hpp
#define probability_matrix_hpp

#include <iterator>
#include <utility>
#include <vector>
#include <functional>

#include "config/common.hpp"
#include "matrix_map.hpp"

namespace octopus {

template <typename T,
          typename Probability = double,
          typename HashT = std::hash<T>,
          typename EqualT = std::equal_to<T>>
using ProbabilityMatrix = MatrixMap<SampleName, T, Probability, std::hash<SampleName>, HashT, std::equal_to<SampleName>, EqualT>;

template <typename T,
          typename Probability,
          typename HashT,
          typename EqualT>
using SampleProbabilities = typename ProbabilityMatrix<T, Probability, HashT, EqualT>::ZipSlice;

template <typename T,
          typename Probability,
          typename HashT,
          typename EqualT>
auto num_samples(const ProbabilityMatrix<T, Probability, HashT, EqualT>& matrix)
{
    return matrix.size1();
}

template <typename T,
          typename Probability,
          typename HashT,
          typename EqualT>
auto extract_keys(const ProbabilityMatrix<T, Probability, HashT, EqualT>& matrix)
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

template <typename T,
          typename Probability,
          typename HashT,
          typename EqualT>
auto extract_key_refs(const ProbabilityMatrix<T, Probability, HashT, EqualT>& matrix)
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

template <typename Container,
          typename T,
          typename Probability,
          typename HashT,
          typename EqualT>
void assign_keys(const Container& keys, ProbabilityMatrix<T, Probability, HashT, EqualT>& matrix)
{
    matrix.assign_keys(std::cbegin(keys), std::cend(keys));
}

template <typename String,
          typename Container,
          typename T,
          typename Probability,
          typename HashT,
          typename EqualT>
void insert_sample(String&& sample, const Container& probabilities, ProbabilityMatrix<T, Probability, HashT, EqualT>& matrix)
{
    matrix.insert_at(std::forward<String>(sample), std::cbegin(probabilities), std::cend(probabilities));
}

template <typename Map,
          typename T,
          typename Probability,
          typename HashT,
          typename EqualT>
void insert_samples(const Map& probabilities, ProbabilityMatrix<T, Probability, HashT, EqualT>& matrix)
{
    matrix.reserve1(probabilities.size());
    for (const auto& s : probabilities) {
        insert_sample(s.first, s.second, matrix);
    }
}

} // namespace octopus

#endif
