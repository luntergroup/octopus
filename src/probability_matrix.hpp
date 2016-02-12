//
//  probability_matrix.hpp
//  Octopus
//
//  Created by Daniel Cooke on 04/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef probability_matrix_hpp
#define probability_matrix_hpp

#include <iterator>
#include <utility>

#include "common.hpp"
#include "matrix_map.hpp"

namespace Octopus
{
    template <typename T>
    using ProbabilityMatrix = MatrixMap<SampleIdType, T, double>;
    
    template <typename T>
    using SampleProbabilities = typename ProbabilityMatrix<T>::ZipSlice;
    
    template <typename Container, typename T>
    void assign_genotypes(const Container& genotypes,
                          ProbabilityMatrix<T>& matrix)
    {
        matrix.assign_keys(std::cbegin(genotypes), std::cend(genotypes));
    }
    
    template <typename S, typename Container, typename T>
    void insert_sample(S&& sample, const Container& probabilities,
                       ProbabilityMatrix<T>& matrix)
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
} // namespace Octopus

#endif /* probability_matrix_hpp */
