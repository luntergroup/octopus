// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef random_select_hpp
#define random_select_hpp

#include <random>
#include <iterator>
#include <cstddef>

#include <boost/random/uniform_int_distribution.hpp>

namespace octopus {

template <typename ForwardIt, typename RandomGenerator>
ForwardIt random_select(ForwardIt first, ForwardIt last, RandomGenerator& g)
{
    if (first == last) return first;
    const auto max = static_cast<std::size_t>(std::distance(first, last));
    if (max == 1) return first;
    boost::random::uniform_int_distribution<std::size_t> dist {0, max - 1};
    std::advance(first, dist(g));
    return first;
}

template <typename ForwardIt>
ForwardIt random_select(ForwardIt first, ForwardIt last)
{
    static thread_local std::mt19937 generator {42};
    return random_select(first, last, generator);
}

template <typename Range>
decltype(auto) random_select(const Range& values)
{
    assert(!values.empty());
    return *random_select(std::cbegin(values), std::cend(values));
}

} // namespace octopus

#endif
