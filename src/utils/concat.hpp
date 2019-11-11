// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef concat_hpp
#define concat_hpp

#include <vector>
#include <iterator>

namespace octopus {

template <typename T>
std::vector<T> concat(const std::vector<T>& lhs, const std::vector<T>& rhs)
{
    if (lhs.empty()) return rhs;
    if (rhs.empty()) return lhs;
    std::vector<T> result {};
    result.reserve(lhs.size() + rhs.size());
    result.insert(result.cend(), lhs.cbegin(), lhs.cend());
    result.insert(result.cend(), rhs.cbegin(), rhs.cend());
    return result;
}

template <typename T>
std::vector<T> concat(std::vector<T>&& lhs, const std::vector<T>& rhs)
{
    lhs.insert(lhs.cend(), rhs.cbegin(), rhs.cend());
    return std::move(lhs);
}

template <typename T>
std::vector<T> concat(const std::vector<T>& lhs, std::vector<T>&& rhs)
{
    rhs.insert(rhs.cbegin(), lhs.cbegin(), lhs.cend());
    return std::move(rhs);
}

template <typename T>
std::vector<T> concat(std::vector<T>&& lhs, std::vector<T>&& rhs)
{
    if (lhs.empty()) return std::move(rhs);
    lhs.insert(lhs.cend(), std::make_move_iterator(rhs.begin()), std::make_move_iterator(rhs.end()));
    return std::move(lhs);
}

} // namespace octopus

#endif
