// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef append_hpp
#define append_hpp

#include <vector>
#include <deque>
#include <iterator>
#include <utility>

#include "containers/mappable_block.hpp"

namespace octopus { namespace utils {

template <typename T>
auto append(std::vector<T>& v, typename std::vector<T>::size_type count, const T& value)
{
    return v.insert(std::cend(v), count, value);
}

template <typename T>
auto append(const std::vector<T>& src, std::vector<T>& dest)
{
    typename std::vector<T>::iterator result;
    
    if (dest.empty()) {
        dest   = src;
        result = std::begin(dest);
    } else {
        result = dest.insert(std::cend(dest), std::cbegin(src), std::cend(src));
    }
    
    return result;
}

template <typename T>
auto append(std::vector<T>&& src, std::vector<T>& dest)
{
    typename std::vector<T>::iterator result;
    
    if (dest.empty()) {
        dest   = std::move(src);
        result = std::begin(dest);
    } else {
        result = dest.insert(std::cend(dest),
                             std::make_move_iterator(std::begin(src)),
                             std::make_move_iterator(std::end(src)));
    }
    
    src.clear();
    src.shrink_to_fit();
    
    return result;
}

template <typename T>
auto append(const std::deque<T>& src, std::deque<T>& dest)
{
    typename std::deque<T>::iterator result;
    
    if (dest.empty()) {
        dest   = src;
        result = std::begin(dest);
    } else {
        result = dest.insert(std::cend(dest), std::cbegin(src), std::cend(src));
    }
    
    return result;
}

template <typename T>
auto append(std::deque<T>&& src, std::deque<T>& dest)
{
    typename std::deque<T>::iterator result;
    
    if (dest.empty()) {
        dest   = std::move(src);
        result = std::begin(dest);
    } else {
        result = dest.insert(std::cend(dest),
                             std::make_move_iterator(std::begin(src)),
                             std::make_move_iterator(std::end(src)));
    }
    
    src.clear();
    src.shrink_to_fit();
    
    return result;
}

template <typename T>
auto append(const std::deque<T>& src, std::vector<T>& dest)
{
    return dest.insert(std::cend(dest), std::cbegin(src), std::cend(src));
}

template <typename T>
auto append(std::deque<T>&& src, std::vector<T>& dest)
{
    const auto result = dest.insert(std::cend(dest),
                                    std::make_move_iterator(std::begin(src)),
                                    std::make_move_iterator(std::end(src)));
    
    src.clear();
    src.shrink_to_fit();
    
    return result;
}

template <typename T>
auto append(const std::vector<T>& src, std::deque<T>& dest)
{
    return dest.insert(std::cend(dest), std::cbegin(src), std::cend(src));
}

template <typename T>
auto append(std::vector<T>&& src, std::deque<T>& dest)
{
    const auto result = dest.insert(std::cend(dest),
                                    std::make_move_iterator(std::begin(src)),
                                    std::make_move_iterator(std::end(src)));
    
    src.clear();
    src.shrink_to_fit();
    
    return result;
}

template <typename T>
auto append(const MappableBlock<T>& src, MappableBlock<T>& dest)
{
    typename std::vector<T>::iterator result;
    if (dest.empty()) {
        dest   = src;
        result = std::begin(dest);
    } else {
        result = dest.insert(std::cend(dest), std::cbegin(src), std::cend(src));
    }
    return result;
}

template <typename T>
auto append(MappableBlock<T>&& src, MappableBlock<T>& dest)
{
    typename std::vector<T>::iterator result;
    if (dest.empty()) {
        dest   = std::move(src);
        result = std::begin(dest);
    } else {
        result = dest.insert(std::cend(dest),
                             std::make_move_iterator(std::begin(src)),
                             std::make_move_iterator(std::end(src)));
    }
    src.clear();
    src.shrink_to_fit();
    return result;
}

} // namespace utils
} // namespace octopus

#endif
