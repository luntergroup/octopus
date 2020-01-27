// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef erase_if_hpp
#define erase_if_hpp

#include <vector>
#include <deque>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <utility>

namespace octopus {

template <typename T, typename Alloc, typename U>
void erase(std::vector<T, Alloc>& c, const U& value)
{
    c.erase(std::remove(std::begin(c), std::end(c), value), std::end(c));
}

template <typename T, typename Alloc, typename UnaryPredicate>
void erase_if(std::vector<T, Alloc>& c, UnaryPredicate&& pred)
{
    c.erase(std::remove_if(std::begin(c), std::end(c), std::forward<UnaryPredicate>(pred)), std::end(c));
}

template <typename T, typename Alloc, typename U>
void erase(std::deque<T, Alloc>& c, const U& value)
{
    c.erase(std::remove(std::begin(c), std::end(c), value), std::end(c));
}

template <typename T, typename Alloc, typename UnaryPredicate>
void erase_if(std::deque<T, Alloc>& c, UnaryPredicate&& pred)
{
    c.erase(std::remove_if(std::begin(c), std::end(c), std::forward<UnaryPredicate>(pred)), std::end(c));
}

template <typename Key, typename T, typename Compare, typename Alloc, typename Pred>
void erase_if(std::map<Key, T ,Compare, Alloc>& c, Pred pred)
{
    for (auto i = std::cbegin(c), last = std::cend(c); i != last;) {
        if (pred(*i)) {
            i = c.erase(i);
        } else {
            ++i;
        }
    }
}

template <typename Key, typename T, typename Hash, typename KeyEqual, typename Alloc, typename Pred>
void erase_if(std::unordered_map<Key, T, Hash, KeyEqual, Alloc>& c, Pred pred)
{
    for (auto i = std::cbegin(c), last = std::cend(c); i != last;) {
        if (pred(*i)) {
            i = c.erase(i);
        } else {
            ++i;
        }
    }
}

} // namespace

#endif //erase_if_hpp
