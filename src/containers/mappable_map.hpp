// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef mappable_map_hpp
#define mappable_map_hpp

#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <stdexcept>
#include <cstddef>
#include <functional>
#include <utility>

#include "mappable_flat_set.hpp"
#include "mappable_flat_multi_set.hpp"

namespace octopus {

template <typename KeyType,
          typename MappableType,
          typename Container = MappableFlatMultiSet<MappableType>>
using MappableMap = std::unordered_map<KeyType, Container>;

template <typename KeyType, typename MappableType>
using MappableSetMap = MappableMap<KeyType, MappableType, MappableFlatSet<MappableType>>;

template <typename Map>
MappableMap<typename Map::key_type, typename Map::mapped_type::value_type>
make_mappable_map(Map map)
{
    using std::make_move_iterator; using std::begin; using std::end;
    using MappableTp = typename Map::mapped_type::value_type;
    
    MappableMap<typename Map::key_type, MappableTp> result {map.size()};
    
    std::for_each(make_move_iterator(begin(map)), make_move_iterator(end(map)),
                  [&result] (auto&& p) {
                      result.emplace(std::piecewise_construct,
                                     std::forward_as_tuple(std::move(p.first)),
                                     std::forward_as_tuple(make_move_iterator(begin(p.second)),
                                                           make_move_iterator(end(p.second))));
                  });
    
    return result;
}

template <typename KeyType, typename Container>
auto sum_region_sizes(const MappableMap<KeyType, typename Container::value_type, Container>& mappables)
{
    using SizeType = typename RegionType<typename Container::value_type>::Position;
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), SizeType {0},
                           [] (const auto curr, const auto& p) {
                               return curr + sum_region_sizes(p.second);
                           });
}

template <typename KeyType, typename Container>
auto encompassing_region(const MappableMap<KeyType, typename Container::value_type, Container>& mappables)
{
    if (mappables.empty()) {
        throw std::runtime_error {"encompassing_region called with empty MappableMap"};
    }
    
    RegionType<typename Container::value_type> result;
    
    const auto it = std::find_if_not(std::cbegin(mappables), std::cend(mappables),
                                     [] (const auto& p) { return p.second.empty(); });
    
    if (it == std::cend(mappables)) {
        throw std::runtime_error {"encompassing_region called with map of empty mappables"};
    }
    
    result = encompassing_region(it->second);
    
    std::for_each(std::next(it), std::cend(mappables),
                  [&result] (const auto& p) {
                      if (!p.second.empty()) {
                          result = encompassing_region(result, encompassing_region(p.second));
                      }
                  });
    
    return result;
}

template <typename KeyType, typename Container>
std::size_t count_mappables(const MappableMap<KeyType, typename Container::value_type, Container>& mappables)
{
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), std::size_t {0},
                           [] (const auto prev, const auto& v) { return prev + v.second.size(); });
}

template <typename KeyType, typename Container, typename MappableType2>
bool
has_overlapped(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
               const MappableType2& mappable)
{
    return std::any_of(std::cbegin(mappables), std::cend(mappables),
                       [&mappable] (const auto& p) {
                           return has_overlapped(p.second, mappable);
                       });
}

template <typename KeyType, typename Container, typename MappableType2>
typename Container::size_type
count_overlapped(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                 const MappableType2& mappable)
{
    using SizeType = typename Container::size_type;
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), SizeType {0},
                           [&mappable] (const auto curr, const auto& p) {
                               return curr + count_overlapped(p.second, mappable);
                           });
}

template <typename KeyType, typename Container, typename MappableType2>
bool
has_contained(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
              const MappableType2& mappable)
{
    return std::any_of(std::cbegin(mappables), std::cend(mappables),
                       [&mappable] (const auto& p) {
                           return has_contained(p.second, mappable);
                       });
}

template <typename KeyType, typename Container, typename MappableType2>
typename Container::size_type
count_contained(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                const MappableType2& mappable)
{
    using SizeType = typename Container::size_type;
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), SizeType {0},
                           [&mappable] (const auto curr, const auto& p) {
                               return curr + count_contained(p.second, mappable);
                           });
}

template <typename KeyType, typename Container, typename MappableType2, typename MappableType3>
bool has_shared(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                const MappableType2& mappable1, const MappableType3& mappable2)
{
    return std::any_of(std::cbegin(mappables), std::cend(mappables),
                       [&mappable1, &mappable2] (const auto& p) {
                           return p.second.has_shared(mappable1, mappable2);
                       });
}

template <typename KeyType, typename Container, typename MappableType2, typename MappableType3>
typename Container::size_type
count_shared(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
             const MappableType2& mappable1, const MappableType3& mappable2)
{
    using SizeType = typename Container::size_type;
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), SizeType {0},
                           [&mappable1, &mappable2] (const auto curr, const auto& p) {
                               return curr + p.second.count_shared(mappable1, mappable2);
                           });
}

template <typename ForwardIterator, typename KeyType, typename Container, typename MappableType2>
ForwardIterator
find_first_shared(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                  ForwardIterator first, ForwardIterator last,
                  const MappableType2& mappable)
{
    using std::cbegin; using std::cend; using std::begin; using std::end;
    
    if (mappables.empty()) return last;
    
    if (mappables.size() == 1) {
        return find_first_shared(cbegin(mappables)->second, first, last, mappable);
    }
    
    std::vector<ForwardIterator> smallest(mappables.size());
    
    std::transform(cbegin(mappables), cend(mappables), begin(smallest),
                   [first, last, &mappable] (const auto& p) {
                       return find_first_shared(p.second, first, last, mappable);
                   });
    
    smallest.erase(std::remove_if(begin(smallest), end(smallest),
                                  [last] (auto it) { return it == last; }),
                   end(smallest));
    
    if (smallest.empty()) return last;
    
    return *std::min_element(cbegin(smallest), cend(smallest), []
                             (auto lhs, auto rhs) { return *lhs < *rhs; });
}

template <typename KeyType, typename Container, typename ForwardIterator>
std::size_t
max_count_if_shared_with_first(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                               ForwardIterator first, ForwardIterator last)
{
    std::size_t maximum {0};
    std::size_t count {0};
    
    for (const auto& map_pair : mappables) {
        count = count_if_shared_with_first(map_pair.second, first, last);
        if (count > maximum) {
            maximum = count;
        }
    }
    
    return maximum;
}

template <typename KeyType, typename Container, typename MappableType2>
typename Container::const_iterator
leftmost_overlapped(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                    const MappableType2& mappable)
{
    if (mappables.empty()) {
        throw std::logic_error {"leftmost_overlapped called with empty MappableMap"};
    }
    
    auto first = std::cbegin(mappables);
    auto last  = std::cend(mappables);
    
    auto result = std::cbegin(first->second);
    
    while (first != last) {
        if (!first->second.empty()) {
            const auto overlapped = overlap_range(first->second, mappable);
            if (!overlapped.empty()) {
                result = std::cbegin(overlapped).base();
                ++first;
                break;
            }
        }
        ++first;
    }
    
    std::for_each(first, last, [&mappable, &result] (const auto& p) {
        const auto overlapped = overlap_range(p.second, mappable);
        if (!overlapped.empty() && begins_before(overlapped.front(), *result)) {
            result = std::cbegin(overlapped).base();
        }
    });
    
    return result;
}

template <typename KeyType, typename Container, typename MappableType2>
typename Container::const_iterator
rightmost_overlapped(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                     const MappableType2& mappable)
{
    if (mappables.empty()) {
        throw std::logic_error {"rightmost_overlapped called with empty MappableMap"};
    }
    
    auto first = std::cbegin(mappables);
    
    const auto last = std::cend(mappables);
    
    auto result = std::cend(first->second);
    
    while (first != last) {
        if (!first->second.empty()) {
            const auto overlapped = overlap_range(first->second, mappable);
            if (!overlapped.empty()) {
                result = rightmost_mappable(overlapped).base();
                ++first;
                break;
            }
        }
        ++first;
    }
    
    std::for_each(first, last, [&mappable, &result] (const auto& p) {
        const auto overlapped = overlap_range(p.second, mappable);
        if (!overlapped.empty()) {
            const auto it = rightmost_mappable(overlapped).base();
            if (ends_before(*result, *it)) {
                result = it;
            }
        }
    });
    
    return result;
}

template <typename KeyType, typename Container>
std::vector<unsigned>
calculate_positional_coverage(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                              const GenomicRegion& region)
{
    std::vector<unsigned> result(region_size(region), 0);
    
    for (const auto& p : mappables) {
        const auto pcoverage = calculate_positional_coverage(p.second, region);
        
        std::transform(std::cbegin(result), std::cend(result), std::cbegin(pcoverage),
                       std::begin(result), std::plus<void> {});
    }
    
    return result;
}

template <typename KeyType, typename Container>
std::vector<unsigned>
calculate_positional_coverage(const MappableMap<KeyType, typename Container::value_type, Container>& mappables)
{
    return calculate_positional_coverage(mappables, encompassing_region(mappables));
}

template <typename KeyType, typename Container, typename MappableType2>
auto
copy_overlapped(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                const MappableType2& mappable)
{
    MappableMap<KeyType, typename Container::value_type, Container> result {mappables.size()};
    
    for (const auto& p : mappables) {
        result.emplace(p.first, copy_overlapped(p.second, mappable));
    }
    
    return result;
}

template <typename KeyType, typename Container, typename MappableType2>
auto
copy_contained(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
               const MappableType2& mappable)
{
    MappableMap<KeyType, typename Container::value_type, Container> result {mappables.size()};
    
    for (const auto& p : mappables) {
        result.emplace(p.first, copy_contained(p.second, mappable));
    }
    
    return result;
}

template <typename KeyType, typename Container, typename MappableType2>
void erase_overlapped(MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                      const MappableType2& mappable)
{
    for (auto& p : mappables) {
        p.second.erase_overlapped(mappable);
    }
}

template <typename KeyType, typename MappableType1, typename MappableType2, typename Container>
void erase_contained(MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                     const MappableType2& mappable)
{
    for (auto& p : mappables) {
        p.second.erase_contained(mappable);
    }
}

template <typename KeyType, typename Container, typename MappableType2>
std::size_t count_spanning(const MappableMap<KeyType, typename Container::value_type, Container>& mappables,
                           const MappableType2& mappable)
{
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), std::size_t {0},
                           [&mappable] (const auto curr, const auto& p) {
                               return curr + count_spanning(p.second, mappable);
                           });
}

} // namespace octopus

#endif
