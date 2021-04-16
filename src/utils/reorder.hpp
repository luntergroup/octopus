// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef reorder_hpp
#define reorder_hpp

#include <vector>
#include <iterator>
#include <algorithm>
#include <utility>

namespace octopus {

template <typename OrderIterator, typename ValueIterator>
ValueIterator reorder(OrderIterator order_begin, OrderIterator order_end, ValueIterator v)
{
    using ValueType = typename std::iterator_traits<ValueIterator>::value_type;
    const auto n = std::distance(order_begin, order_end);
    std::vector<ValueType> buffer {};
    buffer.reserve(n);
    std::transform(order_begin, order_end, std::back_inserter(buffer), [&] (auto idx) { return std::move(*std::next(v, idx)); });
    return std::copy_n(std::make_move_iterator(std::begin(buffer)), n, v);
}

} // namespace octopus

#endif
