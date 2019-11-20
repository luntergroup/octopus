// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef reorder_hpp
#define reorder_hpp

namespace octopus {

template <typename order_iterator, typename value_iterator>
void reorder(order_iterator order_begin, order_iterator order_end, value_iterator v)
{
    // See https://stackoverflow.com/a/1267878/2970186
    using value_t = typename std::iterator_traits< value_iterator >::value_type;
    using index_t = typename std::iterator_traits< order_iterator >::value_type;
    using diff_t  = typename std::iterator_traits< order_iterator >::difference_type;
    diff_t remaining = order_end - 1 - order_begin;
    for (index_t s = index_t {}, d; remaining > 0; ++s) {
        for (d = order_begin[s]; d > s; d = order_begin[d]);
        if (d == s) {
            --remaining;
            value_t temp = std::move(v[s]);
            while (d = order_begin[d], d != s) {
                std::swap(temp, v[d]);
                --remaining;
            }
            v[s] = std::move(temp);
        }
    }
}

} // namespace octopus

#endif
