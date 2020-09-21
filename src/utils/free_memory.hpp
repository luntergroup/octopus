// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef free_memory_hpp
#define free_memory_hpp

#include <utility>

namespace octopus {

template <typename Container>
void free_memory(Container& c)
{
    Container tmp {};
    using std::swap;
    swap(tmp, c);
}

} // namespace octopus

#endif
