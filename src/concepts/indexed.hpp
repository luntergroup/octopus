// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef indexed_hpp
#define indexed_hpp

#include <cstddef>

namespace octopus {

template <typename T, typename IndexType = std::size_t>
class Indexed {};

template <typename T, typename IndexType>
IndexType index_of(const Indexed<T, IndexType>& indexed) noexcept
{
    return static_cast<const T&>(indexed).index();
}

} // namespace octopus

#endif
