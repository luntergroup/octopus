// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "facet.hpp"

namespace octopus { namespace csr {

bool operator==(const Facet& lhs, const Facet& rhs) noexcept
{
    return lhs.name() == rhs.name();
}

bool operator==(const FacetWrapper& lhs, const FacetWrapper& rhs) noexcept
{
    return *lhs.base() == *rhs.base();
}

} // namespace csr
} // namespace octopus
