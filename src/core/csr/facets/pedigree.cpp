// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "pedigree.hpp"

#include <utility>

namespace octopus { namespace csr {

const std::string Pedigree::name_ {"Pedigree"};

Pedigree::Pedigree(octopus::Pedigree pedigree)
: pedigree_ {std::move(pedigree)}
{}

Facet::ResultType Pedigree::do_get() const
{
    return std::cref(pedigree_);
}

} // namespace csr
} // namespace octopus
