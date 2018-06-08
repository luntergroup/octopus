// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "genotypes.hpp"

#include <utility>

namespace octopus { namespace csr {

const std::string Genotypes::name_ {"Genotypes"};

Genotypes::Genotypes(GenotypeMap genotypes) : genotypes_ {std::move(genotypes)} {}

Facet::ResultType Genotypes::do_get() const
{
    return std::cref(genotypes_);
}

} // namespace csr
} // namespace octopus
