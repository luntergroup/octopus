// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "reference_context.hpp"

#include <utility>

namespace octopus { namespace csr {

const std::string ReferenceContext::name_ {"ReferenceContext"};

ReferenceContext::ReferenceContext(const ReferenceGenome& reference, GenomicRegion region)
: result_ {std::make_unique<Haplotype>(region, reference)}
{}

Facet::ResultType ReferenceContext::do_get() const
{
    return std::cref(*result_);
}

} // namespace csr
} // namespace octopus
