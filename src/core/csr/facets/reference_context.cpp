// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "reference_context.hpp"

#include <utility>

namespace octopus { namespace csr {

const std::string ReferenceContext::name_ {"ReferenceContext"};

ReferenceContext::ReferenceContext(const ReferenceGenome& reference, GenomicRegion region)
: reference_ {reference}
, region_ {std::move(region)}
, result_ {}
{}

Facet::ResultType ReferenceContext::do_get() const
{
    if (!result_) {
        result_ = Haplotype {region_, reference_};
    }
    return *result_;
}

} // namespace csr
} // namespace octopus
