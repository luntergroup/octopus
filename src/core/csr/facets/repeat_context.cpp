// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "repeat_context.hpp"

#include <utility>

#include "utils/repeat_finder.hpp"

namespace octopus { namespace csr {

const std::string RepeatContext::name_ {"RepeatContext"};

RepeatContext::RepeatContext(const ReferenceGenome& reference, GenomicRegion region)
: result_ {find_exact_tandem_repeats(reference.fetch_sequence(region), region, 1, 20)}
{}

Facet::ResultType RepeatContext::do_get() const
{
    return std::cref(result_);
}

} // namespace csr
} // namespace octopus
