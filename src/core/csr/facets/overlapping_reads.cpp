// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "overlapping_reads.hpp"

#include <utility>

namespace octopus { namespace csr {

const std::string OverlappingReads::name_ {"OverlappingReads"};

OverlappingReads::OverlappingReads(ReadMap reads) : reads_ {std::move(reads)} {}

Facet::ResultType OverlappingReads::do_get() const
{
    return reads_;
}

} // namespace csr
} // namespace octopus
