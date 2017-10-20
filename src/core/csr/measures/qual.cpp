// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "qual.hpp"

#include "io/variant/vcf_record.hpp"

namespace octopus { namespace csr {

Measure::ResultType Qual::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (call.qual()) {
        return static_cast<double>(*call.qual());
    } else {
        return 0.0;
    }
}

std::string Qual::do_name() const
{
    return "qual";
}
    
} // namespace csr
} // namespace octopus
