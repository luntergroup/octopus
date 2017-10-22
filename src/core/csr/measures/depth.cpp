// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "depth.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"

namespace octopus { namespace csr {

Measure::ResultType Depth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    return std::stod(call.info_value(vcfspec::info::combinedReadDepth).front());
}

std::string Depth::do_name() const
{
    return "strand bias";
}

} // namespace csr
} // namespace octopus