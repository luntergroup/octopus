// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "model_posterior.hpp"

#include "io/variant/vcf_record.hpp"

namespace octopus { namespace csr {

Measure::ResultType ModelPosterior::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (call.has_info("MP")) {
        return std::stod(call.info_value("MP").front());
    } else {
        return boost::none;
    }
}

std::string ModelPosterior::do_name() const
{
    return "model posterior";
}

} // namespace csr
} // namespace octopus
