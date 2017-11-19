// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "gc_content.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/sequence_utils.hpp"
#include "../facets/reference_context.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> GCContent::do_clone() const
{
    return std::make_unique<GCContent>(*this);
}

Measure::ResultType GCContent::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto reference = boost::get<ReferenceContext::ResultType>(facets.at("ReferenceContext").get());
    return utils::gc_content(reference);
}

std::string GCContent::do_name() const
{
    return "GC";
}

std::vector<std::string> GCContent::do_requirements() const
{
    return {"ReferenceContext"};
}

} // namespace csr
} // namespace octopus
