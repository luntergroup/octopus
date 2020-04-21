// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "gc_content.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/sequence_utils.hpp"
#include "../facets/reference_context.hpp"

namespace octopus { namespace csr {

const std::string GCContent::name_ = "GC";

std::unique_ptr<Measure> GCContent::do_clone() const
{
    return std::make_unique<GCContent>(*this);
}

Measure::ResultType GCContent::get_default_result() const
{
    return double {};
}

Measure::ResultType GCContent::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& reference = get_value<ReferenceContext>(facets.at("ReferenceContext"));
    return utils::gc_content(reference.sequence(expand(mapped_region(call), 50)));
}

Measure::ResultCardinality GCContent::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& GCContent::do_name() const
{
    return name_;
}

std::string GCContent::do_describe() const
{
    return "GC bias of the reference in a window centred on the call";
}

std::vector<std::string> GCContent::do_requirements() const
{
    return {"ReferenceContext"};
}

} // namespace csr
} // namespace octopus
