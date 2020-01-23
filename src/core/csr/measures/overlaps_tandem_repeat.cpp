// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "overlaps_tandem_repeat.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/repeat_finder.hpp"
#include "../facets/reference_context.hpp"

namespace octopus { namespace csr {

const std::string OverlapsTandemRepeat::name_ = "STRC";

std::unique_ptr<Measure> OverlapsTandemRepeat::do_clone() const
{
    return std::make_unique<OverlapsTandemRepeat>(*this);
}

Measure::ResultType OverlapsTandemRepeat::get_default_result() const
{
    return bool {};
}

Measure::ResultType OverlapsTandemRepeat::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& reference = get_value<ReferenceContext>(facets.at("ReferenceContext"));
    const auto repeats = find_exact_tandem_repeats(reference.sequence(), reference.mapped_region(), 1, 6);
    return has_overlapped(repeats, call);
}

Measure::ResultCardinality OverlapsTandemRepeat::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& OverlapsTandemRepeat::do_name() const
{
    return name_;
}

std::string OverlapsTandemRepeat::do_describe() const
{
    return "Is the variant in a tandem repeat";
}

std::vector<std::string> OverlapsTandemRepeat::do_requirements() const
{
    return {"ReferenceContext"};
}

} // namespace csr
} // namespace octopus
