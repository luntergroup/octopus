// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "str_period.hpp"

#include <iterator>
#include <algorithm>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/repeat_finder.hpp"
#include "../facets/reference_context.hpp"

namespace octopus { namespace csr {

const std::string STRPeriod::name_ = "STR_PERIOD";

std::unique_ptr<Measure> STRPeriod::do_clone() const
{
    return std::make_unique<STRPeriod>(*this);
}

namespace {

struct PeriodLess
{
    bool operator()(const TandemRepeat& lhs, const TandemRepeat& rhs) const noexcept
    {
        return lhs.period < rhs.period;
    }
};

} // namespace

Measure::ResultType STRPeriod::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    int result {0};
    const auto& reference = get_value<ReferenceContext>(facets.at("ReferenceContext"));
    const auto repeats = find_exact_tandem_repeats(reference.sequence(), reference.mapped_region(), 1, 6);
    const auto overlapping_repeats = overlap_range(repeats, call);
    if (!empty(overlapping_repeats)) {
        result = std::max_element(std::cbegin(overlapping_repeats), std::cend(overlapping_repeats), PeriodLess {})->period;
    }
    return result;
}

Measure::ResultCardinality STRPeriod::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& STRPeriod::do_name() const
{
    return name_;
}

std::string STRPeriod::do_describe() const
{
    return "Length of overlapping STR";
}

std::vector<std::string> STRPeriod::do_requirements() const
{
    return {"ReferenceContext"};
}

} // namespace csr
} // namespace octopus
