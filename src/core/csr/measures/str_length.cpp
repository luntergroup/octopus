// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "str_length.hpp"

#include <iterator>
#include <algorithm>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/repeat_finder.hpp"
#include "../facets/reference_context.hpp"

namespace octopus { namespace csr {

const std::string STRLength::name_ = "STR_LENGTH";

std::unique_ptr<Measure> STRLength::do_clone() const
{
    return std::make_unique<STRLength>(*this);
}

namespace {

auto num_periods(const TandemRepeat& repeat) noexcept
{
    return region_size(repeat) / repeat.period;
}

struct PeriodCountLess
{
    bool operator()(const TandemRepeat& lhs, const TandemRepeat& rhs) const noexcept
    {
        return num_periods(lhs) < num_periods(rhs);
    }
};

} // namespace

Measure::ResultType STRLength::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    int result {0};
    const auto& reference = get_value<ReferenceContext>(facets.at("ReferenceContext"));
    const auto repeats = find_exact_tandem_repeats(reference.sequence(), reference.mapped_region(), 1, 6);
    const auto overlapping_repeats = overlap_range(repeats, call);
    if (!empty(overlapping_repeats)) {
        result = num_periods(*std::max_element(std::cbegin(overlapping_repeats), std::cend(overlapping_repeats), PeriodCountLess {}));
    }
    return result;
}

Measure::ResultCardinality STRLength::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& STRLength::do_name() const
{
    return name_;
}

std::string STRLength::do_describe() const
{
    return "Length of overlapping STR";
}

std::vector<std::string> STRLength::do_requirements() const
{
    return {"ReferenceContext"};
}

} // namespace csr
} // namespace octopus
