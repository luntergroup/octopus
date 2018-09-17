// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "str_period.hpp"

#include <iterator>
#include <algorithm>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "basics/tandem_repeat.hpp"
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

struct RepeatContextLess
{
    bool operator()(const TandemRepeat& lhs, const TandemRepeat& rhs) const noexcept
    {
        // Expand to discount possible reference pad
        const auto expanded_lhs_region = expand_lhs(mapped_region(lhs), 1);
        const auto expanded_rhs_region = expand_lhs(mapped_region(rhs), 1);
        if (overlap_size(expanded_lhs_region, call_) != overlap_size(expanded_rhs_region, call_)) {
            return overlap_size(expanded_lhs_region, call_) < overlap_size(expanded_rhs_region, call_);
        }
        return ends_before(lhs, rhs);
    }
    RepeatContextLess(const VcfRecord& call) : call_ {call} {};
private:
    const VcfRecord& call_;
};

bool could_contain(const TandemRepeat& repeat, const VcfRecord& call)
{
    return contains(expand(mapped_region(repeat), 1), call);
}

boost::optional<TandemRepeat> find_repeat_context(const VcfRecord& call, const Haplotype& reference)
{
    const auto repeats = find_exact_tandem_repeats(reference.sequence(), reference.mapped_region(), 1, 20);
    const auto overlapping_repeats = overlap_range(repeats, expand(mapped_region(call), 1));
    boost::optional<TandemRepeat> result {};
    if (!empty(overlapping_repeats)) {
        for (const auto& repeat : repeats) {
            if (could_contain(repeat, call)) {
                if (result) {
                    result = std::max(repeat, *result, RepeatContextLess {call});
                } else {
                    result = repeat;
                }
            }
        }
        if (!result) {
            result = *max_overlapped(repeats, call);
        }
    }
    return result;
}

} // namespace

Measure::ResultType STRPeriod::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    int result {0};
    const auto& reference = get_value<ReferenceContext>(facets.at("ReferenceContext"));
    const auto repeat_context = find_repeat_context(call, reference);
    if (repeat_context) result = repeat_context->period();
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
