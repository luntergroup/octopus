// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "posterior_probability_by_depth.hpp"

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "io/variant/vcf_record.hpp"
#include "posterior_probability.hpp"

namespace octopus { namespace csr {

const std::string PosteriorProbabilityByDepth::name_ = "PPD";

PosteriorProbabilityByDepth::PosteriorProbabilityByDepth(bool recalculate) : depth_ {recalculate, true} {}

std::unique_ptr<Measure> PosteriorProbabilityByDepth::do_clone() const
{
    return std::make_unique<PosteriorProbabilityByDepth>(*this);
}

Measure::ValueType PosteriorProbabilityByDepth::get_value_type() const
{
    return double {};
}

Measure::ResultType PosteriorProbabilityByDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    Optional<ValueType> result {};
    const auto posterior_probability = boost::get<Optional<ValueType>>(PosteriorProbability().evaluate(call, facets));
    if (posterior_probability) {
        const auto depth = boost::get<std::size_t>(boost::get<ValueType>(depth_.evaluate(call, facets)));
        if (depth > 0) {
            result = boost::get<double>(*posterior_probability) / depth;
        }
    }
    return result;
}

Measure::ResultCardinality PosteriorProbabilityByDepth::do_cardinality() const noexcept
{
    return depth_.cardinality();
}

const std::string& PosteriorProbabilityByDepth::do_name() const
{
    return name_;
}

std::string PosteriorProbabilityByDepth::do_describe() const
{
    return "PP divided by DP";
}

std::vector<std::string> PosteriorProbabilityByDepth::do_requirements() const
{
    return depth_.requirements();
}

bool PosteriorProbabilityByDepth::is_equal(const Measure& other) const noexcept
{
    return depth_ == static_cast<const PosteriorProbabilityByDepth&>(other).depth_;
}

} // namespace csr
} // namespace octopus
