// Copyright (c) 2015-2019 Daniel Cooke
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

Measure::ResultType PosteriorProbabilityByDepth::get_default_result() const
{
    return boost::optional<double> {};
}

Measure::ResultType PosteriorProbabilityByDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    boost::optional<double> result {};
    const auto posterior_probability = boost::get<boost::optional<double>>(PosteriorProbability().evaluate(call, facets));
    if (posterior_probability) {
        auto depth = boost::get<std::size_t>(depth_.evaluate(call, facets));
        if (depth > 0) {
            result = *posterior_probability / depth;
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
