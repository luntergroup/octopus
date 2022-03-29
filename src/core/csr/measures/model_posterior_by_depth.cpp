// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "model_posterior_by_depth.hpp"

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/concat.hpp"
#include "model_posterior.hpp"

namespace octopus { namespace csr {

const std::string ModelPosteriorByDepth::name_ = "MPD";

ModelPosteriorByDepth::ModelPosteriorByDepth(bool recalculate) : depth_ {recalculate, true} {}

std::unique_ptr<Measure> ModelPosteriorByDepth::do_clone() const
{
    return std::make_unique<ModelPosteriorByDepth>(*this);
}

Measure::ValueType ModelPosteriorByDepth::get_value_type() const
{
    return double {};
}

Measure::ResultType ModelPosteriorByDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto depths = boost::get<Array<ValueType>>(depth_.evaluate(call, facets));
    const auto model_posteriors = boost::get<Array<Optional<ValueType>>>(ModelPosterior().evaluate(call, facets));
    assert(depths.size() == model_posteriors.size());
    Array<Optional<ValueType>> result(model_posteriors.size());
    const auto model_posterior_by_depth = [] (const auto gq, const ValueType depth) {
        Optional<ValueType> result {};
        if (gq && boost::get<std::size_t>(depth) > 0) {
            result = boost::get<double>(*gq) / boost::get<std::size_t>(depth);
        }
        return result;
    };
    std::transform(std::cbegin(model_posteriors), std::cend(model_posteriors), 
                   std::cbegin(depths), std::begin(result), model_posterior_by_depth);
    return result;
}

Measure::ResultCardinality ModelPosteriorByDepth::do_cardinality() const noexcept
{
    return depth_.cardinality();
}

const std::string& ModelPosteriorByDepth::do_name() const
{
    return name_;
}

std::string ModelPosteriorByDepth::do_describe() const
{
    return "MP divided by DP";
}

std::vector<std::string> ModelPosteriorByDepth::do_requirements() const
{
    return depth_.requirements();
}

bool ModelPosteriorByDepth::is_equal(const Measure& other) const noexcept
{
    return depth_ == static_cast<const ModelPosteriorByDepth&>(other).depth_;
}

} // namespace csr
} // namespace octopus
