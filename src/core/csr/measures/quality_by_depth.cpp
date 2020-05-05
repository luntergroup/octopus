// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "quality_by_depth.hpp"

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "io/variant/vcf_record.hpp"

namespace octopus { namespace csr {

const std::string QualityByDepth::name_ = "QD";

QualityByDepth::QualityByDepth(bool recalculate) : depth_ {recalculate, true} {}

std::unique_ptr<Measure> QualityByDepth::do_clone() const
{
    return std::make_unique<QualityByDepth>(*this);
}

Measure::ValueType QualityByDepth::get_value_type() const
{
    return double {};
}

Measure::ResultType QualityByDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto depth = boost::get<std::size_t>(boost::get<ValueType>(depth_.evaluate(call, facets)));
    Optional<ValueType> result {};
    if (depth > 0) {
        result = static_cast<double>(*call.qual()) / depth;
    }
    return result;
}

Measure::ResultCardinality QualityByDepth::do_cardinality() const noexcept
{
    return depth_.cardinality();
}

const std::string& QualityByDepth::do_name() const
{
    return name_;
}

std::string QualityByDepth::do_describe() const
{
    return "QUAL divided by DP";
}

std::vector<std::string> QualityByDepth::do_requirements() const
{
    return depth_.requirements();
}

bool QualityByDepth::is_equal(const Measure& other) const noexcept
{
    return depth_ == static_cast<const QualityByDepth&>(other).depth_;
}

} // namespace csr
} // namespace octopus
