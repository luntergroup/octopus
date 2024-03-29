// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "genotype_quality_by_depth.hpp"

#include <cassert>
#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "io/variant/vcf_record.hpp"
#include "utils/concat.hpp"
#include "genotype_quality.hpp"

namespace octopus { namespace csr {

const std::string GenotypeQualityByDepth::name_ = "GQD";

GenotypeQualityByDepth::GenotypeQualityByDepth(bool recalculate) : depth_ {recalculate, false} {}

std::unique_ptr<Measure> GenotypeQualityByDepth::do_clone() const
{
    return std::make_unique<GenotypeQualityByDepth>(*this);
}

Measure::ValueType GenotypeQualityByDepth::get_value_type() const
{
    return double {};
}

Measure::ResultType GenotypeQualityByDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto depths = boost::get<Array<ValueType>>(depth_.evaluate(call, facets));
    const auto genotype_qualities = boost::get<Array<Optional<ValueType>>>(GenotypeQuality().evaluate(call, facets));
    assert(depths.size() == genotype_qualities.size());
    Array<Optional<ValueType>> result(genotype_qualities.size());
    const auto genotype_quality_by_depth = [] (const auto gq, const ValueType depth) {
        Optional<ValueType> result {};
        if (gq && boost::get<std::size_t>(depth) > 0) {
            result = boost::get<double>(*gq) / boost::get<std::size_t>(depth);
        }
        return result;
    };
    std::transform(std::cbegin(genotype_qualities), std::cend(genotype_qualities), 
                   std::cbegin(depths), std::begin(result), genotype_quality_by_depth);
    return result;
}

Measure::ResultCardinality GenotypeQualityByDepth::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& GenotypeQualityByDepth::do_name() const
{
    return name_;
}

std::string GenotypeQualityByDepth::do_describe() const
{
    return "GQ divided by DP";
}

std::vector<std::string> GenotypeQualityByDepth::do_requirements() const
{
    return concat(depth_.requirements(), GenotypeQuality().requirements());
}

bool GenotypeQualityByDepth::is_equal(const Measure& other) const noexcept
{
    return depth_ == static_cast<const GenotypeQualityByDepth&>(other).depth_;
}

void GenotypeQualityByDepth::do_set_parameters(std::vector<std::string> params)
{
    depth_.set_parameters(std::move(params));
}

std::vector<std::string> GenotypeQualityByDepth::do_parameters() const
{
    return depth_.parameters();
}

} // namespace csr
} // namespace octopus
