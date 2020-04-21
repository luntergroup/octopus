// Copyright (c) 2015-2020 Daniel Cooke
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

Measure::ResultType GenotypeQualityByDepth::get_default_result() const
{
    return std::vector<boost::optional<double>> {};
}

Measure::ResultType GenotypeQualityByDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto depths = boost::get<std::vector<std::size_t>>(depth_.evaluate(call, facets));
    const auto genotype_qualities = boost::get<std::vector<boost::optional<double>>>(GenotypeQuality().evaluate(call, facets));
    assert(depths.size() == genotype_qualities.size());
    std::vector<boost::optional<double>> result(genotype_qualities.size());
    std::transform(std::cbegin(genotype_qualities), std::cend(genotype_qualities), 
                   std::cbegin(depths), std::begin(result),
                   [] (auto gq, auto depth) -> boost::optional<double> {
                        if (gq && depth > 0) { return *gq / depth; } else { return boost::none; } });
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
    
} // namespace csr
} // namespace octopus
