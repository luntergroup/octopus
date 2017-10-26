// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "threshold_filter.hpp"

#include <utility>
#include <numeric>
#include <functional>

#include "io/variant/vcf_header.hpp"

namespace octopus { namespace csr {

auto extract_measures(std::vector<ThresholdVariantCallFilter::Condition>& conditions)
{
    std::vector<MeasureWrapper> result {};
    result.reserve(conditions.size());
    for (auto& condition : conditions) {
        result.push_back(std::move(condition.measure));
    }
    return result;
}

auto extract_thresholds(std::vector<ThresholdVariantCallFilter::Condition>& conditions)
{
    std::vector<ThresholdVariantCallFilter::ThresholdWrapper> result {};
    result.reserve(conditions.size());
    for (auto& condition : conditions) {
        result.push_back(std::move(condition.threshold));
    }
    return result;
}

ThresholdVariantCallFilter::ThresholdVariantCallFilter(FacetFactory facet_factory,
                                                       std::vector<Condition> conditions,
                                                       OutputOptions output_config,
                                                       boost::optional<ProgressMeter&> progress)
: VariantCallFilter {std::move(facet_factory), extract_measures(conditions), output_config, progress}
, thresholds_ {extract_thresholds(conditions)}
{
    if (measures_.size() != thresholds_.size()) {
        throw;
    }
}

void ThresholdVariantCallFilter::annotate(VcfHeader& header) const
{
    // TODO
}

VariantCallFilter::Classification ThresholdVariantCallFilter::classify(const MeasureVector& measures) const
{
    if (passes_all_filters(measures)) {
        return Classification {Classification::Category::unfiltered, {}, boost::none};
    } else {
        return Classification {Classification::Category::filtered, {}, boost::none};
    }
}

bool ThresholdVariantCallFilter::passes_all_filters(const MeasureVector& measures) const
{
    return std::inner_product(std::cbegin(measures), std::cend(measures), std::cbegin(thresholds_), true, std::multiplies<> {},
                              [] (const auto& measure, const auto& threshold) -> bool { return threshold(measure); });
}

} // namespace csr
} // namespace octopus
