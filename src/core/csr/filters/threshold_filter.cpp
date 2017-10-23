// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "threshold_filter.hpp"

#include <utility>
#include <numeric>
#include <functional>

#include "io/variant/vcf_header.hpp"

namespace octopus { namespace csr {

ThresholdVariantCallFilter::ThresholdVariantCallFilter(FacetFactory facet_factory,
                                                       std::vector<MeasureWrapper> measures,
                                                       std::vector<std::unique_ptr<Threshold>> thresholds,
                                                       boost::optional<ProgressMeter&> progress)
: VariantCallFilter {std::move(facet_factory), std::move(measures), progress}
, thresholds_ {std::move(thresholds)}
{
    if (measures.size() != thresholds.size()) {
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
                              [] (const auto& measure, const auto& threshold) -> bool { return (*threshold)(measure); });
}

} // namespace csr
} // namespace octopus
