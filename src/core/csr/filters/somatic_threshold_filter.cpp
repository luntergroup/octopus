// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_threshold_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <cassert>

#include <boost/variant.hpp>

#include "../measures/is_somatic.hpp"

namespace octopus { namespace csr {

SomaticThresholdVariantCallFilter::SomaticThresholdVariantCallFilter(FacetFactory facet_factory,
                                                                     ConditionVectorPair germline,
                                                                     ConditionVectorPair somatic,
                                                                     OutputOptions output_config,
                                                                     ConcurrencyPolicy threading,
                                                                     boost::optional<ProgressMeter&> progress)
: ConditionalThresholdVariantCallFilter {
    std::move(facet_factory),
    {std::move(germline), std::move(somatic)},
    {make_wrapped_measure<IsSomatic>(true)},
    [] (const auto& measures) -> std::size_t { return boost::get<bool>(measures.back()); },
    output_config, threading, progress
} {}

bool SomaticThresholdVariantCallFilter::is_soft_filtered(const ClassificationList& sample_classifications,
                                                         const MeasureVector& measures) const
{
    return std::any_of(std::cbegin(sample_classifications), std::cend(sample_classifications),
                       [] (const auto& c) { return c.category != Classification::Category::unfiltered; });
}

} // namespace csr
} // namespace octopus
