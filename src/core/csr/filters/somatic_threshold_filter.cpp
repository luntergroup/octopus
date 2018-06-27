// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_threshold_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <cassert>

#include <boost/variant.hpp>

#include "../measures/is_somatic.hpp"
#include "../measures/is_refcall.hpp"

namespace octopus { namespace csr {

SomaticThresholdVariantCallFilter::SomaticThresholdVariantCallFilter(FacetFactory facet_factory,
                                                                     ConditionVectorPair germline,
                                                                     ConditionVectorPair somatic,
                                                                     ConditionVectorPair reference,
                                                                     OutputOptions output_config,
                                                                     ConcurrencyPolicy threading,
                                                                     boost::optional<ProgressMeter&> progress)
: ConditionalThresholdVariantCallFilter {
    std::move(facet_factory),
    {std::move(germline), std::move(somatic), std::move(reference)},
    {make_wrapped_measure<IsSomatic>(true), make_wrapped_measure<IsRefcall>(true)},
    [] (const MeasureVector& measures) -> std::size_t {
        assert(measures.size() == 2);
        if (boost::get<bool>(measures.front())) {
            return 1;
        } else if (boost::get<bool>(measures.back())) {
            return 2;
        } else {
            return 0;
        }},
    output_config, threading, progress
} {}

SomaticThresholdVariantCallFilter::SomaticThresholdVariantCallFilter(FacetFactory facet_factory,
                                                                     ConditionVectorPair somatic,
                                                                     ConditionVectorPair reference,
                                                                     OutputOptions output_config,
                                                                     ConcurrencyPolicy threading,
                                                                     boost::optional<ProgressMeter&> progress)
: ConditionalThresholdVariantCallFilter {
    std::move(facet_factory),
    {{{{make_wrapped_measure<IsSomatic>(false), make_wrapped_threshold<EqualThreshold<bool>>(false)}}, {}}, std::move(somatic), std::move(reference)},
    {make_wrapped_measure<IsSomatic>(true), make_wrapped_measure<IsRefcall>(true)},
    [] (const MeasureVector& measures) -> std::size_t {
        assert(measures.size() == 2);
        if (boost::get<bool>(measures.front())) {
            return 1;
        } else if (boost::get<bool>(measures.back())) {
            return 2;
        } else {
            return 0;
        }},
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
