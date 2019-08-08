// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_threshold_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <cassert>

#include <boost/variant.hpp>

#include "../measures/is_denovo.hpp"
#include "../measures/is_refcall.hpp"

namespace octopus { namespace csr {

DeNovoThresholdVariantCallFilter::DeNovoThresholdVariantCallFilter(FacetFactory facet_factory,
                                                                   ConditionVectorPair germline,
                                                                   ConditionVectorPair denovo,
                                                                   ConditionVectorPair reference,
                                                                   OutputOptions output_config,
                                                                   ConcurrencyPolicy threading,
                                                                   boost::optional<ProgressMeter&> progress)
: ConditionalThresholdVariantCallFilter {
    std::move(facet_factory),
    {std::move(germline), std::move(denovo), std::move(reference)},
    {make_wrapped_measure<IsRefcall>(true), make_wrapped_measure<IsDenovo>(false)},
    [] (const MeasureVector& measures) -> std::size_t {
        assert(measures.size() == 2);
        if (boost::get<bool>(measures[0])) {
            return 2; // REFCALL sample
        } else if (boost::get<bool>(measures[1])) {
            return 1; // DENOVO call
        } else {
            return 0; // germline variant sample
        }},
    output_config, threading, progress
} {}

DeNovoThresholdVariantCallFilter::DeNovoThresholdVariantCallFilter(FacetFactory facet_factory,
                                                                   ConditionVectorPair denovo,
                                                                   ConditionVectorPair reference,
                                                                   OutputOptions output_config,
                                                                   ConcurrencyPolicy threading,
                                                                   boost::optional<ProgressMeter&> progress)
: ConditionalThresholdVariantCallFilter {
    std::move(facet_factory),
    {{{{make_wrapped_measure<IsDenovo>(false), make_wrapped_threshold<EqualThreshold<bool>>(false)}}, {}}, std::move(denovo), std::move(reference)},
    {make_wrapped_measure<IsRefcall>(true), make_wrapped_measure<IsDenovo>(true), make_wrapped_measure<IsDenovo>(false)},
    [] (const MeasureVector& measures) -> std::size_t {
        assert(measures.size() == 3);
        if (!boost::get<bool>(measures[2])) {
            return 0; // Not DENOVO call
        } else if (boost::get<bool>(measures[1])) {
            return 1; // DENOVO sample
        } else if (boost::get<bool>(measures[0])) {
            return 2; // DENOVO call REFCALL parent
        } else {
            return 1; // DENOVO call non-REFCALL parent
        }},
    output_config, threading, progress
} {}

bool DeNovoThresholdVariantCallFilter::is_soft_filtered(const ClassificationList& sample_classifications,
                                                        const boost::optional<Phred<double>> joint_quality,
                                                        const MeasureVector& measures,
                                                        std::vector<std::string>& reasons) const
{
    if (boost::get<bool>(measures.back())) {
        return std::any_of(std::cbegin(sample_classifications), std::cend(sample_classifications),
                           [] (const auto& c) { return c.category != Classification::Category::unfiltered; });
    } else {
        return std::all_of(std::cbegin(sample_classifications), std::cend(sample_classifications),
                           [] (const auto& c) { return c.category != Classification::Category::unfiltered; });
    }
}

} // namespace csr
} // namespace octopus
