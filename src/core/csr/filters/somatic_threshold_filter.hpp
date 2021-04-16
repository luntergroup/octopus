// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef somatic_threshold_filter_hpp
#define somatic_threshold_filter_hpp

#include <vector>
#include <string>

#include <boost/optional.hpp>

#include "threshold_filter.hpp"
#include "conditional_threshold_filter.hpp"
#include "logging/progress_meter.hpp"
#include "../facets/facet_factory.hpp"
#include "../measures/measure.hpp"

namespace octopus { namespace csr {

class SomaticThresholdVariantCallFilter : public ConditionalThresholdVariantCallFilter
{
public:
    SomaticThresholdVariantCallFilter() = delete;
    
    SomaticThresholdVariantCallFilter(FacetFactory facet_factory,
                                      ConditionVectorPair germline,
                                      ConditionVectorPair somatic,
                                      ConditionVectorPair reference,
                                      OutputOptions output_config,
                                      ConcurrencyPolicy threading,
                                      boost::optional<ProgressMeter&> progress = boost::none);
    // Hard filter germline
    SomaticThresholdVariantCallFilter(FacetFactory facet_factory,
                                      ConditionVectorPair somatic,
                                      ConditionVectorPair reference,
                                      OutputOptions output_config,
                                      ConcurrencyPolicy threading,
                                      boost::optional<ProgressMeter&> progress = boost::none);
    
    SomaticThresholdVariantCallFilter(const SomaticThresholdVariantCallFilter&)            = delete;
    SomaticThresholdVariantCallFilter& operator=(const SomaticThresholdVariantCallFilter&) = delete;
    SomaticThresholdVariantCallFilter(SomaticThresholdVariantCallFilter&&)                 = delete;
    SomaticThresholdVariantCallFilter& operator=(SomaticThresholdVariantCallFilter&&)      = delete;
    
    virtual ~SomaticThresholdVariantCallFilter() override = default;

private:
    virtual bool is_soft_filtered(const ClassificationList& sample_classifications, boost::optional<Phred<double>> joint_quality,
                                  const MeasureVector& measures, std::vector<std::string>& reasons) const override;
};

} // namespace csr
} // namespace octopus

#endif
