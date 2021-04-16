// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef denovo_threshold_filter_hpp
#define denovo_threshold_filter_hpp

#include <vector>
#include <string>

#include <boost/optional.hpp>

#include "threshold_filter.hpp"
#include "conditional_threshold_filter.hpp"
#include "logging/progress_meter.hpp"
#include "../facets/facet_factory.hpp"
#include "../measures/measure.hpp"

namespace octopus { namespace csr {

class DeNovoThresholdVariantCallFilter : public ConditionalThresholdVariantCallFilter
{
public:
    DeNovoThresholdVariantCallFilter() = delete;
    
    DeNovoThresholdVariantCallFilter(FacetFactory facet_factory,
                                     ConditionVectorPair germline,
                                     ConditionVectorPair denovo,
                                     ConditionVectorPair reference,
                                     OutputOptions output_config,
                                     ConcurrencyPolicy threading,
                                     boost::optional<ProgressMeter&> progress = boost::none);
    
    // Hard filter germline
    DeNovoThresholdVariantCallFilter(FacetFactory facet_factory,
                                     ConditionVectorPair denovo,
                                     ConditionVectorPair reference,
                                     OutputOptions output_config,
                                     ConcurrencyPolicy threading,
                                     boost::optional<ProgressMeter&> progress = boost::none);
    
    DeNovoThresholdVariantCallFilter(const DeNovoThresholdVariantCallFilter&)            = delete;
    DeNovoThresholdVariantCallFilter& operator=(const DeNovoThresholdVariantCallFilter&) = delete;
    DeNovoThresholdVariantCallFilter(DeNovoThresholdVariantCallFilter&&)                 = delete;
    DeNovoThresholdVariantCallFilter& operator=(DeNovoThresholdVariantCallFilter&&)      = delete;
    
    virtual ~DeNovoThresholdVariantCallFilter() override = default;

private:
    virtual bool is_soft_filtered(const ClassificationList& sample_classifications, boost::optional<Phred<double>> joint_quality,
                                  const MeasureVector& measures, std::vector<std::string>& reasons) const override;
};

} // namespace csr
} // namespace octopus

#endif
