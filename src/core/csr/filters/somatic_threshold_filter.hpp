// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef somatic_threshold_filter_hpp
#define somatic_threshold_filter_hpp

#include <vector>
#include <string>

#include <boost/optional.hpp>

#include "threshold_filter.hpp"
#include "logging/progress_meter.hpp"
#include "../facets/facet_factory.hpp"
#include "../measures/measure.hpp"

namespace octopus { namespace csr {

class SomaticThresholdVariantCallFilter : public ThresholdVariantCallFilter
{
public:
    SomaticThresholdVariantCallFilter() = delete;
    
    SomaticThresholdVariantCallFilter(FacetFactory facet_factory,
                                      std::vector<Condition> germline_hard_conditions,
                                      std::vector<Condition> germline_soft_conditions,
                                      std::vector<Condition> somatic_hard_conditions,
                                      std::vector<Condition> somatic_soft_conditions,
                                      OutputOptions output_config,
                                      ConcurrencyPolicy threading,
                                      boost::optional<ProgressMeter&> progress = boost::none);
    
    SomaticThresholdVariantCallFilter(const SomaticThresholdVariantCallFilter&)            = delete;
    SomaticThresholdVariantCallFilter& operator=(const SomaticThresholdVariantCallFilter&) = delete;
    SomaticThresholdVariantCallFilter(SomaticThresholdVariantCallFilter&&)                 = default;
    SomaticThresholdVariantCallFilter& operator=(SomaticThresholdVariantCallFilter&&)      = default;
    
    virtual ~SomaticThresholdVariantCallFilter() override = default;

private:
    std::size_t num_germline_hard_conditions_, num_germline_soft_conditions_;
    bool all_unique_germline_filter_keys_, all_unique_somatic_filter_keys_;
    
    virtual bool passes_all_hard_filters(const MeasureVector& measures) const override;
    virtual bool passes_all_soft_filters(const MeasureVector& measures) const override;
    virtual std::vector<std::string> get_failing_vcf_filter_keys(const MeasureVector& measures) const override;
    
    bool passes_all_germline_hard_filters(const MeasureVector& measures) const;
    bool passes_all_somatic_hard_filters(const MeasureVector& measures) const;
    bool passes_all_germline_soft_filters(const MeasureVector& measures) const;
    bool passes_all_somatic_soft_filters(const MeasureVector& measures) const;
    std::vector<std::string> get_failing_germline_vcf_filter_keys(const MeasureVector& measures) const;
    std::vector<std::string> get_failing_somatic_vcf_filter_keys(const MeasureVector& measures) const;
};

} // namespace csr
} // namespace octopus

#endif
