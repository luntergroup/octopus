// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef threshold_filter_hpp
#define threshold_filter_hpp

#include <boost/optional.hpp>

#include "variant_call_filter.hpp"
#include "logging/progress_meter.hpp"
#include "../facets/facet_factory.hpp"

namespace octopus {

class VcfHeader;

namespace csr {

class ThresholdVariantCallFilter : public VariantCallFilter
{
public:
    using VariantCallFilter::OutputOptions;
    
    struct Threshold
    {
        virtual bool operator()(Measure::ResultType value) const noexcept { return true; }
        virtual ~Threshold() = default;
    };
    
    ThresholdVariantCallFilter() = delete;
    
    ThresholdVariantCallFilter(FacetFactory facet_factory,
                               std::vector<MeasureWrapper> measures,
                               std::vector<std::unique_ptr<Threshold>> thresholds,
                               OutputOptions output_config,
                               boost::optional<ProgressMeter&> progress = boost::none);
    
    ThresholdVariantCallFilter(const ThresholdVariantCallFilter&)            = delete;
    ThresholdVariantCallFilter& operator=(const ThresholdVariantCallFilter&) = delete;
    ThresholdVariantCallFilter(ThresholdVariantCallFilter&&)                 = default;
    ThresholdVariantCallFilter& operator=(ThresholdVariantCallFilter&&)      = default;
    
    virtual ~ThresholdVariantCallFilter() = default;

private:
    std::vector<std::unique_ptr<Threshold>> thresholds_;
    
    virtual void annotate(VcfHeader& dest) const override;
    virtual Classification classify(const MeasureVector& measures) const override;
    
    bool passes_all_filters(const MeasureVector& measures) const;
};

} // namespace csr
} // namespace octopus

#endif /* threshold_filter_hpp */
