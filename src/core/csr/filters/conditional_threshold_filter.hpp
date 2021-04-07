// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef conditional_threshold_filter_hpp
#define conditional_threshold_filter_hpp

#include <vector>
#include <string>
#include <functional>
#include <cstddef>

#include <boost/optional.hpp>

#include "threshold_filter.hpp"
#include "logging/progress_meter.hpp"
#include "../facets/facet_factory.hpp"
#include "../measures/measure.hpp"

namespace octopus { namespace csr {

class ConditionalThresholdVariantCallFilter : public ThresholdVariantCallFilter
{
public:
    ConditionalThresholdVariantCallFilter() = delete;
    
    ConditionalThresholdVariantCallFilter(FacetFactory facet_factory,
                                          std::vector<ConditionVectorPair> conditions,
                                          std::vector<MeasureWrapper> chooser_measures,
                                          std::function<std::size_t(std::vector<Measure::ResultType>)> chooser,
                                          OutputOptions output_config,
                                          ConcurrencyPolicy threading,
                                          boost::optional<ProgressMeter&> progress = boost::none);
    
    ConditionalThresholdVariantCallFilter(const ConditionalThresholdVariantCallFilter&)            = delete;
    ConditionalThresholdVariantCallFilter& operator=(const ConditionalThresholdVariantCallFilter&) = delete;
    ConditionalThresholdVariantCallFilter(ConditionalThresholdVariantCallFilter&&)                 = delete;
    ConditionalThresholdVariantCallFilter& operator=(ConditionalThresholdVariantCallFilter&&)      = delete;
    
    virtual ~ConditionalThresholdVariantCallFilter() override = default;

private:
    struct MeasureIndexRange
    {
        std::size_t measure_begin, measure_end, threshold_begin;
    };
    
    std::vector<MeasureIndexRange> hard_ranges_, soft_ranges_;
    std::function<std::size_t(std::vector<Measure::ResultType>)> chooser_;
    std::vector<bool> unique_filter_keys_;
    std::size_t first_chooser_measure_index_, num_chooser_measures_;
    
    virtual bool passes_all_hard_filters(const MeasureVector& measures) const override;
    virtual bool passes_all_soft_filters(const MeasureVector& measures) const override;
    virtual std::vector<std::string> get_failing_vcf_filter_keys(const MeasureVector& measures) const override;
    
    std::size_t choose_filter(const MeasureVector& measures) const;
    bool passes_all_hard_filters(const MeasureVector& measures, MeasureIndexRange range) const;
    bool passes_all_soft_filters(const MeasureVector& measures, MeasureIndexRange range) const;
    std::vector<std::string> get_failing_vcf_filter_keys(const MeasureVector& measures, MeasureIndexRange range) const;
};

} // namespace csr
} // namespace octopus

#endif
