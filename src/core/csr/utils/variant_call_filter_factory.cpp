// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_call_filter_factory.hpp"

#include "../filters/threshold_filter.hpp"

#include "../measures/qual.hpp"
#include "../measures/mean_mapping_quality.hpp"

namespace octopus { namespace csr {

struct LessThreshold : public ThresholdVariantCallFilter::Threshold
{
    explicit LessThreshold(double target) : target_ {target} {}
    bool operator()(double value) const noexcept
    {
        return value >= target_;
    }
private:
    double target_;
};

std::unique_ptr<VariantCallFilter> VariantCallFilterFactory::make(const ReferenceGenome& reference, const ReadPipe& read_pipe) const
{
    using TP = std::unique_ptr<ThresholdVariantCallFilter::Threshold>;
    std::vector<MeasureWrapper> measures {};
    measures.push_back(make_wrapped_measure<Qual>());
    measures.push_back(make_wrapped_measure<MeanMappingQuality>());
    std::vector<TP> thresholds {};
    thresholds.push_back(std::make_unique<LessThreshold>(20.0));
    thresholds.push_back(std::make_unique<LessThreshold>(40.0));
    return std::make_unique<ThresholdVariantCallFilter>(reference, read_pipe, std::move(measures), std::move(thresholds));
}

} // namespace csr
} // namespace octopus
