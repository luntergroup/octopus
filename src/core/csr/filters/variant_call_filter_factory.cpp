// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_call_filter_factory.hpp"

#include <boost/variant.hpp>

#include "threshold_filter.hpp"
#include "core/csr/facets/facet_factory.hpp"
#include "core/csr/measures/qual.hpp"
#include "core/csr/measures/mean_mapping_quality.hpp"
#include "core/csr/measures/allele_frequency.hpp"
#include "core/csr/measures/strand_bias.hpp"

namespace octopus { namespace csr {

struct LessThreshold : public ThresholdVariantCallFilter::Threshold
{
    explicit LessThreshold(Measure::ResultType target) : target_ {target} {}
    bool operator()(Measure::ResultType value) const noexcept
    {
        return value >= target_;
    }
private:
    Measure::ResultType target_;
};

std::unique_ptr<VariantCallFilter> VariantCallFilterFactory::make(const ReferenceGenome& reference,
                                                                  BufferedReadPipe read_pipe,
                                                                  boost::optional<ProgressMeter&> progress) const
{
    FacetFactory facet_factory {reference, std::move(read_pipe)};
    std::vector<MeasureWrapper> measures {};
    measures.push_back(make_wrapped_measure<Qual>());
    measures.push_back(make_wrapped_measure<MeanMappingQuality>());
    measures.push_back(make_wrapped_measure<AlleleFrequency>());
    measures.push_back(make_wrapped_measure<StrandBias>());
    using TP = std::unique_ptr<ThresholdVariantCallFilter::Threshold>;
    std::vector<TP> thresholds {};
    thresholds.push_back(std::make_unique<LessThreshold>(20.0));
    thresholds.push_back(std::make_unique<LessThreshold>(40.0));
    thresholds.push_back(std::make_unique<LessThreshold>(0.2));
    thresholds.push_back(std::make_unique<LessThreshold>(0.1));
    return std::make_unique<ThresholdVariantCallFilter>(std::move(facet_factory), std::move(measures), std::move(thresholds));
    return std::make_unique<ThresholdVariantCallFilter>(std::move(facet_factory), std::move(measures),
                                                        std::move(thresholds), progress);
}

} // namespace csr
} // namespace octopus
