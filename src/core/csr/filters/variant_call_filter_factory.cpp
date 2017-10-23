// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_call_filter_factory.hpp"

#include <boost/variant.hpp>

#include "threshold_filter.hpp"
#include "core/csr/facets/facet_factory.hpp"
#include "core/csr/measures/qual.hpp"
#include "core/csr/measures/mean_mapping_quality.hpp"
#include "core/csr/measures/model_posterior.hpp"
#include "core/csr/measures/allele_frequency.hpp"
#include "core/csr/measures/strand_bias.hpp"
#include "core/csr/measures/mapping_quality_divergence.hpp"

namespace octopus { namespace csr {

struct LessThreshold : public ThresholdVariantCallFilter::Threshold
{
    explicit LessThreshold(double target) : visitor_ {target} {}
    bool operator()(Measure::ResultType value) const noexcept
    {
        return boost::apply_visitor(visitor_, value);
    }
private:
    struct LessVisitor : public boost::static_visitor<bool>
    {
        explicit LessVisitor(double target) : target {target} {}
        bool operator()(double value) const noexcept { return value >= target; }
        bool operator()(boost::optional<double> value) const noexcept
        {
            return !value || *value >= target;
        }
        double target;
    };
    LessVisitor visitor_;
};

struct GreaterThreshold : public ThresholdVariantCallFilter::Threshold
{
    explicit GreaterThreshold(double target) : visitor_ {target} {}
    bool operator()(Measure::ResultType value) const noexcept
    {
        return boost::apply_visitor(visitor_, value);
    }
private:
    struct GreaterVisitor : public boost::static_visitor<bool>
    {
        explicit GreaterVisitor(double target) : target {target} {}
        bool operator()(double value) const noexcept { return value <= target; }
        bool operator()(boost::optional<double> value) const noexcept
        {
            return !value || *value <= target;
        }
        double target;
    };
    GreaterVisitor visitor_;
};

std::unique_ptr<VariantCallFilter> VariantCallFilterFactory::make(const ReferenceGenome& reference,
                                                                  BufferedReadPipe read_pipe,
                                                                  boost::optional<ProgressMeter&> progress) const
{
    FacetFactory facet_factory {reference, std::move(read_pipe)};
    std::vector<MeasureWrapper> measures {};
    measures.push_back(make_wrapped_measure<Qual>());
    measures.push_back(make_wrapped_measure<MeanMappingQuality>());
    measures.push_back(make_wrapped_measure<ModelPosterior>());
    //measures.push_back(make_wrapped_measure<AlleleFrequency>());
    //measures.push_back(make_wrapped_measure<StrandBias>());
    //measures.push_back(make_wrapped_measure<MappingQualityDivergence>());
    using TP = std::unique_ptr<ThresholdVariantCallFilter::Threshold>;
    std::vector<TP> thresholds {};
    thresholds.push_back(std::make_unique<LessThreshold>(10.0));
    thresholds.push_back(std::make_unique<LessThreshold>(30.0));
    thresholds.push_back(std::make_unique<LessThreshold>(20.0));
    //thresholds.push_back(std::make_unique<LessThreshold>(0.2));
    //thresholds.push_back(std::make_unique<LessThreshold>(0.1));
    //thresholds.push_back(std::make_unique<GreaterThreshold>(20.0));
    return std::make_unique<ThresholdVariantCallFilter>(std::move(facet_factory), std::move(measures),
                                                        std::move(thresholds), progress);
}

} // namespace csr
} // namespace octopus
