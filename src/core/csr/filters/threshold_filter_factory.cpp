// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "threshold_filter_factory.hpp"

#include "core/csr/facets/facet_factory.hpp"
#include "core/csr/measures/qual.hpp"
#include "core/csr/measures/mean_mapping_quality.hpp"
#include "core/csr/measures/model_posterior.hpp"
#include "core/csr/measures/allele_frequency.hpp"
#include "core/csr/measures/strand_bias.hpp"
#include "core/csr/measures/mapping_quality_divergence.hpp"

namespace octopus { namespace csr {

ThresholdFilterFactory::ThresholdFilterFactory()
{
    using std::make_unique;
    conditions_.push_back({make_wrapped_measure<Qual>(), make_wrapped_threshold<LessThreshold>(10.0)});
    conditions_.push_back({make_wrapped_measure<MeanMappingQuality>(), make_wrapped_threshold<LessThreshold>(30.0)});
    conditions_.push_back({make_wrapped_measure<ModelPosterior>(), make_wrapped_threshold<LessThreshold>(20.0)});
    conditions_.push_back({make_wrapped_measure<AlleleFrequency>(), make_wrapped_threshold<LessThreshold>(0.2)});
    conditions_.push_back({make_wrapped_measure<StrandBias>(), make_wrapped_threshold<LessThreshold>(0.1)});
    conditions_.push_back({make_wrapped_measure<MappingQualityDivergence>(), make_wrapped_threshold<GreaterThreshold>(20.0)});
}

std::unique_ptr<VariantCallFilterFactory> ThresholdFilterFactory::do_clone() const
{
    return std::make_unique<ThresholdFilterFactory>(*this);
}

std::unique_ptr<VariantCallFilter> ThresholdFilterFactory::do_make(FacetFactory facet_factory,
                                                                   VariantCallFilter::OutputOptions output_config,
                                                                   boost::optional<ProgressMeter&> progress) const
{
    return std::make_unique<ThresholdVariantCallFilter>(std::move(facet_factory), conditions_, output_config, progress);
}

} // namespace csr
} // namespace octopus
