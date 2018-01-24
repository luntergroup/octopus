// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "training_filter_factory.hpp"

#include <iterator>
#include <algorithm>

#include "passing_filter.hpp"
#include "../measures/measure_factory.hpp"

namespace octopus { namespace csr {

namespace {

std::vector<MeasureWrapper> parse_measures(const std::set<std::string>& measure_names)
{
    std::vector<MeasureWrapper> result {};
    result.reserve(measure_names.size());
    std::transform(std::cbegin(measure_names), std::cend(measure_names), std::back_inserter(result), make_measure);
    return result;
}

} // namespace

TrainingFilterFactory::TrainingFilterFactory(const std::set<std::string>& measure_names)
: measures_ {parse_measures(measure_names)}
{}

std::unique_ptr<VariantCallFilterFactory> TrainingFilterFactory::do_clone() const
{
    return std::make_unique<TrainingFilterFactory>(*this);
}

std::unique_ptr<VariantCallFilter> TrainingFilterFactory::do_make(FacetFactory facet_factory,
                                                                  VariantCallFilter::OutputOptions output_config,
                                                                  boost::optional<ProgressMeter&> progress,
                                                                  VariantCallFilter::ConcurrencyPolicy threading) const
{
    output_config.annotate_measures = true;
    output_config.clear_info = true;
    output_config.clear_existing_filters = true;
    return std::make_unique<PassingVariantCallFilter>(std::move(facet_factory), std::move(measures_),
                                                      output_config, threading, progress);
}

} // namespace csr
} // namespace octopus
