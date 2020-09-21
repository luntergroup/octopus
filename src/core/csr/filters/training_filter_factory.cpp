// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "training_filter_factory.hpp"

#include <iterator>
#include <algorithm>

#include "passing_filter.hpp"
#include "../measures/measure_factory.hpp"

namespace octopus { namespace csr {

namespace {

bool is_named_measure_set(const std::string& name)
{
    return name == "all";
}

std::vector<MeasureWrapper> make_measure_set(const std::string& name)
{
    if (name == "all") {
        return make_measures(get_all_measure_names());
    } else {
        return {};
    }
}

std::vector<MeasureWrapper> make_measures(const std::set<std::string>& measures)
{
    std::vector<MeasureWrapper> result {};
    result.reserve(measures.size());
    std::transform(std::cbegin(measures), std::cend(measures), std::back_inserter(result), make_measure);
    return result;
}

} // namespace

TrainingFilterFactory::TrainingFilterFactory(const std::set<std::string>& measures)
{
    if (measures.size() == 1 && is_named_measure_set(*std::cbegin(measures))) {
        measures_ = make_measure_set(*std::cbegin(measures));
    } else {
        measures_ = make_measures(measures);
    }
}

std::unique_ptr<VariantCallFilterFactory> TrainingFilterFactory::do_clone() const
{
    return std::make_unique<TrainingFilterFactory>(*this);
}

std::unique_ptr<VariantCallFilter> TrainingFilterFactory::do_make(FacetFactory facet_factory,
                                                                  VariantCallFilter::OutputOptions output_config,
                                                                  boost::optional<ProgressMeter&> progress,
                                                                  VariantCallFilter::ConcurrencyPolicy threading) const
{
    output_config.annotations.clear();
    output_config.annotations.reserve(measures_.size());
    for (const auto& measure : measures_) output_config.annotations.insert(measure.name());
    output_config.clear_info = true;
    output_config.clear_existing_filters = true;
    output_config.aggregate_allele_annotations = true;
    return std::make_unique<PassingVariantCallFilter>(std::move(facet_factory), std::move(measures_),
                                                      output_config, threading, progress);
}

} // namespace csr
} // namespace octopus
