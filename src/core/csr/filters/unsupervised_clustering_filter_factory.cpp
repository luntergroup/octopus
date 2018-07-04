// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "unsupervised_clustering_filter_factory.hpp"

#include <iterator>
#include <algorithm>

#include "../measures/measure_factory.hpp"

namespace octopus { namespace csr {

namespace {

std::vector<MeasureWrapper> parse_measures(const std::set<std::string>& measure_names)
{
    std::vector<MeasureWrapper> result{};
    result.reserve(measure_names.size());
    std::transform(std::cbegin(measure_names), std::cend(measure_names), std::back_inserter(result), make_measure);
    return result;
}

} // namespace

UnsupervisedClusteringFilterFactory::UnsupervisedClusteringFilterFactory(const std::set<std::string>& measure_names)
: measures_ {parse_measures(measure_names)}
{}

std::unique_ptr<VariantCallFilterFactory> UnsupervisedClusteringFilterFactory::do_clone() const
{
    return std::make_unique<UnsupervisedClusteringFilterFactory>(*this);
}

std::unique_ptr<VariantCallFilter> UnsupervisedClusteringFilterFactory::do_make(FacetFactory facet_factory,
                                                                                VariantCallFilter::OutputOptions output_config,
                                                                                boost::optional<ProgressMeter&> progress,
                                                                                VariantCallFilter::ConcurrencyPolicy threading) const
{
    return std::make_unique<UnsupervisedClusteringFilter>(std::move(facet_factory), measures_, output_config, threading, progress);
}

} // namespace csr
} // namespace octopus
