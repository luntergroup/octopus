// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "random_forest_filter_factory.hpp"

#include "../measures/measure_factory.hpp"

namespace octopus { namespace csr {

namespace {

std::vector<MeasureWrapper> parse_measures(const std::vector<std::string>& measure_names)
{
    std::vector<MeasureWrapper> result{};
    result.reserve(measure_names.size());
    std::transform(std::cbegin(measure_names), std::cend(measure_names), std::back_inserter(result), make_measure);
    return result;
}

} // namespace

RandomForestFilterFactory::RandomForestFilterFactory(Path ranger_forest, Path temp_directory)
: ranger_forest_ {std::move(ranger_forest)}
, temp_directory_ {std::move(temp_directory)}
{
    measures_ = parse_measures({"AF", "CRF", "DP", "FRF", "GC", "GQ", "MQ0", "MQ", "BQ", "QUAL", "QD", "SB", "ARF"});
}

std::unique_ptr<VariantCallFilterFactory> RandomForestFilterFactory::do_clone() const
{
    return std::make_unique<RandomForestFilterFactory>(*this);
}

std::unique_ptr<VariantCallFilter> RandomForestFilterFactory::do_make(FacetFactory facet_factory,
                                                                      VariantCallFilter::OutputOptions output_config,
                                                                      boost::optional<ProgressMeter&> progress,
                                                                      VariantCallFilter::ConcurrencyPolicy threading) const
{
    return std::make_unique<RandomForestFilter>(std::move(facet_factory), measures_, output_config, threading,
                                                ranger_forest_, temp_directory_, progress);
}

} // namespace csr
} // namespace octopus
