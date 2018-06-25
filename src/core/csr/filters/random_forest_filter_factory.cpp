// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "random_forest_filter_factory.hpp"

#include "../measures/measure_factory.hpp"
#include "somatic_random_forest_filter.hpp"

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
: ranger_forests_ {std::move(ranger_forest)}
, temp_directory_ {std::move(temp_directory)}
{
    measures_ = parse_measures({"AF", "CRF", "DP", "FRF", "GC", "GQ", "MQ0", "MQ", "BQ", "QUAL", "QD", "SB", "ARF"});
}

RandomForestFilterFactory::RandomForestFilterFactory(Path germline_ranger_forest, Path somatic_ranger_forest,
                                                     Path refcall_ranger_forest, Path temp_directory)
: ranger_forests_ {std::move(germline_ranger_forest), std::move(somatic_ranger_forest), std::move(refcall_ranger_forest)}
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
    if (ranger_forests_.size() == 1) {
        return std::make_unique<RandomForestFilter>(std::move(facet_factory), measures_, output_config, threading,
                                                    ranger_forests_[0], temp_directory_, progress);
    } else {
        assert(ranger_forests_.size() == 3);
        return std::make_unique<SomaticRandomForestVariantCallFilter>(std::move(facet_factory), measures_,
                                                                      ranger_forests_[0], ranger_forests_[1], ranger_forests_[2],
                                                                      output_config, threading, temp_directory_, progress);
    }
}

} // namespace csr
} // namespace octopus
