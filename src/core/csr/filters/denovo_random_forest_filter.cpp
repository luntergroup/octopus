// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_random_forest_filter.hpp"

#include <utility>

#include "../measures/is_denovo.hpp"
#include "../measures/is_refcall.hpp"

namespace octopus { namespace csr {

DeNovoRandomForestVariantCallFilter::DeNovoRandomForestVariantCallFilter(FacetFactory facet_factory,
                                                                         std::vector<MeasureWrapper> measures,
                                                                         Path germline_forest, Path denovo_forest,
                                                                         OutputOptions output_config,
                                                                         ConcurrencyPolicy threading,
                                                                         Path temp_directory,
                                                                         boost::optional<ProgressMeter&> progress)
: ConditionalRandomForestFilter {
    std::move(facet_factory),
    std::move(measures),
    {make_wrapped_measure<IsDenovo>(true)},
    [] (const MeasureVector& measures) -> std::int8_t { return !boost::get<bool>(measures.front()); },
    {std::move(germline_forest), std::move(denovo_forest)},
    std::move(output_config),
    std::move(threading),
    std::move(temp_directory),
    progress
} {}

DeNovoRandomForestVariantCallFilter::DeNovoRandomForestVariantCallFilter(FacetFactory facet_factory,
                                                                         std::vector<MeasureWrapper> measures,
                                                                         Path germline_forest, Path denovo_forest,
                                                                         Phred<double> min_forest_quality,
                                                                         OutputOptions output_config,
                                                                         ConcurrencyPolicy threading,
                                                                         Path temp_directory,
                                                                         boost::optional<ProgressMeter&> progress)
: ConditionalRandomForestFilter {
    std::move(facet_factory),
    std::move(measures),
    {make_wrapped_measure<IsDenovo>(true)},
    [] (const MeasureVector& measures) -> std::int8_t { return !boost::get<bool>(measures.front()); },
    {std::move(germline_forest), std::move(denovo_forest)},
    min_forest_quality,
    std::move(output_config),
    std::move(threading),
    std::move(temp_directory),
    progress
} {}

DeNovoRandomForestVariantCallFilter::DeNovoRandomForestVariantCallFilter(FacetFactory facet_factory,
                                                                         std::vector<MeasureWrapper> measures,
                                                                         Path denovo_forest,
                                                                         OutputOptions output_config,
                                                                         ConcurrencyPolicy threading,
                                                                         Path temp_directory,
                                                                         boost::optional<ProgressMeter&> progress)
: ConditionalRandomForestFilter {
    std::move(facet_factory),
    std::move(measures),
    {make_wrapped_measure<IsDenovo>(false)},
    [] (const MeasureVector& measures) -> std::int8_t { return !boost::get<bool>(measures.front()); },
    {std::move(denovo_forest)},
    std::move(output_config),
    std::move(threading),
    std::move(temp_directory),
    progress
} {}

DeNovoRandomForestVariantCallFilter::DeNovoRandomForestVariantCallFilter(FacetFactory facet_factory,
                                                                         std::vector<MeasureWrapper> measures,
                                                                         Path denovo_forest,
                                                                         Phred<double> min_forest_quality,
                                                                         OutputOptions output_config,
                                                                         ConcurrencyPolicy threading,
                                                                         Path temp_directory,
                                                                         boost::optional<ProgressMeter&> progress)
: ConditionalRandomForestFilter {
    std::move(facet_factory),
    std::move(measures),
    {make_wrapped_measure<IsDenovo>(false)},
    [] (const MeasureVector& measures) -> std::int8_t { return !boost::get<bool>(measures.front()); },
    {std::move(denovo_forest)},
    min_forest_quality,
    std::move(output_config),
    std::move(threading),
    std::move(temp_directory),
    progress
} {}

bool DeNovoRandomForestVariantCallFilter::is_soft_filtered(const ClassificationList& sample_classifications,
                                                           const MeasureVector& measures) const
{
    return std::any_of(std::cbegin(sample_classifications), std::cend(sample_classifications),
                       [] (const auto& c) { return c.category != Classification::Category::unfiltered; });
}

} // namespace csr
} // namespace octopus
