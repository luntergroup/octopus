// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_random_forest_filter.hpp"

#include <utility>

#include "../measures/is_denovo.hpp"
#include "../measures/is_refcall.hpp"

namespace octopus { namespace csr {

DeNovoRandomForestVariantCallFilter::DeNovoRandomForestVariantCallFilter(FacetFactory facet_factory,
                                                                         Path germline_forest, Path denovo_forest,
                                                                         OutputOptions output_config,
                                                                         ConcurrencyPolicy threading,
                                                                         Path temp_directory,
                                                                         Options options,
                                                                         boost::optional<ProgressMeter&> progress)
: RandomForestFilter {
    std::move(facet_factory),
    {make_wrapped_measure<IsDenovo>(true)},
    [] (const MeasureVector& measures) -> std::int8_t { return !get_value_type<bool>(measures.front()); },
    {std::move(germline_forest), std::move(denovo_forest)},
    std::move(output_config),
    std::move(threading),
    std::move(temp_directory),
    std::move(options),
    progress
} {}

DeNovoRandomForestVariantCallFilter::DeNovoRandomForestVariantCallFilter(FacetFactory facet_factory,
                                                                         Path denovo_forest,
                                                                         OutputOptions output_config,
                                                                         ConcurrencyPolicy threading,
                                                                         Path temp_directory,
                                                                         Options options,
                                                                         boost::optional<ProgressMeter&> progress)
: RandomForestFilter {
    std::move(facet_factory),
    {make_wrapped_measure<IsDenovo>(false)},
    [] (const MeasureVector& measures) -> std::int8_t { return !get_value_type<bool>(measures.front()); },
    {std::move(denovo_forest)},
    std::move(output_config),
    std::move(threading),
    std::move(temp_directory),
    std::move(options),
    progress
} {}

} // namespace csr
} // namespace octopus
