// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_random_forest_filter.hpp"

#include <utility>

#include "../measures/is_somatic.hpp"
#include "../measures/is_refcall.hpp"

namespace octopus { namespace csr {

SomaticRandomForestVariantCallFilter::SomaticRandomForestVariantCallFilter(FacetFactory facet_factory,
                                                                           Path germline_forest, Path somatic_forest,
                                                                           OutputOptions output_config,
                                                                           ConcurrencyPolicy threading,
                                                                           Path temp_directory,
                                                                           Options options,
                                                                           boost::optional<ProgressMeter&> progress)
: RandomForestFilter {
    std::move(facet_factory),
    {make_wrapped_measure<IsSomatic>(true), make_wrapped_measure<IsRefcall>(true)},
    [] (const MeasureVector& measures) -> std::int8_t {
        assert(measures.size() == 2);
        if (boost::get<bool>(measures.front())) {
            return 1;
        } else if (boost::get<bool>(measures.back())) {
            return 1;
        } else {
            return 0;
        }},
    {std::move(germline_forest), std::move(somatic_forest)},
    std::move(output_config),
    std::move(threading),
    std::move(temp_directory),
    std::move(options),
    progress
} {}

SomaticRandomForestVariantCallFilter::SomaticRandomForestVariantCallFilter(FacetFactory facet_factory,
                                                                           Path somatic_forest,
                                                                           OutputOptions output_config,
                                                                           ConcurrencyPolicy threading,
                                                                           Path temp_directory,
                                                                           Options options,
                                                                           boost::optional<ProgressMeter&> progress)
: RandomForestFilter {
    std::move(facet_factory),
    {make_wrapped_measure<IsSomatic>(false)},
    [] (const MeasureVector& measures) -> std::int8_t {
        assert(measures.size() == 1);
        return !boost::get<bool>(measures.front());
        },
    {std::move(somatic_forest)},
    std::move(output_config),
    std::move(threading),
    std::move(temp_directory),
    std::move(options),
    progress
} {}

} // namespace csr
} // namespace octopus
