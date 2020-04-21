// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_random_forest_filter.hpp"

#include <utility>

#include "../measures/is_somatic.hpp"
#include "../measures/is_refcall.hpp"

namespace octopus { namespace csr {

SomaticRandomForestVariantCallFilter::Options::Options(RandomForestFilter::Options common, bool use_somatic_forest_for_refcalls)
: RandomForestFilter::Options {common}
, use_somatic_forest_for_refcalls {use_somatic_forest_for_refcalls}
{}

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
    [use_somatic_forest_for_refcalls = options.use_somatic_forest_for_refcalls] (const MeasureVector& measures) -> std::int8_t {
        assert(measures.size() == 2);
        if (boost::get<bool>(measures.front())) {
            return 1;
        } else if (boost::get<bool>(measures.back())) {
            return use_somatic_forest_for_refcalls;
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
