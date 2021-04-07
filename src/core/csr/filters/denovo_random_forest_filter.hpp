// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef denovo_random_forest_filter_hpp
#define denovo_random_forest_filter_hpp

#include <vector>
#include <string>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "threshold_filter.hpp"
#include "random_forest_filter.hpp"
#include "logging/progress_meter.hpp"
#include "../facets/facet_factory.hpp"
#include "../measures/measure.hpp"

namespace octopus { namespace csr {

class DeNovoRandomForestVariantCallFilter : public RandomForestFilter
{
public:
    using RandomForestFilter::Options;
    
    DeNovoRandomForestVariantCallFilter() = delete;
    
    DeNovoRandomForestVariantCallFilter(FacetFactory facet_factory,
                                        Path germline_forest, Path denovo_forest,
                                        OutputOptions output_config,
                                        ConcurrencyPolicy threading,
                                        Path temp_directory,
                                        Options options,
                                        boost::optional<ProgressMeter&> progress = boost::none);
    // De novo only
    DeNovoRandomForestVariantCallFilter(FacetFactory facet_factory,
                                        Path denovo_forest,
                                        OutputOptions output_config,
                                        ConcurrencyPolicy threading,
                                        Path temp_directory,
                                        Options options,
                                        boost::optional<ProgressMeter&> progress = boost::none);
    
    DeNovoRandomForestVariantCallFilter(const DeNovoRandomForestVariantCallFilter&)            = delete;
    DeNovoRandomForestVariantCallFilter& operator=(const DeNovoRandomForestVariantCallFilter&) = delete;
    DeNovoRandomForestVariantCallFilter(DeNovoRandomForestVariantCallFilter&&)                 = delete;
    DeNovoRandomForestVariantCallFilter& operator=(DeNovoRandomForestVariantCallFilter&&)      = delete;
    
    virtual ~DeNovoRandomForestVariantCallFilter() override = default;
};
    
} // namespace csr
} // namespace octopus

#endif
