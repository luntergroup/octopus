// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef somatic_random_forest_filter_hpp
#define somatic_random_forest_filter_hpp

#include <vector>
#include <string>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "threshold_filter.hpp"
#include "conditional_random_forest_filter.hpp"
#include "logging/progress_meter.hpp"
#include "../facets/facet_factory.hpp"
#include "../measures/measure.hpp"

namespace octopus { namespace csr {

class SomaticRandomForestVariantCallFilter : public ConditionalRandomForestFilter
{
public:
    SomaticRandomForestVariantCallFilter() = delete;
    
    SomaticRandomForestVariantCallFilter(FacetFactory facet_factory,
                                         std::vector<MeasureWrapper> measures,
                                         Path germline_forest, Path somatic_forest,
                                         OutputOptions output_config,
                                         ConcurrencyPolicy threading,
                                         Path temp_directory = "/tmp",
                                         boost::optional<ProgressMeter&> progress = boost::none);
    // Somatics only
    SomaticRandomForestVariantCallFilter(FacetFactory facet_factory,
                                         std::vector<MeasureWrapper> measures,
                                         Path somatic_forest,
                                         OutputOptions output_config,
                                         ConcurrencyPolicy threading,
                                         Path temp_directory = "/tmp",
                                         boost::optional<ProgressMeter&> progress = boost::none);
    
    SomaticRandomForestVariantCallFilter(const SomaticRandomForestVariantCallFilter&)            = delete;
    SomaticRandomForestVariantCallFilter& operator=(const SomaticRandomForestVariantCallFilter&) = delete;
    SomaticRandomForestVariantCallFilter(SomaticRandomForestVariantCallFilter&&)                 = default;
    SomaticRandomForestVariantCallFilter& operator=(SomaticRandomForestVariantCallFilter&&)      = default;
    
    virtual ~SomaticRandomForestVariantCallFilter() override = default;

private:
    virtual bool is_soft_filtered(const ClassificationList& sample_classifications, const MeasureVector& measures) const override;
};

} // namespace csr
} // namespace octopus

#endif
