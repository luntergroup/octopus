// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef passing_filter_hpp
#define passing_filter_hpp

#include <vector>

#include <boost/optional.hpp>

#include "single_pass_variant_call_filter.hpp"
#include "logging/progress_meter.hpp"
#include "../facets/facet_factory.hpp"
#include "../measures/measure.hpp"

namespace octopus { namespace csr {

class PassingVariantCallFilter : public SinglePassVariantCallFilter
{
public:
    PassingVariantCallFilter() = delete;
    
    PassingVariantCallFilter(FacetFactory facet_factory,
                             std::vector<MeasureWrapper> measures,
                             OutputOptions output_config,
                             ConcurrencyPolicy threading,
                             boost::optional<ProgressMeter&> progress = boost::none);
    
    PassingVariantCallFilter(const PassingVariantCallFilter&)            = delete;
    PassingVariantCallFilter& operator=(const PassingVariantCallFilter&) = delete;
    PassingVariantCallFilter(PassingVariantCallFilter&&)                 = delete;
    PassingVariantCallFilter& operator=(PassingVariantCallFilter&&)      = delete;
    
    virtual ~PassingVariantCallFilter() override = default;

private:
    std::string do_name() const override;
    void annotate(VcfHeader::Builder& header) const override;
    Classification classify(const MeasureVector& measures) const override;
};

} // namespace csr
} // namespace octopus

#endif
