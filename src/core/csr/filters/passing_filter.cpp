// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "passing_filter.hpp"

namespace octopus { namespace csr {

PassingVariantCallFilter::PassingVariantCallFilter(FacetFactory facet_factory,
                                                   std::vector<MeasureWrapper> measures,
                                                   OutputOptions output_config,
                                                   ConcurrencyPolicy threading,
                                                   boost::optional<ProgressMeter&> progress)
: SinglePassVariantCallFilter {std::move(facet_factory), std::move(measures), output_config, threading, progress}
{}

void PassingVariantCallFilter::annotate(VcfHeader::Builder& header) const {}

VariantCallFilter::Classification PassingVariantCallFilter::classify(const MeasureVector& measures) const
{
    return {Classification::Category::unfiltered};
}

} // namespace csr
} // namespace octopus
