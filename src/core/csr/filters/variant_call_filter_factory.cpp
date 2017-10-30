// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_call_filter_factory.hpp"

#include <utility>

#include "../facets/facet_factory.hpp"

namespace octopus { namespace csr {

std::unique_ptr<VariantCallFilterFactory> VariantCallFilterFactory::clone() const
{
    return do_clone();
}

std::unique_ptr<VariantCallFilter> VariantCallFilterFactory::make(const ReferenceGenome& reference,
                                                                  BufferedReadPipe read_pipe,
                                                                  VariantCallFilter::OutputOptions output_config,
                                                                  boost::optional<ProgressMeter&> progress) const
{
    FacetFactory facet_factory {reference, std::move(read_pipe)};
    return do_make(std::move(facet_factory), output_config, progress);
}

} // namespace csr
} // namespace octopus
