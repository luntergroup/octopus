// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_call_filter_factory.hpp"

#include <utility>

#include "../facets/facet_factory.hpp"

namespace octopus { namespace csr {

VariantCallFilterFactory::VariantCallFilterFactory(VariantCallFilter::OutputOptions output_options)
: output_options_ {std::move(output_options)}
{}

std::unique_ptr<VariantCallFilterFactory> VariantCallFilterFactory::clone() const
{
    return do_clone();
}

void VariantCallFilterFactory::set_output_options(VariantCallFilter::OutputOptions output_options)
{
    output_options_ = std::move(output_options);
}

std::unique_ptr<VariantCallFilter>
VariantCallFilterFactory::make(const ReferenceGenome& reference,
                               BufferedReadPipe read_pipe,
                               VcfHeader input_header,
                               PloidyMap ploidies,
                               HaplotypeLikelihoodModel likelihood_model,
                               boost::optional<Pedigree> pedigree,
                               VariantCallFilter::OutputOptions output_config,
                               boost::optional<ProgressMeter&> progress,
                               boost::optional<unsigned> max_threads) const
{
    if (pedigree) {
        FacetFactory facet_factory {std::move(input_header), reference, std::move(read_pipe), std::move(ploidies), std::move(likelihood_model), std::move(*pedigree)};
        return do_make(std::move(facet_factory), output_config, progress, {max_threads});
    } else {
        FacetFactory facet_factory {std::move(input_header), reference, std::move(read_pipe), std::move(ploidies), std::move(likelihood_model)};
        return do_make(std::move(facet_factory), output_config, progress, {max_threads});
    }
}

std::unique_ptr<VariantCallFilter>
VariantCallFilterFactory::make(const ReferenceGenome& reference,
                               BufferedReadPipe read_pipe,
                               VcfHeader input_header,
                               PloidyMap ploidies,
                               HaplotypeLikelihoodModel likelihood_model,
                               boost::optional<Pedigree> pedigree,
                               boost::optional<ProgressMeter&> progress,
                               boost::optional<unsigned> max_threads) const
{
    return make(reference, std::move(read_pipe), std::move(input_header), std::move(ploidies),
                std::move(likelihood_model), std::move(pedigree), output_options_, progress, max_threads);
}

} // namespace csr
} // namespace octopus
