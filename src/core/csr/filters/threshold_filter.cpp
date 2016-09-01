// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "threshold_filter.hpp"

#include "io/variant/vcf_header.hpp"

namespace octopus { namespace csr {

ThresholdVariantCallFilter::ThresholdVariantCallFilter(const ReferenceGenome& reference,
                                                       const ReadPipe& read_pipe,
                                                       std::vector<MeasureWrapper> measures,
                                                       std::size_t max_read_buffer_size)
: VariantCallFilter {reference, read_pipe, std::move(measures), max_read_buffer_size}
{}

void ThresholdVariantCallFilter::annotate(VcfHeader& header) const
{
    // TODO
}

VariantCallFilter::Classification ThresholdVariantCallFilter::classify(const MeasureVector& call_measures) const
{
    return VariantCallFilter::Classification {};
}

} // namespace csr
} // namespace octopus
