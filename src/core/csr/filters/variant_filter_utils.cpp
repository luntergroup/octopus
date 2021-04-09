// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_filter_utils.hpp"

#include <memory>

#include "../facets/facet_factory.hpp"
#include "../measures/is_somatic.hpp"
#include "../measures/is_denovo.hpp"
#include "threshold_filter.hpp"

namespace octopus { namespace csr {

namespace {

template <typename Flag>
void copy_flag(const VcfReader& source, VcfWriter& dest, Flag flag)
{
    FacetFactory facet_factory {source.fetch_header()};
    ThresholdVariantCallFilter::OutputOptions output_config {};
    VariantCallFilter::ConcurrencyPolicy thread_policy {};
    thread_policy.max_threads = 1;
    ThresholdVariantCallFilter::ConditionVectorPair conditions {};
    conditions.hard.push_back({std::move(flag), make_wrapped_threshold<EqualThreshold<bool>>(false)});
    std::unique_ptr<VariantCallFilter> filter = std::make_unique<ThresholdVariantCallFilter>(std::move(facet_factory), conditions, output_config, thread_policy);
    filter->filter(source, dest);
}

} // namesapce

void copy_somatics(const VcfReader& source, VcfWriter& dest)
{
    copy_flag(source, dest, make_wrapped_measure<IsSomatic>(false));
}

void copy_somatics(VcfReader::Path source, VcfReader::Path dest)
{
    VcfReader src {std::move(source)};
    VcfWriter dst {std::move(dest)};
    copy_somatics(src, dst);
}

void copy_denovos(const VcfReader& source, VcfWriter& dest)
{
    copy_flag(source, dest, make_wrapped_measure<IsDenovo>());
}

void copy_denovos(VcfReader::Path source, VcfReader::Path dest)
{
    VcfReader src {std::move(source)};
    VcfWriter dst {std::move(dest)};
    copy_denovos(src, dst);
}

} // namespace csr
} // namespace octopus
