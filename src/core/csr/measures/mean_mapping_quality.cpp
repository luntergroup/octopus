// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mean_mapping_quality.hpp"

#include "io/variant/vcf_record.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> MeanMappingQuality::do_clone() const
{
    return std::make_unique<MeanMappingQuality>(*this);
}

Measure::ResultType MeanMappingQuality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    return std::stod(call.info_value("MQ").front());
}

std::string MeanMappingQuality::do_name() const
{
    return "MQ";
}

} // namespace csr
} // namespace octopus
