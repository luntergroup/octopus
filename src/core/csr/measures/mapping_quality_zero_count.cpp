// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mapping_quality_zero_count.hpp"

#include "io/variant/vcf_record.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> MappingQualityZeroCount::do_clone() const
{
    return std::make_unique<MappingQualityZeroCount>(*this);
}

Measure::ResultType MappingQualityZeroCount::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    return static_cast<std::size_t>(std::stoull(call.info_value("MQ0").front()));
}

std::string MappingQualityZeroCount::do_name() const
{
    return "MQ0";
}

} // namespace csr
} // namespace octopus
