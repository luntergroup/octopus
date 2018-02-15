// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "is_denovo.hpp"

#include "io/variant/vcf_record.hpp"
#include "config/octopus_vcf.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> IsDenovo::do_clone() const
{
    return std::make_unique<IsDenovo>(*this);
}

Measure::ResultType IsDenovo::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    return call.has_info(vcf::spec::info::denovo);
}

Measure::ResultCardinality IsDenovo::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

std::string IsDenovo::do_name() const
{
    return "DENOVO";
}

} // namespace csr
} // namespace octopus
