// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "is_somatic.hpp"

#include "io/variant/vcf_record.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> IsSomatic::do_clone() const
{
    return std::make_unique<IsSomatic>(*this);
}

Measure::ResultType IsSomatic::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    return is_somatic(call);
}

Measure::ResultCardinality IsSomatic::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

std::string IsSomatic::do_name() const
{
    return "SOMATIC";
}

std::string IsSomatic::do_describe() const
{
    return "Is the call marked SOMATIC";
}

} // namespace csr
} // namespace octopus
