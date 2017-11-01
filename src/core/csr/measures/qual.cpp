// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "qual.hpp"

#include "io/variant/vcf_record.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> Qual::do_clone() const
{
    return std::make_unique<Qual>(*this);
}

Measure::ResultType Qual::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    boost::optional<double> result {};
    if (call.qual()) {
        result = static_cast<double>(*call.qual());
    }
    return result;
}

std::string Qual::do_name() const
{
    return "QUAL";
}
    
} // namespace csr
} // namespace octopus
