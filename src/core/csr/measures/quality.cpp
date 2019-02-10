// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "quality.hpp"

#include "io/variant/vcf_record.hpp"

namespace octopus { namespace csr {

const std::string Quality::name_ = "QUAL";

std::unique_ptr<Measure> Quality::do_clone() const
{
    return std::make_unique<Quality>(*this);
}

Measure::ResultType Quality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    boost::optional<double> result {};
    if (call.qual()) {
        result = static_cast<double>(*call.qual());
    }
    return result;
}

Measure::ResultCardinality Quality::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& Quality::do_name() const
{
    return name_;
}

std::string Quality::do_describe() const
{
    return "Call QUAL";
}
    
} // namespace csr
} // namespace octopus
