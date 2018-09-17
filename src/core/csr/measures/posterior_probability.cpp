// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "posterior_probability.hpp"

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"

namespace octopus { namespace csr {

const std::string PosteriorProbability::name_ = "PP";

std::unique_ptr<Measure> PosteriorProbability::do_clone() const
{
    return std::make_unique<PosteriorProbability>(*this);
}

Measure::ResultType PosteriorProbability::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    boost::optional<double> result {};
    if (call.has_info("PP")) {
        const auto& pp = call.info_value("PP");
        if (pp.size() == 1 && pp.front() != vcfspec::missingValue) {
            result = std::stod(pp.front());
        }
    }
    return result;
}

Measure::ResultCardinality PosteriorProbability::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& PosteriorProbability::do_name() const
{
    return name_;
}

std::string PosteriorProbability::do_describe() const
{
    return "Call posterior probability";
}

} // namespace csr
} // namespace octopus
