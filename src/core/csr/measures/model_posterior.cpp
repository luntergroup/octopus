// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "model_posterior.hpp"

#include "io/variant/vcf_record.hpp"
#include "config/octopus_vcf.hpp"

namespace octopus { namespace csr {

const std::string ModelPosterior::name_ = "MP";

std::unique_ptr<Measure> ModelPosterior::do_clone() const
{
    return std::make_unique<ModelPosterior>(*this);
}

Measure::ResultType ModelPosterior::get_default_result() const
{
    return boost::optional<double> {};
}

Measure::ResultType ModelPosterior::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    namespace ovcf = octopus::vcf::spec;
    boost::optional<double> result {};
    if (!is_info_missing(ovcf::info::modelPosterior, call)) {
        result = std::stod(call.info_value(ovcf::info::modelPosterior).front());
    }
    return result;
}

Measure::ResultCardinality ModelPosterior::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& ModelPosterior::do_name() const
{
    return name_;
}

std::string ModelPosterior::do_describe() const
{
    return "Model posterior for this haplotype block";
}

} // namespace csr
} // namespace octopus
