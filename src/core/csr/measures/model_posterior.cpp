// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "model_posterior.hpp"

#include <boost/lexical_cast.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "config/octopus_vcf.hpp"
#include "../facets/samples.hpp"

namespace octopus { namespace csr {

const std::string ModelPosterior::name_ = "MP";

std::unique_ptr<Measure> ModelPosterior::do_clone() const
{
    return std::make_unique<ModelPosterior>(*this);
}

Measure::ValueType ModelPosterior::get_value_type() const
{
    return double {};
}

Measure::ResultType ModelPosterior::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    namespace ovcf = octopus::vcf::spec;
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    Array<Optional<ValueType>> result(samples.size());
    if (call.has_format(ovcf::format::modelPosterior)) {
        for (std::size_t s {0}; s < samples.size(); ++s) {
            const auto mp = call.get_sample_value(samples[s], ovcf::format::modelPosterior);
            assert(mp.size() == 1);
            if (mp[0] != vcfspec::missingValue) {
                result[s] = boost::lexical_cast<double>(mp[0]);
            }
        }
    } else if (!is_info_missing(ovcf::info::modelPosterior, call)) {
        auto mp = std::stod(call.info_value(ovcf::info::modelPosterior).front());
        std::fill(std::begin(result), std::end(result), mp);
    }
    return result;
}

Measure::ResultCardinality ModelPosterior::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& ModelPosterior::do_name() const
{
    return name_;
}

std::string ModelPosterior::do_describe() const
{
    return "Model posterior for this haplotype block";
}

std::vector<std::string> ModelPosterior::do_requirements() const
{
    return {"Samples"};
}

} // namespace csr
} // namespace octopus
