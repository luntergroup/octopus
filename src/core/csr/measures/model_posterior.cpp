// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "model_posterior.hpp"

#include "io/variant/vcf_record.hpp"
#include "config/octopus_vcf.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> ModelPosterior::do_clone() const
{
    return std::make_unique<ModelPosterior>(*this);
}

Measure::ResultType ModelPosterior::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    namespace ovcf = octopus::vcf::spec;
    if (call.has_info(ovcf::info::modelPosterior)) {
        return std::stod(call.info_value(ovcf::info::modelPosterior).front());
    } else {
        return boost::none;
    }
}

std::string ModelPosterior::do_name() const
{
    return "MP";
}

} // namespace csr
} // namespace octopus
