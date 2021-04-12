// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "phylogeny_posterior.hpp"

#include "io/variant/vcf_record.hpp"

namespace octopus { namespace csr {

const std::string PhylogenyPosterior::name_ = "PPP";

std::unique_ptr<Measure> PhylogenyPosterior::do_clone() const
{
    return std::make_unique<PhylogenyPosterior>(*this);
}

Measure::ValueType PhylogenyPosterior::get_value_type() const
{
    return double {};
}

Measure::ResultType PhylogenyPosterior::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    Optional<ValueType> result {};
    if (!is_info_missing(this->name(), call)) {
        result = std::stod(call.info_value(this->name()).front());
    }
    return result;
}

Measure::ResultCardinality PhylogenyPosterior::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& PhylogenyPosterior::do_name() const
{
    return name_;
}

std::string PhylogenyPosterior::do_describe() const
{
    return "Phylogeny posterior probability";
}
    
} // namespace csr
} // namespace octopus
