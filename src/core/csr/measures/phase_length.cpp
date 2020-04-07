// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "phase_length.hpp"

#include "io/variant/vcf_record.hpp"
#include "config/octopus_vcf.hpp"
#include "../facets/samples.hpp"
#include "../facets/genotypes.hpp"

namespace octopus { namespace csr {

const std::string PhaseLength::name_ = "PLN";

std::unique_ptr<Measure> PhaseLength::do_clone() const
{
    return std::make_unique<PhaseLength>(*this);
}

Measure::ResultType PhaseLength::get_default_result() const
{
    return std::vector<std::size_t> {};
}

Measure::ResultType PhaseLength::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& genotypes = get_value<Genotypes>(facets.at("Genotypes"));
    std::vector<std::size_t> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        const auto& sample_genotypes = genotypes.at(sample);
        result.push_back(!sample_genotypes.empty() ? size(encompassing_region(sample_genotypes)) : 0);
    }
    return result;
}

Measure::ResultCardinality PhaseLength::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& PhaseLength::do_name() const
{
    return name_;
}

std::string PhaseLength::do_describe() const
{
    return "Length of the phase block for the call";
}

std::vector<std::string> PhaseLength::do_requirements() const
{
    return {"Samples", "Genotypes"};
}

} // namespace csr
} // namespace octopus
