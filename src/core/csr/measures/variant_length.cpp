// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_length.hpp"

#include <algorithm>

#include "io/variant/vcf_record.hpp"
#include "utils/genotype_reader.hpp"
#include "../facets/samples.hpp"

namespace octopus { namespace csr {

const std::string VariantLength::name_ = "VL";

std::unique_ptr<Measure> VariantLength::do_clone() const
{
    return std::make_unique<VariantLength>(*this);
}

Measure::ResultType VariantLength::get_default_result() const
{
    return std::vector<int> {};
}

Measure::ResultType VariantLength::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    std::vector<int> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        std::vector<Allele> alleles; bool has_ref;
        std::tie(alleles, has_ref) = get_called_alleles(call, sample);
        int sample_result {0};
        for (const auto& allele : alleles) {
            sample_result = std::max({sample_result, static_cast<int>(region_size(allele)), static_cast<int>(sequence_size(allele))});
        }
        result.push_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality VariantLength::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& VariantLength::do_name() const
{
    return name_;
}

std::string VariantLength::do_describe() const
{
    return "Maximum length of called alleles";
}

std::vector<std::string> VariantLength::do_requirements() const
{
    return {"Samples"};
}

} // namespace csr
} // namespace octopus
