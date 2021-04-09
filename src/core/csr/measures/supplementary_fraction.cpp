// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "supplementary_fraction.hpp"

#include <iterator>
#include <algorithm>
#include <cassert>

#include "io/variant/vcf_record.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "utils/read_stats.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string SupplementaryFraction::name_ = "SF";

std::unique_ptr<Measure> SupplementaryFraction::do_clone() const
{
    return std::make_unique<SupplementaryFraction>(*this);
}

Measure::ValueType SupplementaryFraction::get_value_type() const
{
    return double {};
}

Measure::ResultType SupplementaryFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    Array<Optional<ValueType>> result(samples.size());
    for (std::size_t s {0}; s < samples.size(); ++s) {
        const auto& sample = samples[s];
        const auto sample_alleles = get_called_alt(alleles, call, sample);
        if (!sample_alleles.empty()) {
            double sample_result {0};
            const auto& support = assignments.at(sample);
            for (const auto& allele : sample_alleles) {
                const auto support_set_itr = support.find(allele);
                if (support_set_itr != std::cend(support) && !support_set_itr->second.empty()) {
                    auto allele_supplementary_fraction = static_cast<double>(count_supplementary(support_set_itr->second)) / support_set_itr->second.size();
                    sample_result += std::max(sample_result, allele_supplementary_fraction);
                }
            }
            result[s] = sample_result;
        }
    }
    return result;
}

Measure::ResultCardinality SupplementaryFraction::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& SupplementaryFraction::do_name() const
{
    return name_;
}

std::string SupplementaryFraction::do_describe() const
{
    return "Max fraction of reads supporting ALT alleles that are supplementary";
}

std::vector<std::string> SupplementaryFraction::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
