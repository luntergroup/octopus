// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "supplementary_fraction.hpp"

#include <iterator>
#include <algorithm>
#include <cassert>

#include "io/variant/vcf_record.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/read_stats.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string SupplementaryFraction::name_ = "SF";

std::unique_ptr<Measure> SupplementaryFraction::do_clone() const
{
    return std::make_unique<SupplementaryFraction>(*this);
}

Measure::ResultType SupplementaryFraction::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    std::vector<boost::optional<double>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        std::vector<Allele> alleles; bool has_ref;
        std::tie(alleles, has_ref) = get_called_alleles(call, sample);
        if (has_ref) alleles.erase(std::cbegin(alleles));
        boost::optional<double> sample_result {};
        if (!alleles.empty()) {
            const auto sample_allele_support = compute_allele_support(alleles, assignments, sample);
            sample_result = 0;
            for (const auto& p : sample_allele_support) {
                if (!p.second.empty()) {
                    auto allele_supplementary_fraction = static_cast<double>(count_supplementary(p.second)) / p.second.size();
                    *sample_result += std::max(*sample_result, allele_supplementary_fraction);
                }
            }
        }
        result.push_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality SupplementaryFraction::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
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
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
