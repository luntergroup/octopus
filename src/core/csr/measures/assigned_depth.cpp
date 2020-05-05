// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "assigned_depth.hpp"

#include <numeric>
#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string AssignedDepth::name_ = "ADP";

std::unique_ptr<Measure> AssignedDepth::do_clone() const
{
    return std::make_unique<AssignedDepth>(*this);
}

Measure::ValueType AssignedDepth::get_value_type() const
{
    return std::size_t {};
}

namespace {

auto sum_support_counts(const std::vector<Allele>& alleles, const AlleleSupportMap& support)
{
    const auto add_support_count = [&] (auto total, const auto& allele) { return total + support.at(allele).size(); };
    return std::accumulate(std::cbegin(alleles), std::cend(alleles), std::size_t {0}, add_support_count);
}

} // namespace

Measure::ResultType AssignedDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    Array<ValueType> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        const auto sample_alleles = get_all(alleles, call, sample);
        result.emplace_back(sum_support_counts(sample_alleles, assignments.at(sample)));
    }
    return result;
}

Measure::ResultCardinality AssignedDepth::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& AssignedDepth::do_name() const
{
    return name_;
}

std::string AssignedDepth::do_describe() const
{
    return "Number of reads overlapping the position that could be assigned to an allele";
}

std::vector<std::string> AssignedDepth::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}
    
} // namespace csr
} // namespace octopus
