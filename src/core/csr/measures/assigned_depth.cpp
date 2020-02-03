// Copyright (c) 2015-2019 Daniel Cooke
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

Measure::ResultType AssignedDepth::get_default_result() const
{
    return std::vector<std::size_t> {};
}

namespace {

template <typename Map>
std::size_t sum_value_sizes(const Map& map) noexcept
{
    return std::accumulate(std::cbegin(map), std::cend(map), std::size_t {0},
                           [] (auto curr, const auto& p) noexcept { return curr + p.second.size(); });
}

} // namespace

Measure::ResultType AssignedDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    std::vector<std::size_t> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        const auto sample_alleles = copy_overlapped(alleles.at(sample), call);
        const auto allele_support = compute_allele_support(sample_alleles, assignments, sample);
        result.push_back(sum_value_sizes(allele_support));
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
