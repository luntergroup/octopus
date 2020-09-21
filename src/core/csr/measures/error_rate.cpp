// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_rate.hpp"

#include <numeric>
#include <algorithm>
#include <iterator>
#include <map>
#include <cstddef>

#include <boost/variant.hpp>

#include "basics/cigar_string.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/haplotype.hpp"
#include "io/variant/vcf_record.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string ErrorRate::name_ = "ER";

std::unique_ptr<Measure> ErrorRate::do_clone() const
{
    return std::make_unique<ErrorRate>(*this);
}

Measure::ValueType ErrorRate::get_value_type() const
{
    return double {};
}

namespace {

boost::optional<double> 
compute_error_rate(const Facet::SupportMaps::HaplotypeSupportMaps& assignments, const GenomicRegion& region) noexcept
{
    std::size_t error_bases {0}, total_bases {0};
    for (const auto& p : assignments.assigned_wrt_haplotype) {
        for (const auto& read : overlap_range(p.second, region)) {
            error_bases += sum_non_matches(read.cigar());
            total_bases += sequence_size(read);
        }
    }
    for (const auto& read : overlap_range(assignments.ambiguous_wrt_haplotype, region)) {
        error_bases += sum_non_matches(read.read.cigar());
        total_bases += sequence_size(read.read);
    }
    if (total_bases > 0) {
        return static_cast<double>(error_bases) / total_bases;
    } else {
        return boost::none;
    }
}

} // namespace

Measure::ResultType ErrorRate::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).haplotypes;
    Array<Optional<ValueType>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.emplace_back(compute_error_rate(assignments.at(sample), mapped_region(call)));
    }
    return result;
}

Measure::ResultCardinality ErrorRate::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& ErrorRate::do_name() const
{
    return name_;
}

std::string ErrorRate::do_describe() const
{
    return "Error rate in reads overlapping the site";
}

std::vector<std::string> ErrorRate::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}
    
} // namespace csr
} // namespace octopus
