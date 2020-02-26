// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_rate_stdev.hpp"

#include <deque>
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
#include "utils/maths.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string ErrorRateStdev::name_ = "ERS";

std::unique_ptr<Measure> ErrorRateStdev::do_clone() const
{
    return std::make_unique<ErrorRateStdev>(*this);
}

Measure::ResultType ErrorRateStdev::get_default_result() const
{
    return std::vector<boost::optional<double>> {};
}

namespace {


auto calculate_error_rate(const AlignedRead& read) noexcept
{
    return static_cast<double>(sum_non_matches(read.cigar())) / sequence_size(read);
}

boost::optional<double>
compute_error_rate_stdev(const Facet::SupportMaps::HaplotypeSupportMaps& assignments, const GenomicRegion& region) noexcept
{
    std::deque<double> error_rates {};
    for (const auto& p : assignments.assigned_wrt_haplotype) {
        for (const auto& read : overlap_range(p.second, region)) {
            error_rates.push_back(calculate_error_rate(read));
        }
    }
    for (const auto& read : overlap_range(assignments.ambiguous_wrt_haplotype, region)) {
        error_rates.push_back(calculate_error_rate(read.read));
    }
    if (!error_rates.empty()) {
        return maths::stdev(error_rates);
    } else {
        return boost::none;
    }
}

} // namespace

Measure::ResultType ErrorRateStdev::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).haplotypes;
    std::vector<boost::optional<double>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.push_back(compute_error_rate_stdev(assignments.at(sample), mapped_region(call)));
    }
    return result;
}

Measure::ResultCardinality ErrorRateStdev::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& ErrorRateStdev::do_name() const
{
    return name_;
}

std::string ErrorRateStdev::do_describe() const
{
    return "Error rate standard deviation in reads overlapping the site";
}

std::vector<std::string> ErrorRateStdev::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
