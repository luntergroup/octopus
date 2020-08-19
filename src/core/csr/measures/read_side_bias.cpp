// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_side_bias.hpp"

#include <algorithm>
#include <iterator>
#include <random>
#include <functional>
#include <cmath>
#include <cassert>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "utils/maths.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string ReadSideBias::name_ = "RSB";

std::unique_ptr<Measure> ReadSideBias::do_clone() const
{
    return std::make_unique<ReadSideBias>(*this);
}

Measure::ValueType ReadSideBias::get_value_type() const
{
    return double {};
}

namespace {

struct PositionCounts
{
    unsigned head, tail;
};

bool overlaps_rhs(const Allele& allele, const AlignedRead& read)
{
    using D = GenomicRegion::Distance;
    return overlaps(allele, expand_lhs(mapped_region(read), -static_cast<D>(region_size(read) / 2)));
}

void update_counts(const Allele& allele, const AlignedRead& read, PositionCounts& counts)
{
    if (overlaps(allele, read)) {
        if (region_size(allele) >= region_size(read) / 2) {
            ++counts.head;
            ++counts.tail;
        } else if (overlaps_rhs(allele, read)) {
            if (is_forward_strand(read)) {
                ++counts.tail;
            } else {
                ++counts.head;
            }
        } else {
            if (is_forward_strand(read)) {
                ++counts.head;
            } else {
                ++counts.tail;
            }
        }
    }
}

auto compute_allele_positions(const Allele& allele, const ReadRefSupportSet& support)
{
    PositionCounts forward_counts {}, reverse_counts {};
    for (const auto& read : support) {
        if (is_forward_strand(read)) {
            update_counts(allele, read, forward_counts);
        } else {
            update_counts(allele, read, reverse_counts);
        }
    }
    return std::make_pair(forward_counts, reverse_counts);
}

double calculate_position_bias(const PositionCounts forward_counts, const PositionCounts reverse_counts,
                               const double tolerance = 0.5)
{
    assert(tolerance > 0.0 && tolerance < 1.0);
    const auto num_lhs = forward_counts.head + reverse_counts.tail;
    const auto num_rhs = forward_counts.tail + reverse_counts.head;
    const auto prob_lhs_biased = maths::beta_sf(static_cast<double>(num_lhs + 1), static_cast<double>(num_rhs + 1), 0.5 + tolerance / 2);
    const auto prob_rhs_biased = maths::beta_cdf(static_cast<double>(num_lhs + 1), static_cast<double>(num_rhs + 1), 0.5 - tolerance / 2);
    return prob_lhs_biased + prob_rhs_biased;
}

double calculate_position_bias(const Allele& allele, const ReadRefSupportSet& support)
{
    PositionCounts forward_counts, reverse_counts;
    std::tie(forward_counts, reverse_counts) = compute_allele_positions(allele, support);
    return calculate_position_bias(forward_counts, reverse_counts);
}

double calculate_position_bias(const std::vector<Allele>& alleles, const AlleleSupportMap& support)
{
    double result {0};
    for (const auto& allele : alleles) {
        const auto support_set_itr = support.find(allele);
        if (support_set_itr != std::cend(support)) {
            auto bias = calculate_position_bias(allele, support_set_itr->second);
            result = std::max(result, bias);
        }
    }
    return result;
}

} // namespace

Measure::ResultType ReadSideBias::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    Array<ValueType> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.emplace_back(calculate_position_bias(get_called(alleles, call, sample), assignments.at(sample)));
    }
    return result;
}

Measure::ResultCardinality ReadSideBias::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& ReadSideBias::do_name() const
{
    return name_;
}

std::string ReadSideBias::do_describe() const
{
    return "Bias of variant side (head or tail half) in supporting reads";
}

std::vector<std::string> ReadSideBias::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
