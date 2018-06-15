// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_position_bias.hpp"

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
#include "utils/genotype_reader.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string ReadPositionBias::name_ = "RPB";

std::unique_ptr<Measure> ReadPositionBias::do_clone() const
{
    return std::make_unique<ReadPositionBias>(*this);
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

double calculate_position_bias(const AlleleSupportMap& support)
{
    double result {0};
    for (const auto& p : support) {
        auto bias = calculate_position_bias(p.first, p.second);
        result = std::max(result, bias);
    }
    return result;
}

} // namespace

Measure::ResultType ReadPositionBias::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).support;
    std::vector<double> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        std::vector<Allele> alleles; bool has_ref;
        std::tie(alleles, has_ref) = get_called_alleles(call, sample, true);
        assert(!alleles.empty());
        const auto allele_support = compute_allele_support(alleles, assignments.at(sample));
        auto position_bias = calculate_position_bias(allele_support);
        result.push_back(position_bias);
    }
    return result;
}

Measure::ResultCardinality ReadPositionBias::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& ReadPositionBias::do_name() const
{
    return name_;
}

std::string ReadPositionBias::do_describe() const
{
    return "Bias of variant position in supporting reads";
}

std::vector<std::string> ReadPositionBias::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
