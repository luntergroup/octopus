// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_tail_bias.hpp"

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

const std::string ReadTailBias::name_ = "RTB";

std::unique_ptr<Measure> ReadTailBias::do_clone() const
{
    return std::make_unique<ReadTailBias>(*this);
}

Measure::ValueType ReadTailBias::get_value_type() const
{
    return double {};
}

void ReadTailBias::do_set_parameters(std::vector<std::string> params)
{
    if (params.size() != 1) {
        throw BadMeasureParameters {this->name(), "only has one parameter (tail fraction)"};
    }
    try {
        tail_fraction_ = boost::lexical_cast<decltype(tail_fraction_)>(params.front());
    } catch (const boost::bad_lexical_cast&) {
        throw BadMeasureParameters {this->name(), "given parameter \"" + params.front() + "\" cannot be parsed"};
    }
    if (tail_fraction_ < 0 || tail_fraction_ > 1) {
        throw BadMeasureParameters {this->name(), "tail fraction must be between 0 and 1"};
    }
}

namespace {

struct TailDefinition
{
    double tail_fraction = 0.03;
};

struct PositionCounts
{
    unsigned head, tail;
};

bool overlaps_rhs(const Allele& allele, const AlignedRead& read)
{
    using D = GenomicRegion::Distance;
    return overlaps(allele, expand_lhs(mapped_region(read), -static_cast<D>(region_size(read) / 2)));
}

auto tail_region(const AlignedRead& read, const unsigned tail_bases)
{
    if (is_forward_strand(read)) {
        return octopus::tail_region(read, tail_bases);
    } else {
        return head_region(read, tail_bases);
    }
}

bool is_in_tail(const Allele& allele, const AlignedRead& read, const TailDefinition tail_def = TailDefinition {})
{
    assert(contains(read, allele));
    const auto tail_bases = std::max(static_cast<unsigned>(sequence_size(read) * tail_def.tail_fraction), 1u);
    return contains(tail_region(read, tail_bases), allele);
}

auto compute_head_tail_counts(const Allele& allele, const ReadRefSupportSet& support, const TailDefinition tail_def = TailDefinition {})
{
    PositionCounts result {};
    for (const auto& read : support) {
        if (contains(read, allele)) {
            if (is_in_tail(allele, read)) {
                ++result.tail;
            } else {
                ++result.head;
            }
        }
    }
    return result;
}

double calculate_tail_bias(const PositionCounts counts, const double tolerance = 0.5, const double prior = 1.0)
{
    assert(tolerance > 0.0 && tolerance < 1.0);
    return maths::beta_cdf(counts.head + prior, counts.tail + prior, tolerance);
}

double calculate_tail_bias(const Allele& allele, const ReadRefSupportSet& support, const TailDefinition tail_def = TailDefinition {})
{
    const auto counts = compute_head_tail_counts(allele, support);
    return calculate_tail_bias(counts, std::min(3 * tail_def.tail_fraction, 0.5));
}

double calculate_max_tail_bias(const std::vector<Allele>& alleles, const AlleleSupportMap& support, const TailDefinition tail_def = TailDefinition {})
{
    double result {0};
    for (const auto& allele : alleles) {
        const auto support_set_itr = support.find(allele);
        if (support_set_itr != std::cend(support)) {
            auto bias = calculate_tail_bias(allele, support_set_itr->second, tail_def);
            result = std::max(result, bias);
        }
    }
    return result;
}

} // namespace

Measure::ResultType ReadTailBias::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    Array<ValueType> result {};
    result.reserve(samples.size());
    const TailDefinition tail_def {tail_fraction_};
    for (const auto& sample : samples) {
        result.emplace_back(calculate_max_tail_bias(get_called(alleles, call, sample), assignments.at(sample), tail_def));
    }
    return result;
}

Measure::ResultCardinality ReadTailBias::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& ReadTailBias::do_name() const
{
    return name_;
}

std::string ReadTailBias::do_describe() const
{
    return "Probability allele occurs in the tail of supporting reads";
}

std::vector<std::string> ReadTailBias::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}

bool ReadTailBias::is_equal(const Measure& other) const noexcept
{
    return tail_fraction_ == static_cast<const ReadTailBias&>(other).tail_fraction_;
}

} // namespace csr
} // namespace octopus
