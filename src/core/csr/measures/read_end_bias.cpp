// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_end_bias.hpp"

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

const std::string ReadEndBias::name_ = "REB";

std::unique_ptr<Measure> ReadEndBias::do_clone() const
{
    return std::make_unique<ReadEndBias>(*this);
}

Measure::ResultType ReadEndBias::get_default_result() const
{
    return std::vector<double> {};
}

void ReadEndBias::do_set_parameters(std::vector<std::string> params)
{
    if (params.size() != 1) {
        throw BadMeasureParameters {this->name(), "only has one parameter (end fraction)"};
    }
    try {
        end_fraction_ = boost::lexical_cast<decltype(end_fraction_)>(params.front());
    } catch (const boost::bad_lexical_cast&) {
        throw BadMeasureParameters {this->name(), "given parameter \"" + params.front() + "\" cannot be parsed"};
    }
    if (end_fraction_ < 0 || end_fraction_ > 1) {
        throw BadMeasureParameters {this->name(), "end fraction must be between 0 and 1"};
    }
}

namespace {

struct EndDefinition
{
    double end_fraction = 0.03;
};

struct PositionCounts
{
    unsigned end, middle;
};

bool overlaps_rhs(const Allele& allele, const AlignedRead& read)
{
    using D = GenomicRegion::Distance;
    return overlaps(allele, expand_lhs(mapped_region(read), -static_cast<D>(region_size(read) / 2)));
}

auto head_region(const AlignedRead& read, const unsigned tail_bases)
{
    if (is_forward_strand(read)) {
        return octopus::head_region(read, tail_bases);
    } else {
        return octopus::tail_region(read, tail_bases);
    }
}


auto tail_region(const AlignedRead& read, const unsigned tail_bases)
{
    if (is_forward_strand(read)) {
        return octopus::tail_region(read, tail_bases);
    } else {
        return  octopus::head_region(read, tail_bases);
    }
}

bool is_in_end(const Allele& allele, const AlignedRead& read, const EndDefinition end_def = EndDefinition {})
{
    assert(contains(read, allele));
    const auto end_bases = std::max(static_cast<unsigned>(sequence_size(read) * end_def.end_fraction), 1u);
    return contains(head_region(read, end_bases), allele) || contains(tail_region(read, end_bases), allele);
}

auto compute_head_tail_counts(const Allele& allele, const ReadRefSupportSet& support, const EndDefinition end_def = EndDefinition {})
{
    PositionCounts result {};
    for (const auto& read : support) {
        if (contains(read, allele)) {
            if (is_in_end(allele, read)) {
                ++result.end;
            } else {
                ++result.middle;
            }
        }
    }
    return result;
}

double calculate_tail_bias(const PositionCounts counts, const double tolerance = 0.5, const double prior = 1.0)
{
    assert(tolerance > 0.0 && tolerance < 1.0);
    return maths::beta_cdf(counts.middle + prior, counts.end + prior, tolerance);
}

double calculate_tail_bias(const Allele& allele, const ReadRefSupportSet& support, const EndDefinition end_def = EndDefinition {})
{
    const auto counts = compute_head_tail_counts(allele, support);
    return calculate_tail_bias(counts, std::min(3 * end_def.end_fraction, 0.5));
}

double calculate_max_tail_bias(const AlleleSupportMap& support, const EndDefinition end_def = EndDefinition {})
{
    double result {0};
    for (const auto& p : support) {
        auto bias = calculate_tail_bias(p.first, p.second, end_def);
        result = std::max(result, bias);
    }
    return result;
}

} // namespace

Measure::ResultType ReadEndBias::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    std::vector<double> result {};
    result.reserve(samples.size());
    const EndDefinition end_def {end_fraction_};
    for (const auto& sample : samples) {
        std::vector<Allele> alleles; bool has_ref;
        std::tie(alleles, has_ref) = get_called_alleles(call, sample);
        if (!alleles.empty()) {
            const auto allele_support = compute_allele_support(alleles, assignments, sample);
            result.push_back(calculate_max_tail_bias(allele_support, end_def));
        } else {
            result.push_back(0);
        }
    }
    return result;
}

Measure::ResultCardinality ReadEndBias::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& ReadEndBias::do_name() const
{
    return name_;
}

std::string ReadEndBias::do_describe() const
{
    return "Probability allele occurs at the end (head or tail) of supporting reads";
}

std::vector<std::string> ReadEndBias::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}

bool ReadEndBias::is_equal(const Measure& other) const noexcept
{
    return end_fraction_ == static_cast<const ReadEndBias&>(other).end_fraction_;
}

} // namespace csr
} // namespace octopus
