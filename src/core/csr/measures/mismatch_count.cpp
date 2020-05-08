// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mismatch_count.hpp"

#include <iterator>
#include <algorithm>
#include <cassert>

#include "io/variant/vcf_record.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "core/tools/read_assigner.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string MismatchCount::name_ = "MC";

std::unique_ptr<Measure> MismatchCount::do_clone() const
{
    return std::make_unique<MismatchCount>(*this);
}

Measure::ValueType MismatchCount::get_value_type() const
{
    return int {};
}

namespace {

bool completely_contains(const AlignedRead& read, const Allele& allele)
{
    return contains(read, allele) && !are_adjacent(read, allele);
}

bool mismatches(const AlignedRead& read, const Allele& allele)
{
    if (!overlaps(read, allele)) return false;
    const auto read_section = copy_sequence(read, mapped_region(allele));
    if (completely_contains(read, allele)) {
        return read_section != allele.sequence();
    } else {
        if (read_section.size() > sequence_size(allele)) return true;
        if (begins_before(read, allele)) {
            return !std::equal(std::cbegin(read_section), std::cend(read_section), std::cbegin(allele.sequence()));
        } else {
            return !std::equal(std::crbegin(read_section), std::crend(read_section), std::crbegin(allele.sequence()));
        }
    }
}

} // namespace

Measure::ResultType MismatchCount::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    Array<ValueType> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        int sample_result {0};
        for (const auto& allele : get_called(alleles, call, sample)) {
            for (const auto& read : assignments.at(sample).at(allele)) {
                sample_result += mismatches(read, allele);
            }
        }
        result.emplace_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality MismatchCount::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& MismatchCount::do_name() const
{
    return name_;
}

std::string MismatchCount::do_describe() const
{
    return "Number of mismatches at variant position in reads supporting variant haplotype";
}

std::vector<std::string> MismatchCount::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
