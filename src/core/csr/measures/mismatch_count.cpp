// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mismatch_count.hpp"

#include <iterator>
#include <algorithm>
#include <cassert>

#include "io/variant/vcf_record.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "core/tools/read_assigner.hpp"
#include "utils/genotype_reader.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string MismatchCount::name_ = "MC";

std::unique_ptr<Measure> MismatchCount::do_clone() const
{
    return std::make_unique<MismatchCount>(*this);
}

Measure::ResultType MismatchCount::get_default_result() const
{
    return std::vector<int> {};
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
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    std::vector<int> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        std::vector<Allele> alleles; bool has_ref;
        std::tie(alleles, has_ref) = get_called_alleles(call, sample);
        int sample_result {0};
        if (!alleles.empty()) {
            const auto sample_allele_support = compute_allele_support(alleles, assignments, sample);
            for (const auto& p : sample_allele_support) {
                for (const auto& read : p.second) {
                    sample_result += mismatches(read, p.first);
                }
            }
        }
        result.push_back(sample_result);
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
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
