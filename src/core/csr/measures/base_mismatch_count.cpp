// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "base_mismatch_count.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cassert>

#include "io/variant/vcf_record.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "core/tools/read_assigner.hpp"
#include "utils/genotype_reader.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string BaseMismatchCount::name_ = "BMC";

std::unique_ptr<Measure> BaseMismatchCount::do_clone() const
{
    return std::make_unique<BaseMismatchCount>(*this);
}

namespace {

bool completely_contains(const AlignedRead& read, const Allele& allele)
{
    return contains(read, allele) && !are_adjacent(read, allele);
}

template <typename Iterator1, typename Iterator2>
auto count_mismatches(Iterator1 first1, Iterator1 last1, Iterator2 first2) noexcept
{
    return std::inner_product(first1, last1, first2, 0u, std::plus<> {}, std::not_equal_to<> {});
}

template <typename Range1, typename Range2>
auto count_mismatches_helper(const Range1& lhs, const Range2& rhs) noexcept
{
    return count_mismatches(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}

unsigned count_mismatches(const AlignedRead& read, const Allele& allele)
{
    if (!overlaps(read, allele)) return 0;
    const auto read_section = copy_sequence(read, mapped_region(allele));
    if (completely_contains(read, allele)) {
        if (read_section.size() == sequence_size(allele)) {
            return count_mismatches_helper(read_section, allele.sequence());
        } else if (read_section.size() < sequence_size(allele)) {
            return count_mismatches_helper(read_section, allele.sequence());// + sequence_size(allele) - read_section.size();
        } else {
            return count_mismatches_helper(allele.sequence(), read_section);// + read_section.size() - sequence_size(allele);
        }
    } else if (begins_before(read, allele)) {
        if (read_section.size() <= sequence_size(allele)) {
            return count_mismatches_helper(read_section, allele.sequence());
        } else {
            return count_mismatches_helper(allele.sequence(), read_section);
        }
    } else {
        if (read_section.size() <= sequence_size(allele)) {
            return count_mismatches(std::crbegin(read_section), std::crend(read_section), std::crbegin(allele.sequence()));
        } else {
            return count_mismatches(std::crbegin(allele.sequence()), std::crend(allele.sequence()), std::crbegin(read_section));
        }
    }
}

} // namespace

Measure::ResultType BaseMismatchCount::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
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
                    sample_result += count_mismatches(read, p.first);
                }
            }
        }
        result.push_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality BaseMismatchCount::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& BaseMismatchCount::do_name() const
{
    return name_;
}

std::string BaseMismatchCount::do_describe() const
{
    return "Number of base mismatches at variant position in reads supporting variant haplotype";
}

std::vector<std::string> BaseMismatchCount::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
