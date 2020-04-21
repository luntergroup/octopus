// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "base_mismatch_quality.hpp"

#include <deque>
#include <iterator>
#include <algorithm>
#include <functional>
#include <cassert>

#include "io/variant/vcf_record.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string BaseMismatchQuality::name_ = "BMQ";

std::unique_ptr<Measure> BaseMismatchQuality::do_clone() const
{
    return std::make_unique<BaseMismatchQuality>(*this);
}

Measure::ResultType BaseMismatchQuality::get_default_result() const
{
    return std::vector<int> {};
}

namespace {

bool completely_contains(const AlignedRead& read, const Allele& allele)
{
    return contains(read, allele) && !are_adjacent(read, allele);
}

template <typename Iterator1, typename Iterator2, typename Iterator3, typename OutputIterator>
auto
copy_mismatches(Iterator1 first1, Iterator1 last1, Iterator2 first2,
                Iterator3 first_base_quality, OutputIterator result)
{
    for (; first1 != last1; ++first1, ++first2, ++first_base_quality) {
        if (*first1 != *first2) {
            *result++ = *first_base_quality;
        }
    }
}

template <typename Range1, typename Range2, typename Range3, typename OutputIterator>
auto copy_mismatches(const Range1& lhs, const Range2& rhs,
                     const Range3& base_qualities, OutputIterator result)
{
    return copy_mismatches(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs),
                           std::cbegin(base_qualities), result);
}

template <typename OutputIterator>
void copy_mismatch_base_quality(const AlignedRead& read, const Allele& allele, OutputIterator result)
{
    if (!overlaps(read, allele)) return;
    const auto read_section = copy_sequence(read, mapped_region(allele));
    const auto read_base_qualities = copy_base_qualities(read, mapped_region(allele));
    if (completely_contains(read, allele)) {
        if (read_section.size() == sequence_size(allele)) {
            copy_mismatches(read_section, allele.sequence(), read_base_qualities, result);
        } else if (read_section.size() < sequence_size(allele)) {
            copy_mismatches(read_section, allele.sequence(), read_base_qualities, result);
        } else {
            copy_mismatches(allele.sequence(), read_section, read_base_qualities, result);
        }
    } else if (begins_before(read, allele)) {
        if (read_section.size() <= sequence_size(allele)) {
            copy_mismatches(read_section, allele.sequence(), read_base_qualities, result);
        } else {
            copy_mismatches(allele.sequence(), read_section, read_base_qualities, result);
        }
    } else {
        if (read_section.size() <= sequence_size(allele)) {
            copy_mismatches(std::crbegin(read_section), std::crend(read_section), std::crbegin(allele.sequence()),
                            std::crbegin(read_base_qualities), result);
        } else {
            copy_mismatches(std::crbegin(allele.sequence()), std::crend(allele.sequence()), std::crbegin(read_section),
                            std::crbegin(read_base_qualities), result);
        }
    }
}

} // namespace

Measure::ResultType BaseMismatchQuality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    std::vector<boost::optional<int>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        const auto sample_alleles = get_all(alleles, call, sample);
        boost::optional<int> sample_result {};
        if (!sample_alleles.empty()) {
            std::deque<int> mismatch_qualities {};
            const auto& sample_support = assignments.at(sample);
            for (const auto& allele : sample_alleles) {
                for (const auto& read : sample_support.at(allele)) {
                    copy_mismatch_base_quality(read, allele, std::back_inserter(mismatch_qualities));
                }
            }
            if (!mismatch_qualities.empty()) {
                sample_result = maths::median(mismatch_qualities);
            }
        }
        result.push_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality BaseMismatchQuality::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& BaseMismatchQuality::do_name() const
{
    return name_;
}

std::string BaseMismatchQuality::do_describe() const
{
    return "Median quality of base mismatches in reads assigned to a unique allele";
}

std::vector<std::string> BaseMismatchQuality::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
