// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "duplicate_concordance.hpp"

#include <algorithm>
#include <iterator>
#include <cassert>

#include <boost/variant.hpp>

#include "basics/mappable_reference_wrapper.hpp"
#include "core/tools/read_assigner.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_duplicates.hpp"
#include "utils/maths.hpp"
#include "../facets/samples.hpp"
#include "../facets/overlapping_reads.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string DuplicateConcordance::name_ = "DPC";

std::unique_ptr<Measure> DuplicateConcordance::do_clone() const
{
    return std::make_unique<DuplicateConcordance>(*this);
}

Measure::ValueType DuplicateConcordance::get_value_type() const
{
    return double {};
}

namespace {

bool is_canonical(const VcfRecord::NucleotideSequence& allele) noexcept
{
    const static VcfRecord::NucleotideSequence deleted_allele {vcfspec::deletedBase};
    return !(allele == vcfspec::missingValue || allele == deleted_allele);
}

bool has_called_alt_allele(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    if (!call.has_genotypes()) return true;
    const auto& genotype = get_genotype(call, sample);
    return std::any_of(std::cbegin(genotype), std::cend(genotype),
                       [&] (const auto& allele) { return allele != call.ref() && is_canonical(allele); });
}

auto find_duplicate_overlapped_reads(const ReadContainer& reads, const GenomicRegion& region)
{
    const auto overlapped_reads = overlap_range(reads, region);
    const auto duplicate_itrs = find_duplicate_reads(std::cbegin(overlapped_reads), std::cend(overlapped_reads));
    std::vector<std::vector<AlignedRead>> result {};
    result.reserve(duplicate_itrs.size());
    for (const auto& itrs : duplicate_itrs) {
        std::vector<AlignedRead> dups {};
        dups.reserve(itrs.size());
        std::transform(std::cbegin(itrs), std::cend(itrs), std::back_inserter(dups), [] (auto itr) { return *itr; });
        result.push_back(std::move(dups));
    }
    return result;
}

bool other_segments_equal(const AlignedRead& lhs, const AlignedRead& rhs) noexcept
{
    if (lhs.has_other_segment()) {
        return rhs.has_other_segment() && lhs.next_segment() == rhs.next_segment();
    } else {
        return !rhs.has_other_segment();
    }
}

bool are_realigned_equal(const AlignedRead& lhs, const AlignedRead& rhs) noexcept
{
    return lhs.mapping_quality()    == rhs.mapping_quality()
           && lhs.name()            == rhs.name()
           && lhs.sequence()        == rhs.sequence()
           && lhs.base_qualities()  == rhs.base_qualities()
           && lhs.read_group()      == rhs.read_group()
           && lhs.flags()           == rhs.flags()
           && other_segments_equal(lhs, rhs);
}

bool is_duplicate(const AlignedRead& realigned_read, const std::vector<AlignedRead>& duplicate_reads)
{
    // Duplicate realigned reads may have different mapping position and cigar to the raw duplicate reads
    const auto is_duplicate = [&] (const auto& read) { return are_realigned_equal(read, realigned_read); };
    return std::find_if(std::cbegin(duplicate_reads), std::cend(duplicate_reads), is_duplicate) != std::cend(duplicate_reads);
}

double calculate_support_concordance(const std::vector<AlignedRead>& duplicate_reads, const AlleleSupportMap& allele_support)
{
    assert(duplicate_reads.size() > 1);
    std::vector<unsigned> support_counts {};
    unsigned total_duplicate_support {0};
    support_counts.reserve(allele_support.size());
    for (const auto& p : allele_support) {
        const auto is_duplicate_helper = [&] (const auto& read) { return is_duplicate(read, duplicate_reads); };
        const auto duplicate_support = std::count_if(std::cbegin(p.second), std::cend(p.second), is_duplicate_helper);
        if (duplicate_support > 0) {
            support_counts.push_back(duplicate_support);
            total_duplicate_support += duplicate_support;
        }
    }
    if (support_counts.size() < 2) {
        return 1;
    } else {
        std::vector<double> support_proportions(support_counts.size());
        std::transform(std::cbegin(support_counts), std::cend(support_counts), std::begin(support_proportions),
                       [=] (auto count) { return static_cast<double>(count) / total_duplicate_support; });
        return maths::entropy2(support_proportions);
    }
}

} // namespace

Measure::ResultType DuplicateConcordance::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    Array<Optional<ValueType>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        Optional<ValueType> sample_result {};
        if (has_called_alt_allele(call, sample)) {
            const auto duplicate_reads = find_duplicate_overlapped_reads(reads.at(sample), mapped_region(call));
            if (!duplicate_reads.empty()) {
                const auto allele_support = compute_allele_support(get_called_alleles(call, sample).first, assignments, sample);
                for (const auto& duplicates : duplicate_reads) {
                    const auto concordance = calculate_support_concordance(duplicates, allele_support);
                    if (!sample_result || concordance < boost::get<double>(*sample_result)) {
                        sample_result = concordance;
                    }
                }
            }
        }
        result.push_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality DuplicateConcordance::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& DuplicateConcordance::do_name() const
{
    return name_;
}

std::string DuplicateConcordance::do_describe() const
{
    return "Concordance of allele support from duplicate reads";
}

std::vector<std::string> DuplicateConcordance::do_requirements() const
{
    
    return {"Samples", "OverlappingReads", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
