// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "duplicate_concordance.hpp"

#include <algorithm>
#include <iterator>
#include <unordered_set>
#include <cassert>

#include <boost/variant.hpp>

#include "basics/mappable_reference_wrapper.hpp"
#include "core/tools/read_assigner.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"
#include "utils/erase_if.hpp"
#include "../facets/samples.hpp"
#include "../facets/reads_summary.hpp"
#include "../facets/alleles.hpp"
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
    return !(allele == vcfspec::missingValue || allele == vcfspec::deleteMaskAllele);
}

bool has_called_alt_allele(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    if (!call.has_genotypes()) return true;
    const auto& genotype = get_genotype(call, sample);
    return std::any_of(std::cbegin(genotype), std::cend(genotype),
                       [&] (const auto& allele) { return allele != call.ref() && is_canonical(allele); });
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

using AlleleSupportNameSets = std::vector<std::unordered_set<std::string>>;

double calculate_support_concordance(const std::vector<AlignedRead>& duplicate_reads, AlleleSupportNameSets& support_names)
{
    assert(duplicate_reads.size() > 1);
    std::vector<unsigned> duplicate_allele_support_counts {};
    duplicate_allele_support_counts.reserve(support_names.size());
    unsigned total_duplicate_support {0};
    for (auto& allele_support_names : support_names) {
        unsigned allele_duplicate_support {0};
        for (const auto& read : duplicate_reads) {
            const auto support_itr = allele_support_names.find(read.name());
            if (support_itr != std::cend(allele_support_names)) {
                ++allele_duplicate_support;
                allele_support_names.erase(support_itr);
            }
        }
        if (allele_duplicate_support > 0) {
            duplicate_allele_support_counts.push_back(allele_duplicate_support);
            total_duplicate_support += allele_duplicate_support;
        }
    }
    if (total_duplicate_support > 0) {
        erase_if(support_names, [] (const auto& names) { return names.empty(); });
    }
    if (duplicate_allele_support_counts.size() < 2) {
        return 1;
    } else {
        std::vector<double> support_proportions(duplicate_allele_support_counts.size());
        std::transform(std::cbegin(duplicate_allele_support_counts), std::cend(duplicate_allele_support_counts),
                       std::begin(support_proportions),
                       [=] (auto count) { return static_cast<double>(count) / total_duplicate_support; });
        return 1 - std::min(maths::entropy2(support_proportions), 1.0);
    }
}

} // namespace

Measure::ResultType DuplicateConcordance::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& reads = get_value<ReadsSummary>(facets.at("ReadsSummary"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    Array<Optional<ValueType>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        Optional<ValueType> sample_result {};
        if (has_called_alt_allele(call, sample)) {
            const auto sample_alleles = get_called(alleles, call, sample);
            const auto& allele_support = assignments.at(sample);
            AlleleSupportNameSets supporting_read_name {};
            supporting_read_name.reserve(allele_support.size());
            for (const auto& p : allele_support) {
                std::unordered_set<std::string> names {};
                names.reserve(p.second.size());
                for (const auto& read : p.second) {
                    names.insert(read.get().name());
                }
            }
            for (const auto& duplicates : overlap_range(reads.at(sample).duplicates, call)) {
                const auto concordance = calculate_support_concordance(duplicates.reads, supporting_read_name);
                if (!sample_result || concordance < boost::get<double>(*sample_result)) {
                    sample_result = concordance;
                }
                if (supporting_read_name.size() < 2) break;
                if (concordance < 1e-3) break;
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
    
    return {"Samples", "ReadsSummary", "Alleles", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
