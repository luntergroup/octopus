// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "duplicate_allele_depth.hpp"

#include <algorithm>
#include <iterator>
#include <cassert>

#include <boost/variant.hpp>

#include "basics/mappable_reference_wrapper.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"
#include "../facets/samples.hpp"
#include "../facets/reads_summary.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string DuplicateAlleleDepth::name_ = "DAD";

std::unique_ptr<Measure> DuplicateAlleleDepth::do_clone() const
{
    return std::make_unique<DuplicateAlleleDepth>(*this);
}

Measure::ValueType DuplicateAlleleDepth::get_value_type() const
{
    return std::size_t {};
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

bool is_evaluable(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    return has_called_alt_allele(call, sample);
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

} // namespace

Measure::ResultType DuplicateAlleleDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& reads = get_value<ReadsSummary>(facets.at("ReadsSummary"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    Array<Optional<Array<ValueType>>> result(samples.size());
    for (std::size_t s {0}; s < samples.size(); ++s) {
        const auto& sample = samples[s];
        if (is_evaluable(call, sample)) {
            const auto sample_alleles = get_all(alleles, call, sample);
            result[s] = Array<ValueType>(sample_alleles.size(), std::size_t {0});
            const auto& duplicate_reads = overlap_range(reads.at(sample).duplicates, call);
            if (!duplicate_reads.empty()) {
                const auto& sample_support = assignments.at(sample);
                const auto compute_duplicate_support = [&] (const auto& allele) {
                    const auto& support = sample_support.at(allele);
                    std::size_t result {0};
                    for (const auto& duplicates : duplicate_reads) {
                        const auto is_duplicate_helper = [&] (const auto& read) { return is_duplicate(read, duplicates.reads); };
                        result += std::count_if(std::cbegin(support), std::cend(support), is_duplicate_helper);
                        --result; // One 'duplicate' read is not actually a duplicate
                    }
                    if (result < 0) result = 0;
                    return result;
                };
                std::transform(std::cbegin(sample_alleles), std::cend(sample_alleles), std::begin(*result[s]), compute_duplicate_support);
            }
        }
    }
    return result;
}

Measure::ResultCardinality DuplicateAlleleDepth::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alt_alleles;
}

const std::string& DuplicateAlleleDepth::do_name() const
{
    return name_;
}

std::string DuplicateAlleleDepth::do_describe() const
{
    return "Number of realigned reads supporting ALT alleles identified as duplicates";
}

std::vector<std::string> DuplicateAlleleDepth::do_requirements() const
{
    
    return {"Samples", "ReadsSummary", "Alleles", "ReadAssignments"};
}
    
} // namespace csr
} // namespace octopus
