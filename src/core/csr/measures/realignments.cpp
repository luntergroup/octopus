// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "realignments.hpp"

#include <algorithm>
#include <iterator>
#include <cassert>
#include <array>
#include <iostream>
#include <sstream>

#include <boost/variant.hpp>
#include <boost/any.hpp>

#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "core/tools/read_assigner.hpp"
#include "core/tools/read_realigner.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/overlapping_reads.hpp"
#include "../facets/read_assignments.hpp"
#include "../facets/reference_context.hpp"
#include "utils/append.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> Realignments::do_clone() const
{
    return std::make_unique<Realignments>(*this);
}

namespace {

void remove_partial_alleles(std::vector<VcfRecord::NucleotideSequence>& genotype)
{
    genotype.erase(std::remove_if(std::begin(genotype), std::end(genotype),
                                  [] (const auto& seq) {
                                      static const std::string deleted_sequence {vcfspec::deletedBase};
                                      return seq == deleted_sequence || seq == vcfspec::missingValue;
                                  }), std::end(genotype));
}

auto num_matching_lhs_bases(const VcfRecord::NucleotideSequence& lhs, const VcfRecord::NucleotideSequence& rhs) noexcept
{
    auto p = std::mismatch(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
    return static_cast<int>(std::distance(std::cbegin(lhs), p.first));
}

auto get_called_alleles(const VcfRecord& call, const VcfRecord::SampleName& sample, const bool trim_padding = false)
{
    auto genotype = get_genotype(call, sample);
    remove_partial_alleles(genotype);
    std::sort(std::begin(genotype), std::end(genotype));
    genotype.erase(std::unique(std::begin(genotype), std::end(genotype)), std::end(genotype));
    const auto call_region = mapped_region(call);
    std::vector<Allele> result {};
    result.reserve(genotype.size());
    bool has_ref {false};
    if (trim_padding) {
        auto first_itr = std::begin(genotype);
        const auto ref_itr = std::find(first_itr, std::end(genotype), call.ref());
        if (ref_itr != std::end(genotype)) {
            if (ref_itr != first_itr) std::iter_swap(first_itr, ref_itr);
            ++first_itr;
            has_ref = true;
        }
        auto min_removed_bases = static_cast<int>(call.ref().size());
        std::for_each(first_itr, std::end(genotype), [&] (auto& allele) {
            const auto num_bases_to_remove = num_matching_lhs_bases(call.ref(), allele);
            allele.erase(std::cbegin(allele), std::next(std::cbegin(allele), num_bases_to_remove));
            auto allele_region = expand_lhs(call_region, -num_bases_to_remove);
            result.emplace_back(std::move(allele_region), std::move(allele));
            min_removed_bases = std::min(min_removed_bases, num_bases_to_remove);
        });
        if (has_ref) {
            auto& ref = genotype.front();
            ref.erase(std::cbegin(ref), std::next(std::cbegin(ref), min_removed_bases));
            auto allele_region = expand_lhs(call_region, -min_removed_bases);
            result.emplace_back(std::move(allele_region), std::move(ref));
            std::rotate(std::rbegin(result), std::next(std::rbegin(result)), std::rend(result));
        }
    } else {
        const auto ref_itr = std::find(std::begin(genotype), std::end(genotype), call.ref());
        if (ref_itr != std::end(genotype)) {
            if (ref_itr != std::begin(genotype)) std::iter_swap(std::begin(genotype), ref_itr);
            has_ref = true;
        }
        std::transform(std::cbegin(genotype), std::cend(genotype), std::back_inserter(result),
                       [&] (const auto& alt_seq) { return Allele {call_region, alt_seq}; });
    }
    return std::make_pair(std::move(result), has_ref);
}

using ReadRealignments = std::vector<std::vector<AlignedRead>>;

bool is_forward(const AlignedRead& read) noexcept
{
    return read.direction() == AlignedRead::Direction::forward;
}

struct DirectionCounts
{
    unsigned forward, reverse;
};

DirectionCounts count_directions(const std::vector<AlignedRead>& reads)
{
    auto n_fwd = static_cast<unsigned>(std::count_if(std::cbegin(reads), std::cend(reads),
                                                     [] (const auto& read) { return is_forward(read); }));
    return {n_fwd, static_cast<unsigned>(reads.size()) - n_fwd};
}

using MQHistogram = std::array<unsigned, 6>;

MQHistogram compute_mq_hist(const std::vector<AlignedRead>& reads)
{
    MQHistogram result {};
    for (const auto& read : reads) {
        ++result[std::min(read.mapping_quality(), AlignedRead::MappingQuality {59}) / 10];
    }
    return result;
}

struct Pileup
{
    std::array<unsigned, 4> match_counts, mismatch_counts;
    unsigned insertion_count, deletion_count;
    unsigned match_quality_sum, mismatch_quality_sum, insertion_quality_sum;
};

constexpr std::size_t window_size {61};

using PileupWindow = std::array<Pileup, window_size>;
using ExpandedCigarString = std::vector<CigarOperation::Flag>;

auto expand_cigar(const CigarString& cigar)
{
    ExpandedCigarString result {};
    result.reserve(sum_operation_sizes(cigar));
    for (const auto& op : cigar) {
        utils::append(result, op.size(), op.flag());
    }
    return result;
}

PileupWindow make_pileup(const std::vector<AlignedRead>& reads, const GenomicRegion& region)
{
    assert(size(region) == window_size);
    PileupWindow result {};
    for (const auto& read : reads) {
        if (overlaps(read, region)) {
            const auto read_fragment = copy(read, region);
            const auto expanded_cigar = expand_cigar(read_fragment.cigar());
            const auto fragment_offset = static_cast<unsigned>(begin_distance(region, read_fragment));
            unsigned read_pos {0}, pileup_pos {fragment_offset};
            for (const auto flag : expanded_cigar) {
                switch (flag) {
                    assert(pileup_pos < result.size());
                    using Flag = CigarOperation::Flag;
                    case Flag::sequenceMatch:
                    {
                        assert(read_pos < sequence_size(read_fragment));
                        switch (read_fragment.sequence()[read_pos]) {
                            case 'A': ++result[pileup_pos].match_counts[0]; break;
                            case 'C': ++result[pileup_pos].match_counts[1]; break;
                            case 'G': ++result[pileup_pos].match_counts[2]; break;
                            case 'T': ++result[pileup_pos].match_counts[3]; break;
                        }
                        result[pileup_pos].match_quality_sum += read_fragment.base_qualities()[read_pos];
                        ++pileup_pos; ++read_pos;
                        break;
                    }
                    case Flag::substitution:
                    {
                        assert(read_pos < sequence_size(read_fragment));
                        switch (read_fragment.sequence()[read_pos]) {
                            case 'A': ++result[pileup_pos].mismatch_counts[0]; break;
                            case 'C': ++result[pileup_pos].mismatch_counts[1]; break;
                            case 'G': ++result[pileup_pos].mismatch_counts[2]; break;
                            case 'T': ++result[pileup_pos].mismatch_counts[3]; break;
                        }
                        result[pileup_pos].mismatch_quality_sum += read_fragment.base_qualities()[read_pos];
                        ++pileup_pos; ++read_pos;
                        break;
                    }
                    case Flag::insertion:
                    {
                        ++result[pileup_pos].insertion_count;
                        result[pileup_pos].insertion_quality_sum += read_fragment.base_qualities()[read_pos];
                        ++read_pos;
                        break;
                    }
                    case Flag::deletion:
                    {
                        ++result[pileup_pos].deletion_count;
                        ++pileup_pos;
                        break;
                    }
                    default: continue;
                }
            }
        }
    }
    return result;
}

struct HaplotypeSummary
{
    DirectionCounts strand_counts;
    MQHistogram mq_hist;
    PileupWindow pileups;
};

using RealignmentSummary = std::vector<HaplotypeSummary>;

auto compute_realignment_summary(const ReadRealignments& realignments, const GenomicRegion& region)
{
    RealignmentSummary result {};
    result.reserve(realignments.size());
    for (const auto& reads : realignments) {
        result.push_back({count_directions(reads), compute_mq_hist(reads), make_pileup(reads, region)});
    }
    return result;
}

bool is_padded(const VcfRecord& call) noexcept
{
    const auto& ref = call.ref();
    if (ref.empty()) return false;
    return std::any_of(std::cbegin(call.alt()), std::cend(call.alt()),
                       [&] (const auto& alt) { return !alt.empty() && alt.front() == ref.front(); });
}

auto compute_realignment_summary(const ReadRealignments& realignments, const VcfRecord& call)
{
    auto call_position = head_region(call, 1);
    if (is_padded(call)) {
        call_position = shift(call_position, 1);
    }
    return compute_realignment_summary(realignments, expand(head_position(call), window_size / 2));
}

} // namespace

Measure::ResultType Realignments::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto assignments = boost::get<ReadAssignments::ResultType>(facets.at("ReadAssignments").get());
    assert(assignments.size() == 1);
    const auto& sample = assignments.cbegin()->first;
    const auto& support = assignments.cbegin()->second;
    std::vector<Allele> alleles; bool has_ref;
    std::tie(alleles, has_ref) = get_called_alleles(call, sample, true);
    ReadRealignments result {};
    result.reserve(alleles.size() + 1);
    std::vector<AlignedRead> assigned_reads {};
    for (const auto& allele : alleles) {
        std::vector<AlignedRead> allele_realignments {};
        for (const auto& h : support) {
            const auto& haplotype = h.first;
            const auto& supporting_reads = h.second;
            if (!supporting_reads.empty() && haplotype.contains(allele)) {
                const auto expanded_haplotype = expand_for_realignment(haplotype, supporting_reads);
                auto realignments = realign(supporting_reads, expanded_haplotype);
                for (std::size_t i {0}; i < realignments.size(); ++i) {
                    const auto& original_read = supporting_reads[i];
                    const auto& realigned_read = realignments[i];
                    if (overlaps(realigned_read, allele)) {
                        assigned_reads.push_back(original_read);
                        allele_realignments.push_back(realigned_read);
                    }
                }
            }
        }
        result.push_back(std::move(allele_realignments));
    }
    const auto reads = boost::get<OverlappingReads::ResultType>(facets.at("OverlappingReads").get());
    const auto overlapping_reads = overlap_range(reads.at(sample), call);
    std::vector<AlignedRead> unassigned_reads {};
    if (assigned_reads.size() < size(overlapping_reads)) {
        std::sort(std::begin(assigned_reads), std::end(assigned_reads));
        unassigned_reads.reserve(size(overlapping_reads) - assigned_reads.size());
        std::set_difference(std::cbegin(overlapping_reads), std::cend(overlapping_reads),
                            std::cbegin(assigned_reads), std::cend(assigned_reads),
                            std::back_inserter(unassigned_reads));
    }
    if (!unassigned_reads.empty()) {
        const auto reference = boost::get<ReferenceContext::ResultType>(facets.at("ReferenceContext").get());
        const auto expanded_reference = expand_for_realignment(reference, unassigned_reads);
        auto unassigned_realignments = realign(unassigned_reads, expanded_reference);
        result.push_back(std::move(unassigned_realignments));
    } else {
        result.push_back(std::move(unassigned_reads));
    }
    return boost::any {compute_realignment_summary(result, call)};
}

std::string Realignments::do_name() const
{
    return "RA";
}

std::vector<std::string> Realignments::do_requirements() const
{
    return {"ReferenceContext", "OverlappingReads", "ReadAssignments"};
}

void serialise_reads(const std::vector<AlignedRead>& reads, std::ostringstream& ss)
{
    ss << "{";
    for (const auto& read : reads) {
        ss << "[" << read.name() << "|" << mapped_begin(read) << "|" << read.cigar() << "]";
    }
    ss << "}";
}

namespace {

std::ostream& operator<<(std::ostream& os, const DirectionCounts& counts)
{
    os << counts.forward << "," << counts.reverse;
    return os;
}

std::ostream& operator<<(std::ostream& os, const MQHistogram& hist)
{
    std::copy(std::cbegin(hist), std::prev(std::cend(hist)), std::ostream_iterator<unsigned> {os, ","});
    os << hist.back();
    return os;
}

std::ostream& operator<<(std::ostream& os, const Pileup& pileup)
{
    std::copy(std::cbegin(pileup.match_counts), std::cend(pileup.match_counts), std::ostream_iterator<unsigned> {os, ","});
    std::copy(std::cbegin(pileup.mismatch_counts), std::cend(pileup.mismatch_counts), std::ostream_iterator<unsigned> {os, ","});
    os << pileup.insertion_count << "," << pileup.deletion_count << ","
       << pileup.match_quality_sum << "," << pileup.mismatch_quality_sum << "," << pileup.insertion_quality_sum;
    return os;
}

std::ostream& operator<<(std::ostream& os, const PileupWindow& window)
{
    std::copy(std::cbegin(window), std::prev(std::cend(window)), std::ostream_iterator<Pileup> {os, ":"});
    os << window.back();
    return os;
}

std::ostream& operator<<(std::ostream& os, const HaplotypeSummary& summary)
{
    os << "{" << summary.strand_counts << "|" << summary.mq_hist << "|" << summary.pileups << "}";
    return os;
}

} // namespace

std::string Realignments::do_serialise(const ResultType& value) const
{
    const auto any = boost::get<boost::any>(value);
    const auto& summary = *boost::any_cast<RealignmentSummary>(&any);
    std::ostringstream ss {};
    std::copy(std::cbegin(summary), std::cend(summary), std::ostream_iterator<HaplotypeSummary> {ss});
    return ss.str();
}

} // namespace csr
} // namespace octopus
