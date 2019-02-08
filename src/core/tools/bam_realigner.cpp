// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "bam_realigner.hpp"

#include <deque>
#include <iterator>
#include <algorithm>
#include <utility>
#include <thread>
#include <cassert>

#include "basics/genomic_region.hpp"
#include "basics/cigar_string.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_header.hpp"
#include "io/read/buffered_read_writer.hpp"
#include "io/read/annotated_aligned_read.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/append.hpp"
#include "utils/read_stats.hpp"
#include "read_assigner.hpp"
#include "read_realigner.hpp"

namespace octopus {

namespace {

unsigned get_pool_size(const BAMRealigner::Config& config)
{
    const auto num_cores = std::thread::hardware_concurrency();
    if (config.max_threads) {
        if (*config.max_threads > 1) {
            return num_cores > 0 ? std::min(*config.max_threads, num_cores) : *config.max_threads;
        } else {
            return 0;
        }
    } else {
        return num_cores > 0 ? num_cores : 8;
    }
}

} // namespace

BAMRealigner::BAMRealigner(Config config)
: config_ {std::move(config)}
, workers_ {get_pool_size(config_)}
{}

namespace {

bool is_alignable(const AlignedRead& read) noexcept
{
    return is_valid(read.cigar());
}

auto remove_unalignable_reads(std::vector<AlignedRead>& reads)
{
    const auto bad_read_itr = std::stable_partition(std::begin(reads), std::end(reads), is_alignable);
    std::vector<AlignedRead> result {std::make_move_iterator(bad_read_itr), std::make_move_iterator(std::end(reads))};
    reads.erase(bad_read_itr, std::end(reads));
    return result;
}

bool is_homozygous_nonreference(const Genotype<Haplotype>& genotype)
{
    return genotype.is_homozygous() && !is_reference(genotype[0]);
}

auto copy_reads(AmbiguousReadList&& reads)
{
    std::vector<AlignedRead> result {};
    result.reserve(reads.size());
    std::transform(std::make_move_iterator(std::begin(reads)), std::make_move_iterator(std::end(reads)), std::back_inserter(result),
                   [] (auto&& read) { return std::move(read.read); });
    return result;
}

auto copy_cigars(const std::vector<AlignedRead>& reads)
{
    std::vector<CigarString> result {};
    result.reserve(reads.size());
    for (const auto& read : reads) result.push_back(read.cigar());
    return result;
}

auto to_md_string(const CigarString& cigar, const Haplotype& haplotype)
{
    std::ostringstream ss {};
    auto sequence_itr = std::cbegin(haplotype.sequence());
    std::size_t match_length {0};
    for (const auto& op : cigar) {
        assert(sequence_itr <= std::cend(haplotype.sequence()));
        using Flag = CigarOperation::Flag;
        switch (op.flag()) {
            case Flag::alignmentMatch: // Assumes no mismatches (would need read sequence)
            case Flag::sequenceMatch: {
                match_length += op.size();
                sequence_itr += op.size();
                break;
            }
            case Flag::deletion: {
                if (match_length > 0) {
                    ss << match_length;
                    match_length = 0;
                }
                ss << '^';
            }
            // fall through
            case Flag::substitution: {
                if (match_length > 0) {
                    ss << match_length;
                    match_length = 0;
                }
                std::copy(sequence_itr, sequence_itr + op.size(), std::ostream_iterator<char> {ss});
            }
            // fall through
            case Flag::softClipped:
            case Flag::padding:
                sequence_itr += op.size();
                // fall through
            case Flag::insertion:
            case Flag::hardClipped:
            default: break;
        }
    }
    if (match_length > 0) ss << match_length;
    return ss.str();
}

auto insertion_size_until_span(const CigarString& cigar, int span)
{
    CigarOperation::Size result {0};
    for (const auto& op : cigar) {
        if (is_insertion(op)) {
            result += op.size();
        } else if (is_deletion(op)) {
            if (result > op.size()) {
                result -= op.size();
            } else {
                result = 0;
            }
        } else {
            span -= op.size();
            if (span <= 0) break;
        }
    }
    return result;
}

template <typename Sequence>
void resize_back(Sequence& sequence, std::size_t n)
{
    sequence.erase(std::cbegin(sequence), std::prev(std::cend(sequence), n));
}

Haplotype get_aligned_part(const Haplotype& inferred_haplotype, const AlignedRead& realigned_read,
                           const CigarString& inferred_cigar, const ReferenceGenome& reference)
{
    if (is_insertion(realigned_read.cigar().front())) {
        const auto realigned_insertion_size = realigned_read.cigar().front().size();
        const auto error_insertion_size = std::min(insertion_size_until_span(inferred_cigar, realigned_insertion_size), realigned_insertion_size);
        const auto inferred_insertion_size = realigned_insertion_size - error_insertion_size;
        if (inferred_insertion_size > 0) {
            assert(!begins_before(realigned_read, inferred_haplotype));
            const auto haplotype_head_region = left_overhang_region(inferred_haplotype, realigned_read);
            auto inserted_sequence = inferred_haplotype.sequence(haplotype_head_region);
            auto anchor_sequence = remap(inferred_haplotype, expand_rhs(haplotype_head_region, 1)).sequence();
            resize_back(anchor_sequence, anchor_sequence.size() - inserted_sequence.size());
            resize_back(inserted_sequence, inferred_insertion_size);
            const auto tail_region = right_overhang_region(realigned_read, expand_rhs(haplotype_head_region, 1));
            auto tail_sequence = remap(inferred_haplotype, tail_region).sequence();
            auto aligned_sequence = inserted_sequence + anchor_sequence + tail_sequence;
            return Haplotype {mapped_region(realigned_read), std::move(aligned_sequence), reference};
        }
    }
    return remap(inferred_haplotype, mapped_region(realigned_read));
}

auto realign_and_annotate(const std::vector<AlignedRead>& reads, const Haplotype& haplotype,
                          const ReferenceGenome& reference,
                          boost::optional<int> haplotype_id = boost::none)
{
    auto realignments = safe_realign(reads, haplotype);
    const auto inferred_cigars = copy_cigars(realignments);
    rebase(realignments, haplotype);
    std::vector<AnnotatedAlignedRead> result {};
    result.reserve(realignments.size());
    for (std::size_t n {0}; n < realignments.size(); ++n) {
        result.emplace_back(std::move(realignments[n]));
        const auto& read = result.back().read();
        const Haplotype reference_haplotype {mapped_region(read), reference};
        result.back().annotate("MD", to_md_string(read.cigar(), reference_haplotype));
        result.back().annotate("hc", to_string(inferred_cigars[n]));
        const auto inferred_haplotype = get_aligned_part(haplotype, read, inferred_cigars[n], reference);
        result.back().annotate("md", to_md_string(inferred_cigars[n], inferred_haplotype));
        if (haplotype_id) {
            result.back().annotate("hi", std::to_string(*haplotype_id));
        }
    }
    return result;
}

auto assign_and_realign(const std::vector<AlignedRead>& reads, const Genotype<Haplotype>& genotype,
                        const ReferenceGenome& reference, BAMRealigner::Report& report)
{
    std::vector<AnnotatedAlignedRead> result {};
    if (!reads.empty()) {
        result.reserve(reads.size());
        if (is_homozygous_nonreference(genotype)) {
            utils::append(realign_and_annotate(reads, genotype[0], reference, genotype.ploidy()), result);
        } else {
            AmbiguousReadList unassigned_reads {};
            auto support = compute_haplotype_support(genotype, reads, unassigned_reads);
            int haplotype_id {0};
            for (auto& p : support) {
                if (!p.second.empty()) {
                    report.n_reads_assigned += p.second.size();
                    utils::append(realign_and_annotate(p.second, p.first, reference, haplotype_id), result);
                }
                ++haplotype_id;
            }
            if (!unassigned_reads.empty()) {
                report.n_reads_assigned += unassigned_reads.size();
                utils::append(realign_and_annotate(copy_reads(std::move(unassigned_reads)), genotype[0], reference, genotype.ploidy()), result);
            }
        }
        std::sort(std::begin(result), std::end(result));
    }
    return result;
}

auto to_annotated(std::vector<AlignedRead> reads)
{
    std::vector<AnnotatedAlignedRead> result {};
    result.reserve(reads.size());
    for (auto& read : reads) result.emplace_back(std::move(read));
    return result;
}

template <typename Container, typename T>
auto move_merge(Container&& src, std::vector<T>& dst)
{
    assert(std::is_sorted(std::cbegin(src), std::cend(src)));
    assert(std::is_sorted(std::cbegin(dst), std::cend(dst)));
    auto itr = utils::append(std::move(src), dst);
    std::inplace_merge(std::begin(dst), itr, std::end(dst));
}

} // namespace

BAMRealigner::Report
BAMRealigner::realign(ReadReader& src, VcfReader& variants, ReadWriter& dst,
                      const ReferenceGenome& reference, SampleList samples) const
{
    io::BufferedReadWriter<AnnotatedAlignedRead>::Config writer_config {};
    writer_config.max_buffer_footprint = config_.max_buffer;
    io::BufferedReadWriter<AnnotatedAlignedRead> writer {dst, writer_config};
    Report report {};
    BatchList batch {};
    boost::optional<GenomicRegion> batch_region {};
    for (auto p = variants.iterate(); p.first != p.second;) {
        std::tie(batch, batch_region) = read_next_batch(p.first, p.second, src, reference, samples, batch_region);
        for (auto& sample : batch) {
            std::vector<AlignedRead> genotype_reads {};
            std::vector<AnnotatedAlignedRead> realigned_reads {};
            auto sample_reads_itr = std::begin(sample.reads);
            for (const auto& genotype : sample.genotypes) {
                const auto padded_genotype_region = expand(mapped_region(genotype), 1);
                const auto overlapped_reads = bases(overlap_range(sample_reads_itr, std::end(sample.reads), padded_genotype_region));
                genotype_reads.assign(std::make_move_iterator(overlapped_reads.begin()),
                                      std::make_move_iterator(overlapped_reads.end()));
                sample_reads_itr = sample.reads.erase(overlapped_reads.begin(), overlapped_reads.end());
                auto bad_reads = to_annotated(remove_unalignable_reads(genotype_reads));
                auto realignments = assign_and_realign(genotype_reads, genotype, reference, report);
                report.n_reads_unassigned += bad_reads.size();
                move_merge(bad_reads, realignments);
                move_merge(realignments, realigned_reads);
            }
            move_merge(to_annotated(std::move(sample.reads)), realigned_reads);
            writer << realigned_reads;
        }
        batch_region = encompassing_region(batch.front().genotypes);
    }
    return report;
}

BAMRealigner::Report BAMRealigner::realign(ReadReader& src, VcfReader& variants, ReadWriter& dst,
                                           const ReferenceGenome& reference) const
{
    return realign(src, variants, dst, reference, src.extract_samples());
}

namespace {

auto split_and_realign(const std::vector<AlignedRead>& reads, const Genotype<Haplotype>& genotype,
                       BAMRealigner::Report& report)
{
    std::vector<std::vector<AlignedRead>> result(genotype.zygosity() + 1);
    if (!reads.empty()) {
        if (is_homozygous_nonreference(genotype)) {
            report.n_reads_assigned += reads.size();
            result.back() = safe_realign_to_reference(reads, genotype[0]);
        } else {
            AmbiguousReadList unassigned_reads {};
            auto support = compute_haplotype_support(genotype, reads, unassigned_reads);
            std::size_t result_idx {0};
            for (const auto& haplotype : genotype) {
                auto support_itr = support.find(haplotype);
                if (support_itr != std::cend(support)) {
                    auto& haplotype_support = support_itr->second;
                    if (!haplotype_support.empty()) {
                        report.n_reads_assigned += haplotype_support.size();
                        safe_realign_to_reference(haplotype_support, haplotype);
                        result[result_idx] = std::move(haplotype_support);
                        ++result_idx;
                    }
                    support.erase(support_itr);
                }
            }
            if (!unassigned_reads.empty()) {
                report.n_reads_assigned += unassigned_reads.size();
                move_merge(copy_reads(std::move(unassigned_reads)), result.back());
            }
        }
        for (auto& set : result) std::sort(std::begin(set), std::end(set));
    }
    return result;
}

void move_merge(std::vector<std::vector<AlignedRead>>& src, std::vector<std::vector<AlignedRead>>& dst)
{
    if (dst.empty()) {
        dst = std::move(src);
    } else {
        if (src.size() > dst.size()) dst.resize(src.size());
        for (std::size_t i {0}; i < src.size(); ++i) {
            move_merge(std::move(src[i]), dst[i]);
        }
    }
}

} // namespace

BAMRealigner::Report
BAMRealigner::realign(ReadReader& src, VcfReader& variants, std::vector<ReadWriter>& dsts,
                      const ReferenceGenome& reference, SampleList samples) const
{
    if (dsts.size() == 1) return realign(src, variants, dsts.front(), reference, samples);
    io::BufferedReadWriter<AlignedRead>::Config writer_config {};
    writer_config.max_buffer_footprint = config_.max_buffer.bytes() / dsts.size();
    std::vector<io::BufferedReadWriter<AlignedRead>> writers {};
    writers.reserve(dsts.size());
    for (auto& dst : dsts) writers.emplace_back(dst, writer_config);
    Report report {};
    BatchList batch {};
    boost::optional<GenomicRegion> batch_region {};
    for (auto p = variants.iterate(); p.first != p.second; ) {
        std::tie(batch, batch_region) = read_next_batch(p.first, p.second, src, reference, samples, batch_region);
        for (auto& sample : batch) {
            std::vector<AlignedRead> genotype_reads {}, unassigned_realigned_reads {};
            std::vector<std::vector<AlignedRead>> assigned_realigned_reads {};
            auto sample_reads_itr = std::begin(sample.reads);
            for (const auto& genotype : sample.genotypes) {
                const auto overlapped_reads = bases(overlap_range(sample_reads_itr, std::end(sample.reads), genotype));
                genotype_reads.assign(std::make_move_iterator(overlapped_reads.begin()),
                                      std::make_move_iterator(overlapped_reads.end()));
                sample_reads_itr = sample.reads.erase(overlapped_reads.begin(), overlapped_reads.end());
                auto bad_reads = remove_unalignable_reads(genotype_reads);
                auto realignments = split_and_realign(genotype_reads, genotype, report);
                report.n_reads_unassigned += bad_reads.size();
                move_merge(bad_reads, realignments.back());
                move_merge(realignments.back(), unassigned_realigned_reads); // end is always unassigned, but ploidy can change
                realignments.pop_back();
                move_merge(realignments, assigned_realigned_reads);
            }
            move_merge(unassigned_realigned_reads, sample.reads);
            for (unsigned i {0}; i < assigned_realigned_reads.size(); ++i) {
                writers[i] << assigned_realigned_reads[i];
            }
            writers.back() << sample.reads;
        }
    }
    return report;
}

namespace {

auto get_sample_intersection(const io::ReadReader& bam, const VcfReader& vcf)
{
    auto bam_samples = bam.extract_samples();
    std::sort(std::begin(bam_samples), std::end(bam_samples));
    auto vcf_samples = vcf.fetch_header().samples();
    std::sort(std::begin(vcf_samples), std::end(vcf_samples));
    std::vector<VcfRecord::SampleName> result {};
    result.reserve(std::min(bam_samples.size(), vcf_samples.size()));
    std::set_intersection(std::cbegin(bam_samples), std::cend(bam_samples),
                          std::cbegin(vcf_samples), std::cend(vcf_samples),
                          std::back_inserter(result));
    return result;
}

} // namespace

BAMRealigner::Report BAMRealigner::realign(ReadReader& src, VcfReader& variants, std::vector<ReadWriter>& dsts,
                                           const ReferenceGenome& reference) const
{
    auto samples = get_sample_intersection(src, variants);
    return realign(src, variants, dsts, reference, std::move(samples));
}

// private methods

namespace {

GenomicRegion get_phase_set(const VcfRecord& record, const SampleName& sample)
{
    auto result = get_phase_region(record, sample);
    return result ? *result : mapped_region(record);
}

std::vector<GenomicRegion> get_phase_sets(const VcfRecord& record, const std::vector<SampleName>& samples)
{
    std::vector<GenomicRegion> result{};
    result.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::back_inserter(result),
                   [&record](const auto& sample) { return get_phase_set(record, sample); });
    return result;
}

GenomicRegion get_phase_region(const VcfRecord& record, const std::vector<SampleName>& samples)
{
    return encompassing_region(get_phase_sets(record, samples));
}

template<typename T, typename _>
std::vector<T> copy_each_first(const std::vector<std::pair<T, _>>& items)
{
    std::vector<T> result {};
    result.reserve(items.size());
    std::transform(std::cbegin(items), std::cend(items), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    return result;
}

} // namespace

BAMRealigner::CallBlock
BAMRealigner::read_next_block(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const
{
    std::vector<std::pair<VcfRecord, GenomicRegion>> block {};
    for (; first != last; ++first) {
        const VcfRecord& call {*first};
        auto call_phase_region = get_phase_region(call, samples);
        if (!block.empty() && !overlaps(block.back().second, call_phase_region)) {
            return copy_each_first(block);
        }
        block.emplace_back(call, std::move(call_phase_region));
    }
    return copy_each_first(block);
}

namespace {

void sort(io::ReadReader::SampleReadMap& reads)
{
    for (auto& p : reads) std::sort(std::begin(p.second), std::end(p.second));
}

void filter_primary(std::vector<AlignedRead>& reads)
{
    auto itr = std::remove_if(std::begin(reads), std::end(reads), [&] (const auto& read) { return !is_primary_alignment(read); });
    reads.erase(itr, std::end(reads));
}

void filter_primary(io::ReadReader::SampleReadMap& reads)
{
    for (auto& p : reads) filter_primary(p.second);
}

void erase_overlapped(std::vector<AlignedRead>& reads, const GenomicRegion& region)
{
    auto itr = std::remove_if(std::begin(reads), std::end(reads), [&] (const auto& read) { return overlaps(read, region); });
    reads.erase(itr, std::end(reads));
}

} // namespace

BAMRealigner::BatchListRegionPair
BAMRealigner::read_next_batch(VcfIterator& first, const VcfIterator& last, ReadReader& src,
                              const ReferenceGenome& reference, const SampleList& samples,
                              const boost::optional<GenomicRegion>& prev_batch_region) const
{
    const auto records = read_next_block(first, last, samples);
    BatchList batches {};
    boost::optional<GenomicRegion> batch_region {};
    if (!records.empty()) {
        auto genotypes = extract_genotypes(records, samples, reference);
        batches.reserve(samples.size());
        batch_region = encompassing_region(records);
        if (config_.copy_hom_ref_reads) {
            if (prev_batch_region && is_same_contig(*batch_region, *prev_batch_region)) {
                batch_region = right_overhang_region(*batch_region, *prev_batch_region);
                if (first != last && is_same_contig(*batch_region, *first)) {
                    batch_region = expand_rhs(*batch_region, intervening_region_size(*batch_region, *first) / 2);
                } else {
                    batch_region = closed_region(*batch_region, reference.contig_region(prev_batch_region->contig_name()));
                }
            } else {
                if (first != last && is_same_contig(*batch_region, *first)) {
                    batch_region = expand(*batch_region, batch_region->begin(), intervening_region_size(*batch_region, *first) / 2);
                } else {
                    batch_region = expand_lhs(*batch_region, batch_region->begin());
                }
            }
        }
        auto reads = src.fetch_reads(samples, expand(*batch_region, 1)); // Pad for soft clipped reads
        sort(reads);
        if (config_.primary_only) {
            filter_primary(reads);
        }
        boost::optional<GenomicRegion> reads_region {};
        if (has_coverage(reads)) reads_region = encompassing_region(reads);
        for (const auto& sample : samples) {
            auto& sample_genotypes = genotypes[sample];
            if (prev_batch_region) {
                erase_overlapped(reads[sample], expand_rhs(*prev_batch_region, 1));
            }
            batches.push_back({std::move(sample_genotypes), std::move(reads[sample])});
        }
        if (first != last && reads_region && overlaps(*first, *reads_region)) {
            auto p = read_next_batch(first, last, src, reference, samples, batch_region);
            merge(p.first, batches);
            assert(p.second);
            batch_region = closed_region(*batch_region, *p.second);
        }
    } else if (prev_batch_region && config_.copy_hom_ref_reads) {
        const auto contig_region = reference.contig_region(prev_batch_region->contig_name());
        const auto reads_region = right_overhang_region(contig_region, *prev_batch_region);
        if (!is_empty(reads_region)) {
            auto reads = src.fetch_reads(samples, reads_region);
            for (const auto& sample : samples) {
                erase_overlapped(reads[sample], *prev_batch_region);
                batches.push_back({{}, std::move(reads[sample])});
            }
        }
    }
    return {std::move(batches), std::move(batch_region)};
}

void BAMRealigner::merge(BatchList& src, BatchList& dst) const
{
    assert(src.size() == dst.size());
    for (std::size_t i {0}; i < src.size(); ++i) {
        dst[i].genotypes.insert(std::make_move_iterator(std::begin(src[i].genotypes)),
                                std::make_move_iterator(std::end(src[i].genotypes)));
        move_merge(std::move(src[i].reads), dst[i].reads);
    }
}

// non-member methods

BAMRealigner::Report realign(io::ReadReader::Path src, VcfReader::Path variants, io::ReadWriter::Path dst,
                             const ReferenceGenome& reference)
{
    return realign(std::move(src), std::move(variants), std::move(dst), reference, BAMRealigner::Config {});
}

BAMRealigner::Report realign(io::ReadReader::Path src, VcfReader::Path variants, io::ReadWriter::Path dst,
                             const ReferenceGenome& reference, BAMRealigner::Config config)
{
    io::ReadWriter dst_bam {std::move(dst), src};
    io::ReadReader src_bam {std::move(src)};
    VcfReader vcf {std::move(variants)};
    BAMRealigner realigner {std::move(config)};
    return realigner.realign(src_bam, vcf, dst_bam, reference);
}

BAMRealigner::Report realign(io::ReadReader::Path src, VcfReader::Path variants,
                             std::vector<io::ReadWriter::Path> dsts, const ReferenceGenome& reference)
{
    return realign(std::move(src), std::move(variants), std::move(dsts), reference, BAMRealigner::Config {});
}

BAMRealigner::Report realign(io::ReadReader::Path src, VcfReader::Path variants, std::vector<io::ReadWriter::Path> dsts,
                             const ReferenceGenome& reference, BAMRealigner::Config config)
{
    std::vector<io::ReadWriter> dst_bams {};
    dst_bams.reserve(dsts.size());
    for (auto& dst : dsts) {
        dst_bams.emplace_back(std::move(dst), src);
    }
    io::ReadReader src_bam {std::move(src)};
    VcfReader vcf {std::move(variants)};
    BAMRealigner realigner {std::move(config)};
    return realigner.realign(src_bam, vcf, dst_bams, reference);
}

} // namespace octopus
