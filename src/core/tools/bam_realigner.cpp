// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "bam_realigner.hpp"

#include <deque>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <utility>
#include <thread>
#include <cmath>
#include <cassert>

#include <boost/functional/hash.hpp>

#include "basics/genomic_region.hpp"
#include "basics/cigar_string.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_header.hpp"
#include "io/read/buffered_read_writer.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/append.hpp"
#include "utils/read_stats.hpp"
#include "utils/random_select.hpp"
#include "utils/maths.hpp"
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

auto copy_reads(AmbiguousReadList&& reads)
{
    std::vector<AlignedRead> result {};
    result.reserve(reads.size());
    std::transform(std::make_move_iterator(std::begin(reads)), std::make_move_iterator(std::end(reads)), std::back_inserter(result),
                   [] (auto&& read) { return std::move(read.read); });
    return result;
}

struct Alignment
{
    GenomicRegion region;
    CigarString cigar;
};

auto copy_alignments(const std::vector<AlignedRead>& reads)
{
    std::vector<Alignment> result {};
    result.reserve(reads.size());
    for (const auto& read : reads) result.push_back({read.mapped_region(), read.cigar()});
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

template <typename Container>
auto copy(const Container& sequence, std::size_t pos, std::size_t length)
{
    assert(pos + length <= sequence.size());
    auto begin_itr = std::next(std::cbegin(sequence), pos);
    return Container {begin_itr, std::next(begin_itr, length)};
}

Haplotype get_aligned_part(const Haplotype& inferred_haplotype, const AlignedRead& realigned_read,
                           const Alignment& inferred_alignment, const ReferenceGenome& reference)
{
    assert(contains(inferred_haplotype, realigned_read));
    const auto& sequence = inferred_haplotype.sequence();
    const auto alignment_start = begin_distance(inferred_haplotype, inferred_alignment.region);
    const auto alignment_length = reference_size(inferred_alignment.cigar);
    return Haplotype {mapped_region(realigned_read), copy(sequence, alignment_start, alignment_length), reference};
}

std::string to_string(const std::vector<int>& ids)
{
    assert(!ids.empty());
    auto result = std::to_string(ids.front());
    if (ids.size() > 1) {
        std::for_each(std::next(std::cbegin(ids)), std::cend(ids), [&result] (auto id) {
            result += "," + std::to_string(id);
        });
    }
    return result;
}

auto realign_and_annotate(const std::vector<AlignedRead>& reads,
                          const Haplotype& haplotype,
                          const ReferenceGenome& reference,
                          const HaplotypeLikelihoodModel& alignment_model,
                          std::vector<int> haplotype_ids = {})
{
    std::vector<AlignedRead> result {};
    if (reads.empty()) return result;
    const auto expanded_haplotype = expand_for_realignment(haplotype, reads, alignment_model);
    std::vector<HaplotypeLikelihoodModel::LogProbability> log_likelihoods {};
    auto realignments = realign(reads, expanded_haplotype, alignment_model, log_likelihoods);
    const auto inferred_alignments = copy_alignments(realignments);
    rebase(realignments, haplotype);
    result.reserve(realignments.size());
    for (std::size_t n {0}; n < realignments.size(); ++n) {
        result.emplace_back(std::move(realignments[n]));
        auto& read = result.back();
        const Haplotype reference_haplotype {mapped_region(read), reference};
        result.back().add_annotation({'M','D'}, to_md_string(read.cigar(), reference_haplotype));
        result.back().add_annotation({'h','c'}, to_string(inferred_alignments[n].cigar));
        const auto inferred_haplotype = get_aligned_part(expanded_haplotype, read, inferred_alignments[n], reference);
        result.back().add_annotation({'m','d'}, to_md_string(inferred_alignments[n].cigar, inferred_haplotype));
        if (!haplotype_ids.empty()) {
            result.back().add_annotation({'H','P'}, to_string(haplotype_ids));
        }
        result.back().add_annotation({'P','S'}, to_string(mapped_region(haplotype)));
        result.back().add_annotation({'L','K'}, std::to_string(static_cast<unsigned>(std::abs(log_likelihoods[n] / maths::constants::ln10Div10<>)))); // std::abs to avoid -0.0
    }
    return result;
}

std::size_t count_reads(const std::vector<AlignedTemplate>& templates) noexcept
{
    return std::accumulate(std::cbegin(templates), std::cend(templates), std::size_t {0},
                           [] (auto curr, const auto& reads) { return curr + reads.size(); });
}

HaplotypeSupportMap
compute_haplotype_support_helper(const Genotype<Haplotype>& genotype,
                                 const std::vector<AlignedTemplate>& templates,
                                 AmbiguousReadList& unassigned_reads,
                                 const HaplotypeLikelihoodModel& alignment_model)
{
    AssignmentConfig assigner_config {};
    assigner_config.ambiguous_record = AssignmentConfig::AmbiguousRecord::haplotypes;
    AmbiguousTemplateList unassigned_templates {};
    const auto template_support = compute_haplotype_support(genotype, templates, unassigned_templates, alignment_model, assigner_config);
    HaplotypeSupportMap result {};
    result.reserve(template_support.size());
    for (const auto& p : template_support) {
        result[p.first].reserve(count_reads(p.second));
        for (const auto& template_reads : p.second) {
            std::copy(std::cbegin(template_reads), std::cend(template_reads), std::back_inserter(result[p.first]));
        }
    }
    for (const auto& unassigned_template : unassigned_templates) {
        for (const auto& read : unassigned_template.read_template) {
            assert(unassigned_template.haplotypes);
            unassigned_reads.emplace_back(read, *unassigned_template.haplotypes);
        }
    }
    return result;
}

auto make_read_templates(const std::vector<AlignedRead>& reads, const ReadLinkageConfig& read_linkage)
{
    std::vector<AlignedTemplate> result {};
    make_read_templates(std::cbegin(reads), std::cend(reads), std::back_inserter(result), read_linkage);
    return result;
}

HaplotypeSupportMap
compute_haplotype_support_helper(const Genotype<Haplotype>& genotype,
                                 const std::vector<AlignedRead>& reads,
                                 AmbiguousReadList& unassigned_reads,
                                 const HaplotypeLikelihoodModel& alignment_model,
                                 const ReadLinkageConfig& read_linkage)
{
    if (read_linkage.linkage != ReadLinkageType::none) {
        const auto templates = make_read_templates(reads, read_linkage);
        return compute_haplotype_support_helper(genotype, templates, unassigned_reads, alignment_model);
    } else {
        AssignmentConfig assigner_config {};
        assigner_config.ambiguous_record = AssignmentConfig::AmbiguousRecord::haplotypes;
        return compute_haplotype_support(genotype, reads, unassigned_reads, alignment_model, assigner_config);
    }
}

auto make_haplotype_id_map(const Genotype<Haplotype>& genotype)
{
    std::unordered_map<Haplotype, std::vector<int>> result {};
    result.reserve(genotype.ploidy());
    int id {0};
    for (const auto& haplotype : genotype) {
        result[haplotype].push_back(id++);
    }
    return result;
}

struct HaplotypeIDHasher
{
    auto operator()(const std::vector<int>& ids) const noexcept
    {
        return boost::hash_range(std::cbegin(ids), std::cend(ids));
    }
};

auto assign_and_realign(const std::vector<AlignedRead>& reads,
                        const Genotype<Haplotype>& genotype,
                        const ReferenceGenome& reference,
                        const HaplotypeLikelihoodModel& alignment_model,
                        const ReadLinkageConfig& read_linkage,
                        BAMRealigner::Report& report)
{
    std::vector<AlignedRead> result {};
    if (!reads.empty()) {
        result.reserve(reads.size());
        if (is_homozygous(genotype)) {
            std::vector<int> amiguous_haplotypes_ids(genotype.ploidy());
            std::iota(std::begin(amiguous_haplotypes_ids), std::end(amiguous_haplotypes_ids), 0);
            utils::append(realign_and_annotate(reads, genotype[0], reference, alignment_model, amiguous_haplotypes_ids), result);
        } else {
            AmbiguousReadList unassigned_reads {};
            auto support = compute_haplotype_support_helper(genotype, reads, unassigned_reads, alignment_model, read_linkage);
            const auto haplotype_ids = make_haplotype_id_map(genotype);
            for (auto& p : support) {
                auto& assigned = p.second;
                if (!assigned.empty()) {
                    const auto& haplotype = p.first;
                    assert(haplotype_ids.count(haplotype) == 1 && ! haplotype_ids.at(haplotype).empty());
                    report.n_reads_assigned += assigned.size();
                    utils::append(realign_and_annotate(assigned, p.first, reference, alignment_model, haplotype_ids.at(haplotype)), result);
                }
            }
            if (!unassigned_reads.empty()) {
                std::unordered_map<std::vector<int>, std::vector<std::size_t>, HaplotypeIDHasher> ambiguous_indices {};
                ambiguous_indices.reserve(2 * genotype.ploidy());
                for (std::size_t idx {0}; idx < unassigned_reads.size(); ++idx) {
                    const AmbiguousRead& ambiguous = unassigned_reads[idx];
                    assert(ambiguous.haplotypes && !ambiguous.haplotypes->empty());
                    std::vector<int> possible_haplotype_ids {};
                    possible_haplotype_ids.reserve(genotype.ploidy() - 1);
                    for (const auto& haplotype : *ambiguous.haplotypes) {
                        utils::append(haplotype_ids.at(*haplotype), possible_haplotype_ids);
                    }
                    std::sort(std::begin(possible_haplotype_ids), std::end(possible_haplotype_ids));
                    ambiguous_indices[possible_haplotype_ids].push_back(idx);
                }
                for (const auto& ambiguous_pair : ambiguous_indices) {
                    const auto& possible_haplotype_ids = ambiguous_pair.first;
                    const auto& read_indices = ambiguous_pair.second;
                    // Reads that could not be assigned to a unique haplotype are randomly aligned to any of the
                    // ambiguous haplotypes for realignment
                    std::unordered_map<Haplotype, std::vector<AlignedRead>> random_assigned_reads {};
                    random_assigned_reads.reserve(possible_haplotype_ids.size());
                    for (const auto& read_idx : read_indices) {
                        const auto& haplotype = genotype[random_select(possible_haplotype_ids)];
                        random_assigned_reads[haplotype].push_back(std::move(unassigned_reads[read_idx].read));
                    }
                    for (auto& p : random_assigned_reads) {
                        utils::append(realign_and_annotate(std::move(p.second), p.first, reference, alignment_model, possible_haplotype_ids), result);
                    }
                }
                report.n_reads_assigned += unassigned_reads.size();
            }
        }
        std::sort(std::begin(result), std::end(result));
    }
    return result;
}

auto to_annotated(std::vector<AlignedRead> reads)
{
    std::vector<AlignedRead> result {};
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
    io::BufferedReadWriter<AlignedRead>::Config writer_config {};
    writer_config.max_buffer_footprint = config_.max_buffer;
    io::BufferedReadWriter<AlignedRead> writer {dst, writer_config};
    Report report {};
    BatchList batch {};
    boost::optional<GenomicRegion> batch_region {};
    for (auto p = variants.iterate(); p.first != p.second;) {
        std::tie(batch, batch_region) = read_next_batch(p.first, p.second, src, reference, samples, batch_region);
        for (auto& sample : batch) {
            std::vector<AlignedRead> genotype_reads {};
            std::vector<AlignedRead> realigned_reads {};
            auto sample_reads_itr = std::begin(sample.reads);
            for (const auto& genotype : sample.genotypes) {
                const auto padded_genotype_region = expand(mapped_region(genotype), 1);
                const auto overlapped_reads = bases(overlap_range(sample_reads_itr, std::end(sample.reads), padded_genotype_region));
                genotype_reads.assign(std::make_move_iterator(overlapped_reads.begin()),
                                      std::make_move_iterator(overlapped_reads.end()));
                sample_reads_itr = sample.reads.erase(overlapped_reads.begin(), overlapped_reads.end());
                auto bad_reads = to_annotated(remove_unalignable_reads(genotype_reads));
                auto realignments = assign_and_realign(genotype_reads, genotype, reference, config_.alignment_model, config_.read_linkage, report);
                report.n_reads_unassigned += bad_reads.size();
                move_merge(bad_reads, realignments);
                move_merge(realignments, realigned_reads);
            }
            move_merge(to_annotated(std::move(sample.reads)), realigned_reads);
            writer << realigned_reads;
        }
    }
    return report;
}

BAMRealigner::Report BAMRealigner::realign(ReadReader& src, VcfReader& variants, ReadWriter& dst,
                                           const ReferenceGenome& reference) const
{
    return realign(src, variants, dst, reference, src.extract_samples());
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
        const GenomicRegion::Distance fetch_pad {config_.copy_hom_ref_reads ? 0 : 10}; // Pad for soft clipped reads
        auto reads = src.fetch_reads(samples, expand(*batch_region, fetch_pad));
        sort(reads);
        if (config_.primary_only) {
            filter_primary(reads);
        }
        boost::optional<GenomicRegion> reads_region {};
        if (has_coverage(reads)) reads_region = encompassing_region(reads);
        for (const auto& sample : samples) {
            auto& sample_genotypes = genotypes[sample];
            if (prev_batch_region) {
                erase_overlapped(reads[sample], expand_rhs(*prev_batch_region, fetch_pad));
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

BAMRealigner::Report
realign(io::ReadReader::Path src, VcfReader::Path variants, io::ReadWriter::Path dst,
        const ReferenceGenome& reference)
{
    return realign(std::move(src), std::move(variants), std::move(dst), reference, BAMRealigner::Config {});
}

BAMRealigner::Report
realign(io::ReadReader::Path src, VcfReader::Path variants, io::ReadWriter::Path dst,
        const ReferenceGenome& reference, BAMRealigner::Config config)
{
    io::ReadWriter dst_bam {std::move(dst), src};
    io::ReadReader src_bam {std::move(src)};
    VcfReader vcf {std::move(variants)};
    BAMRealigner realigner {std::move(config)};
    return realigner.realign(src_bam, vcf, dst_bam, reference);
}

} // namespace octopus
