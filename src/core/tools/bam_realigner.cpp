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
#include "utils/genotype_reader.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/append.hpp"
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
    std::deque<AlignedRead> result {std::make_move_iterator(bad_read_itr),
                                    std::make_move_iterator(std::end(reads))};
    reads.erase(bad_read_itr, std::end(reads));
    return result;
}

bool is_homozygous_nonreference(const Genotype<Haplotype>& genotype)
{
    return genotype.is_homozygous() && !is_reference(genotype[0]);
}

auto vectorise(std::deque<AlignedRead>&& reads)
{
    std::vector<AlignedRead> result {};
    utils::append(std::move(reads), result);
    return result;
}

auto assign_and_realign(const std::vector<AlignedRead>& reads, const Genotype<Haplotype>& genotype,
                        BAMRealigner::Report& report)
{
    std::vector<AlignedRead> result {};
    if (!reads.empty()) {
        result.reserve(reads.size());
        if (is_homozygous_nonreference(genotype)) {
            utils::append(safe_realign_to_reference(reads, genotype[0]), result);
        } else {
            std::deque<AlignedRead> unassigned_reads {};
            auto support = compute_haplotype_support(genotype, reads, unassigned_reads);
            for (auto& p : support) {
                if (!p.second.empty()) {
                    report.n_reads_assigned += p.second.size();
                    safe_realign_to_reference(p.second, p.first);
                    utils::append(std::move(p.second), result);
                }
            }
            if (!unassigned_reads.empty()) {
                report.n_reads_assigned += unassigned_reads.size();
                utils::append(safe_realign_to_reference(vectorise(std::move(unassigned_reads)), genotype[0]), result);
            }
        }
        std::sort(std::begin(result), std::end(result));
    }
    return result;
}

auto merge(std::deque<AlignedRead>& src, std::vector<AlignedRead>& dst)
{
    auto itr = utils::append(std::move(src), dst);
    std::inplace_merge(std::begin(dst), itr, std::end(dst));
}

} // namespace

BAMRealigner::Report BAMRealigner::realign(ReadReader& src, VcfReader& variants, ReadWriter& dst,
                                           const ReferenceGenome& reference, SampleList samples) const
{
    Report report {};
    for (auto p = variants.iterate(); p.first != p.second; ) {
        for (auto sample : read_next_batch(p.first, p.second, src, reference, samples)) {
            for (const auto& genotype : sample.genotypes) {
                auto genotype_reads = copy_overlapped(sample.reads, genotype);
                auto bad_reads = remove_unalignable_reads(genotype_reads);
                auto realignments = assign_and_realign(genotype_reads, genotype, report);
                report.n_reads_unassigned += bad_reads.size();
                merge(bad_reads, realignments);
                dst << realignments;
            }
        }
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
            std::deque<AlignedRead> unassigned_reads {};
            auto support = compute_haplotype_support(genotype, reads, unassigned_reads);
            std::size_t result_idx {0};
            for (auto& p : support) {
                if (!p.second.empty()) {
                    report.n_reads_assigned += p.second.size();
                    safe_realign_to_reference(p.second, p.first);
                    result[result_idx] = std::move(p.second);
                    ++result_idx;
                }
            }
            if (!unassigned_reads.empty()) {
                report.n_reads_assigned += unassigned_reads.size();
                utils::append(std::move(unassigned_reads), result.back());
            }
        }
        for (auto& set : result) std::sort(std::begin(set), std::end(set));
    }
    return result;
}

} // namespace

BAMRealigner::Report BAMRealigner::realign(ReadReader& src, VcfReader& variants, std::vector<ReadWriter>& dsts,
                                           const ReferenceGenome& reference, SampleList samples) const
{
    Report report {};
    for (auto p = variants.iterate(); p.first != p.second; ) {
        for (auto sample : read_next_batch(p.first, p.second, src, reference, samples)) {
            for (const auto& genotype : sample.genotypes) {
                auto genotype_reads = copy_overlapped(sample.reads, genotype);
                auto bad_reads = remove_unalignable_reads(genotype_reads);
                auto realignments = split_and_realign(genotype_reads, genotype, report);
                report.n_reads_unassigned += bad_reads.size();
                merge(bad_reads, realignments.back());
                assert(realignments.size() <= dsts.size());
                for (unsigned i {0}; i < realignments.size() - 1; ++i) {
                    dsts[i] << realignments[i];
                }
                dsts.back() << realignments.back(); // end is always unassigned, but ploidy can change
            }
        }
    }
    return report;
}

BAMRealigner::Report BAMRealigner::realign(ReadReader& src, VcfReader& variants, std::vector<ReadWriter>& dsts,
                                           const ReferenceGenome& reference) const
{
    return realign(src, variants, dsts, reference, src.extract_samples());
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

BAMRealigner::BatchList
BAMRealigner::read_next_batch(VcfIterator& first, const VcfIterator& last, ReadReader& src,
                              const ReferenceGenome& reference, const SampleList& samples) const
{
    const auto records = read_next_block(first, last, samples);
    auto genotypes = extract_genotypes(records, samples, reference);
    BatchList result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        auto& sample_genotypes = genotypes[sample];
        const auto sample_genotype_region = encompassing_region(sample_genotypes);
        auto sample_batch_reads = src.fetch_reads(sample, sample_genotype_region);
        result.push_back({std::move(sample_genotypes), std::move(sample_batch_reads)});
    }
    return result;
}

// non-member methods

BAMRealigner::Report realign(io::ReadReader::Path src, VcfReader::Path variants, io::ReadWriter::Path dst,
                             const ReferenceGenome& reference)
{
    io::ReadWriter dst_bam {std::move(dst), src};
    io::ReadReader src_bam {std::move(src)};
    VcfReader vcf {std::move(variants)};
    BAMRealigner realigner {};
    return realigner.realign(src_bam, vcf, dst_bam, reference);
}

BAMRealigner::Report realign(io::ReadReader::Path src, VcfReader::Path variants,
                             std::vector<io::ReadWriter::Path> dsts, const ReferenceGenome& reference)
{
    std::vector<io::ReadWriter> dst_bams {};
    dst_bams.reserve(dsts.size());
    for (auto& dst : dsts) {
        dst_bams.emplace_back(std::move(dst), src);
    }
    io::ReadReader src_bam {std::move(src)};
    VcfReader vcf {std::move(variants)};
    BAMRealigner realigner {};
    return realigner.realign(src_bam, vcf, dst_bams, reference);
}

} // namespace octopus
