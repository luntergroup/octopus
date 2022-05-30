// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_realigner.hpp"

#include <iterator>
#include <algorithm>
#include <utility>
#include <cassert>

#include "basics/genomic_region.hpp"
#include "utils/maths.hpp"
#include "utils/kmer_mapper.hpp"
#include "utils/append.hpp"
#include "utils/parallel_transform.hpp"
#include "core/models/error/error_model_factory.hpp"

namespace octopus {

namespace {

GenomicRegion::Size estimate_max_indel_size(const AlignedRead& read)
{
    const auto p = std::minmax({region_size(read), static_cast<GenomicRegion::Size>(sequence_size(read))});
    return p.second - p.first;
}

auto estimate_max_indel_size(const std::vector<AlignedRead>& reads)
{
    GenomicRegion::Size result {0};
    for (const auto& read : reads) {
        result = std::max(result, estimate_max_indel_size(read));
    }
    return result;
}

auto safe_expand(const GenomicRegion& region, const GenomicRegion::Size n)
{
    if (region.begin() < n) {
        return GenomicRegion {region.contig_name(), 0, region.end() + 2 * n - region.begin()};
    } else {
        return expand(region, n);
    }
}

HaplotypeLikelihoodModel make_default_haplotype_likelihood_model()
{
    HaplotypeLikelihoodModel::Config config {};
    config.max_indel_error = 8;
    config.use_flank_state = false;
    config.use_mapping_quality = false;
    return {config};
}

} // namespace

Haplotype expand_for_realignment(const Haplotype& haplotype, const std::vector<AlignedRead>& reads,
                                 const HaplotypeLikelihoodModel& model)
{
    const auto reads_region = encompassing_region(reads);
    auto required_flank_pad = 2 * model.pad_requirement();
    if (region_size(haplotype) > sequence_size(haplotype)) {
        required_flank_pad += region_size(haplotype) - sequence_size(haplotype);
    }
    required_flank_pad += estimate_max_indel_size(reads);
    const auto required_haplotype_region = safe_expand(reads_region, required_flank_pad);
    if (contains(haplotype, required_haplotype_region)) {
        return haplotype;
    } else {
        return remap(haplotype, encompassing_region(haplotype, required_haplotype_region));
    }
}

Haplotype expand_for_realignment(const Haplotype& haplotype, const std::vector<AlignedRead>& reads)
{
    return expand_for_realignment(haplotype, reads, make_default_haplotype_likelihood_model());
}

namespace {

auto compute_read_hashes(const std::vector<AlignedRead>& reads, boost::optional<ThreadPool&> workers = boost::none)
{
    static constexpr unsigned char mapperKmerSize {6};
    std::vector<KmerPerfectHashes> result(reads.size());
    using octopus::transform;
    transform(std::cbegin(reads), std::cend(reads), std::begin(result),
              [=] (const AlignedRead& read) { return compute_kmer_hashes<mapperKmerSize>(read.sequence()); },
              workers);
    return result;
}

void realign(AlignedRead& read, const Haplotype& haplotype,
             GenomicRegion::Position mapping_position, CigarString alignment)
{
    const auto remapped_read_begin = mapped_begin(haplotype) + mapping_position;
    const auto remapped_read_end = remapped_read_begin + reference_size(alignment);
    read.realign(GenomicRegion {contig_name(read), remapped_read_begin, remapped_read_end}, std::move(alignment));
}

void realign(AlignedRead& read, const Haplotype& haplotype, HaplotypeLikelihoodModel::Alignment alignment)
{
    realign(read, haplotype, alignment.mapping_position, std::move(alignment.cigar));
}

void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype,
             HaplotypeLikelihoodModel model,
             boost::optional<std::vector<HaplotypeLikelihoodModel::LogProbability>&> log_likelihoods,
             boost::optional<ThreadPool&> workers)
{
    if (!reads.empty()) {
        const auto read_hashes = compute_read_hashes(reads, workers);
        static constexpr unsigned char mapperKmerSize {6};
        const auto haplotype_hashes = make_kmer_hash_table<mapperKmerSize>(haplotype.sequence());
        model.reset(haplotype);
        auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
        if (workers && workers->size() > 1 && reads.size() > 1) {
            std::vector<std::future<void>> futures(reads.size());
            for (std::size_t i {0}; i < reads.size(); ++i) {
                workers->try_push([
                    &haplotype,
                    &haplotype_hashes, 
                    &read = reads[i], 
                    &read_hashes = read_hashes[i], 
                    haplotype_mapping_counts,
                    &model,
                    &log_likelihoods,
                    i] () mutable {
                    auto mapping_positions = map_query_to_target(read_hashes, haplotype_hashes, haplotype_mapping_counts);
                    auto alignment = model.align(read, mapping_positions);
                    if (log_likelihoods) (*log_likelihoods)[i] = alignment.likelihood;
                    realign(read, haplotype, std::move(alignment));
                });
            }
            for (auto& f : futures) f.get();
        } else {
            for (std::size_t i {0}; i < reads.size(); ++i) {
                auto mapping_positions = map_query_to_target(read_hashes[i], haplotype_hashes, haplotype_mapping_counts);
                reset_mapping_counts(haplotype_mapping_counts);
                auto alignment = model.align(reads[i], mapping_positions);
                if (log_likelihoods) (*log_likelihoods)[i] = alignment.likelihood;
                realign(reads[i], haplotype, std::move(alignment));
            }
        }
    }
}

} // namespace

void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype,
             HaplotypeLikelihoodModel model,
             std::vector<HaplotypeLikelihoodModel::LogProbability>& log_likelihoods,
             boost::optional<ThreadPool&> workers)
{
    log_likelihoods.resize(reads.size());
    boost::optional<std::vector<HaplotypeLikelihoodModel::LogProbability>&> ll {log_likelihoods};
    realign(reads, haplotype, std::move(model), ll, workers);
}

void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype, 
            HaplotypeLikelihoodModel model, 
            boost::optional<ThreadPool&> workers)
{
    realign(reads, haplotype, std::move(model), boost::none, workers);
}

void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype, boost::optional<ThreadPool&> workers)
{
    realign(reads, haplotype, make_default_haplotype_likelihood_model(), workers);
}

std::vector<AlignedRead> 
realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype,
        HaplotypeLikelihoodModel model, boost::optional<ThreadPool&> workers)
{
    auto result = reads;
    realign(result, haplotype, model, workers);
    return result;
}

std::vector<AlignedRead> realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, boost::optional<ThreadPool&> workers)
{
    return realign(reads, haplotype, make_default_haplotype_likelihood_model(), workers);
}

std::vector<AlignedRead>
realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype,
        HaplotypeLikelihoodModel model,
        std::vector<HaplotypeLikelihoodModel::LogProbability>& log_likelihoods,
        boost::optional<ThreadPool&> workers)
{
    auto result = reads;
    realign(result, haplotype, model, log_likelihoods, workers);
    return result;
}

void safe_realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model, boost::optional<ThreadPool&> workers)
{
    if (!reads.empty()) {
        auto expanded_haplotype = expand_for_realignment(haplotype, reads, model);
        try {
            realign(reads, expanded_haplotype, model, workers);
        } catch (const HaplotypeLikelihoodModel::ShortHaplotypeError& e) {
            expanded_haplotype = expand(expanded_haplotype, e.required_extension());
            realign(reads, expanded_haplotype, model, workers);
        }
    }
}

void safe_realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model,
                  std::vector<HaplotypeLikelihoodModel::LogProbability>& log_likelihoods,
                  boost::optional<ThreadPool&> workers)
{
    if (!reads.empty()) {
        auto expanded_haplotype = expand_for_realignment(haplotype, reads, model);
        try {
            realign(reads, expanded_haplotype, model, log_likelihoods, workers);
        } catch (const HaplotypeLikelihoodModel::ShortHaplotypeError& e) {
            log_likelihoods.clear();
            expanded_haplotype = expand(expanded_haplotype, e.required_extension());
            realign(reads, expanded_haplotype, model, log_likelihoods, workers);
        }
    }
}

void safe_realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype, boost::optional<ThreadPool&> workers)
{
    safe_realign(reads, haplotype, make_default_haplotype_likelihood_model(), workers);
}

std::vector<AlignedRead>
safe_realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model,
             boost::optional<ThreadPool&> workers)
{
    auto result = reads;
    safe_realign(result, haplotype, std::move(model), workers);
    return result;
}

std::vector<AlignedRead> 
safe_realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, boost::optional<ThreadPool&> workers)
{
    return safe_realign(reads, haplotype, make_default_haplotype_likelihood_model(), workers);
}

namespace {

CigarString minimise(const CigarString& cigar)
{
    CigarString result {};
    result.reserve(cigar.size());
    for (auto op_itr = std::cbegin(cigar); op_itr != std::cend(cigar);) {
        const auto next_op_itr = std::find_if_not(std::next(op_itr), std::cend(cigar),
                                                  [=] (const auto& op) { return op.flag() == op_itr->flag();});
        const auto op_size = std::accumulate(op_itr, next_op_itr, CigarOperation::Size {0},
                                             [] (auto curr, const auto& op) { return curr + op.size(); });
        if (op_size > 0) {
            result.emplace_back(op_size, op_itr->flag());
        }
        op_itr = next_op_itr;
    }
    return result;
}

auto get_match_type(const CigarOperation::Flag haplotype, const CigarOperation::Flag read) noexcept
{
    assert(is_match_or_substitution(haplotype) && is_match_or_substitution(read));
    using Flag = CigarOperation::Flag;
    if ((haplotype == Flag::substitution && read == Flag::sequenceMatch)
        || (haplotype == Flag::sequenceMatch && read == Flag::substitution)) {
        return Flag::substitution;
    } else if (haplotype == Flag::sequenceMatch && read == Flag::sequenceMatch) {
        return Flag::sequenceMatch;
    } else {
        return Flag::alignmentMatch;
    }
}

} // namespace

CigarString rebase(const CigarString& read_to_haplotype, const CigarString& haplotype_to_reference)
{
    assert(is_valid(read_to_haplotype) && is_valid(haplotype_to_reference));
    assert(reference_size(read_to_haplotype) <= sequence_size(haplotype_to_reference));
    const auto haplotypes_ops = decompose(haplotype_to_reference);
    CigarString result {};
    result.reserve(haplotypes_ops.size());
    auto hap_flag_itr = std::cbegin(haplotypes_ops);
    for (const auto& read_op : read_to_haplotype) {
        if (is_match_or_substitution(read_op)) {
            for (unsigned n {0}; n < read_op.size();) {
                assert(hap_flag_itr != std::cend(haplotypes_ops));
                if (is_match_or_substitution(*hap_flag_itr)) {
                    result.emplace_back(1, get_match_type(*hap_flag_itr, read_op.flag()));
                    ++n;
                } else {
                    result.emplace_back(1, *hap_flag_itr);
                    if (advances_sequence(*hap_flag_itr)) {
                        ++n;
                    }
                }
                ++hap_flag_itr;
            }
        } else if (is_insertion(read_op)) {
            result.push_back(read_op);
        } else { // deletion
            auto op_size = read_op.size();
            for (unsigned n {0}; n < read_op.size();) {
                assert(hap_flag_itr != std::cend(haplotypes_ops));
                if (is_deletion(*hap_flag_itr)) {
                    result.emplace_back(1, *hap_flag_itr);
                    ++hap_flag_itr;
                } else {
                    if (is_insertion(*hap_flag_itr)) {
                        --op_size;
                    }
                    ++hap_flag_itr;
                    ++n;
                }
            }
            if (op_size > 0) {
                result.emplace_back(op_size, read_op.flag());
            }
        }
    }
    return minimise(result);
}

namespace {

bool is_sequence_match(const CigarOperation& op) noexcept
{
    return op.flag() == CigarOperation::Flag::sequenceMatch;
}

CigarString pad_reference(const GenomicRegion& read_region, const CigarString& read_to_haplotype,
                          const GenomicRegion& haplotype_region, const CigarString& haplotype_to_reference)
{
    assert(overlaps(read_region, haplotype_region));
    assert(!read_to_haplotype.empty() || !haplotype_to_reference.empty());
    CigarString result {};
    using Flag = CigarOperation::Flag;
    if (read_region == haplotype_region) {
        result = haplotype_to_reference;
    } else if (contains(haplotype_region, read_region)) {
        const auto offset = left_overhang_size(haplotype_region, read_region);
        const auto copy_length = std::max(sequence_size(read_to_haplotype), sequence_size(haplotype_to_reference));
        result = copy_sequence(haplotype_to_reference, offset, copy_length);
    } else {
        result.reserve(haplotype_to_reference.size() + 2);
        if (contains(read_region, haplotype_region)) {
            const auto lhs_pad_size = left_overhang_size(read_region, haplotype_region);
            const auto rhs_pad_size = right_overhang_size(read_region, haplotype_region);
            if (is_sequence_match(haplotype_to_reference.front())) {
                if (haplotype_to_reference.size() == 1) {
                    result.emplace_back(lhs_pad_size + haplotype_to_reference.front().size() + rhs_pad_size, Flag::sequenceMatch);
                } else {
                    result.emplace_back(lhs_pad_size + haplotype_to_reference.front().size(), Flag::sequenceMatch);
                    result.insert(result.cend(), haplotype_to_reference.cbegin() + 1, haplotype_to_reference.cend() - 1);
                    if (is_sequence_match(haplotype_to_reference.back())) {
                        result.emplace_back(haplotype_to_reference.back().size() + rhs_pad_size, Flag::sequenceMatch);
                    } else {
                        result.push_back(haplotype_to_reference.back());
                        if (rhs_pad_size > 0) result.emplace_back(rhs_pad_size, Flag::sequenceMatch);
                    }
                }
            } else {
                if (lhs_pad_size > 0) result.emplace_back(lhs_pad_size, Flag::sequenceMatch);
                if (is_sequence_match(haplotype_to_reference.back())) {
                    result.insert(result.cend(), haplotype_to_reference.cbegin(), haplotype_to_reference.cend() - 1);
                    result.emplace_back(haplotype_to_reference.back().size() + rhs_pad_size, Flag::sequenceMatch);
                } else {
                    utils::append(haplotype_to_reference, result);
                    if (rhs_pad_size > 0) result.emplace_back(rhs_pad_size, Flag::sequenceMatch);
                }
            }
        } else if (begins_before(read_region, haplotype_region)) {
            assert(ends_before(read_region, haplotype_region));
            const auto lhs_pad_size = left_overhang_size(read_region, haplotype_region);
            if (is_sequence_match(haplotype_to_reference.front())) {
                result.emplace_back(lhs_pad_size + haplotype_to_reference.front().size(), Flag::sequenceMatch);
                result.insert(result.cend(), haplotype_to_reference.cbegin() + 1, haplotype_to_reference.cend());
            } else {
                assert(lhs_pad_size > 0);
                result.emplace_back(lhs_pad_size, Flag::sequenceMatch);
                utils::append(haplotype_to_reference, result);
            }
        } else {
            assert(begins_before(haplotype_region, read_region) && ends_before(haplotype_region, read_region));
            const auto offset = left_overhang_size(haplotype_region, read_region);
            const auto rhs_pad_size = right_overhang_size(read_region, haplotype_region);
            result = copy_sequence(haplotype_to_reference, offset);
            if (!result.empty() && is_sequence_match(result.back())) {
                increment_size(result.back(), rhs_pad_size);
            } else if (rhs_pad_size > 0) {
                result.emplace_back(rhs_pad_size, Flag::sequenceMatch);
            }
        }
    }
    if (sequence_size(result) < reference_size(read_to_haplotype)) {
        const auto rhs_pad_size = reference_size(read_to_haplotype) - sequence_size(result);
        if (!result.empty() && is_sequence_match(result.back())) {
            increment_size(result.back(), rhs_pad_size);
        } else {
            assert(rhs_pad_size > 0);
            result.emplace_back(rhs_pad_size, Flag::sequenceMatch);
        }
    }
    return result;
}

CigarString copy_tail(const GenomicRegion& haplotype_region, const CigarString& haplotype_to_reference,
                      const GenomicRegion& rebased_read_region)
{
    const auto offset = left_overhang_size(haplotype_region, rebased_read_region);
    auto result = copy_reference(haplotype_to_reference, offset, size(haplotype_region));
    if (sequence_size(result) < size(rebased_read_region)) {
        const auto rhs_pad_size = size(rebased_read_region) - sequence_size(result);
        if (is_sequence_match(result.back())) {
            increment_size(result.back(), rhs_pad_size);
        } else {
            result.emplace_back(rhs_pad_size, CigarOperation::Flag::sequenceMatch);
        }
    }
    return result;
}

auto get_tail_op_size(const GenomicRegion& haplotype_region, const CigarString& haplotype_to_reference,
                      const GenomicRegion& read_region)
{
    auto tail = copy_sequence(haplotype_to_reference, left_overhang_size(haplotype_region, read_region),
                              haplotype_to_reference.back().size());
    return tail.empty() ? 0 : tail.back().size();
}

bool has_indels(const CigarString& cigar) noexcept
{
    return std::any_of(std::cbegin(cigar), std::cend(cigar), [] (const auto& op) { return is_indel(op); });
}

auto calculate_rebase_shift(const AlignedRead& read, const GenomicRegion& haplotype_region,
                            const CigarString& haplotype_to_reference)
{
    GenomicRegion::Distance result {0};
    if (begins_before(haplotype_region, read) && has_indels(haplotype_to_reference)) {
        auto lhs_flank_length = static_cast<int>(left_overhang_size(haplotype_region, read));
        for (const auto& op : haplotype_to_reference) {
            if (lhs_flank_length == 0) {
                if (is_deletion(op)) {
                    result += op.size();
                }
                break;
            } else if (lhs_flank_length < 0) {
                break;
            }
            if (is_insertion(op)) {
                result -= std::min(static_cast<int>(op.size()), lhs_flank_length);
                lhs_flank_length -= op.size();
            } else if (is_deletion(op)) {
                result += op.size();
            } else {
                lhs_flank_length -= op.size();
            }
        }
    }
    return result;
}

void rebase_overlapped(AlignedRead& read, const GenomicRegion& haplotype_region, const CigarString& haplotype_to_reference,
                       const GenomicRegion::Distance rebase_shift)
{
    auto padded_haplotype_cigar = pad_reference(read.mapped_region(), read.cigar(), haplotype_region, haplotype_to_reference);
    assert(!padded_haplotype_cigar.empty());
    if (is_deletion(padded_haplotype_cigar.front())) {
        padded_haplotype_cigar.erase(std::cbegin(padded_haplotype_cigar));
    }
    auto rebased_cigar = rebase(read.cigar(), padded_haplotype_cigar);
    auto rebased_read_region = expand_rhs(shift(head_region(read), rebase_shift), reference_size(rebased_cigar));
    read.realign(std::move(rebased_read_region), std::move(rebased_cigar));
}

void rebase_not_overlapped(AlignedRead& read, const GenomicRegion& haplotype_region, const CigarString& haplotype_to_reference,
                           GenomicRegion rebased_read_region)
{
    if (overlaps(rebased_read_region, haplotype_region)) {
        auto padded_haplotype_cigar = copy_tail(haplotype_region, haplotype_to_reference, rebased_read_region);
        assert(!padded_haplotype_cigar.empty());
        if (is_insertion(padded_haplotype_cigar.front())) {
            const auto called_insertion_size = padded_haplotype_cigar.front().size();
            const auto rebase_shift = static_cast<GenomicRegion::Size>(begin_distance(rebased_read_region, read));
            const auto supported_insertion_size = called_insertion_size - std::min(rebase_shift, called_insertion_size);
            assert(supported_insertion_size <= called_insertion_size);
            if (supported_insertion_size > 0) {
                padded_haplotype_cigar.front().set_size(supported_insertion_size);
                if (is_match_or_substitution(padded_haplotype_cigar.back())) {
                    increment_size(padded_haplotype_cigar.back(), supported_insertion_size);
                } else {
                    padded_haplotype_cigar.emplace_back(supported_insertion_size, CigarOperation::Flag::sequenceMatch);
                }
            } else {
                padded_haplotype_cigar.erase(std::cbegin(padded_haplotype_cigar));
            }
        }
        if (reference_size(read.cigar()) > sequence_size(padded_haplotype_cigar)) {
            const auto pad = reference_size(read.cigar()) - sequence_size(padded_haplotype_cigar);
            if (is_match_or_substitution(padded_haplotype_cigar.back())) {
                increment_size(padded_haplotype_cigar.back(), pad);
            } else {
                padded_haplotype_cigar.emplace_back(pad, CigarOperation::Flag::sequenceMatch);
            }
        }
        auto rebased_cigar = rebase(read.cigar(), padded_haplotype_cigar);
        rebased_read_region = expand_rhs(head_region(rebased_read_region), reference_size(rebased_cigar));
        read.realign(std::move(rebased_read_region), std::move(rebased_cigar));
    } else if (are_adjacent(haplotype_region, rebased_read_region) && is_insertion(haplotype_to_reference.back())) {
        const auto insertion_size = get_tail_op_size(haplotype_region, haplotype_to_reference, mapped_region(read));
        if (insertion_size > 0) {
            using Flag = CigarOperation::Flag;
            const CigarString padded_haplotype_cigar {CigarOperation {insertion_size, Flag::insertion},
                                                      CigarOperation {size(rebased_read_region), Flag::sequenceMatch}};
            auto rebased_cigar = rebase(read.cigar(), padded_haplotype_cigar);
            rebased_read_region = expand_rhs(head_region(rebased_read_region), reference_size(rebased_cigar));
            read.realign(std::move(rebased_read_region), std::move(rebased_cigar));
        } else {
            read.realign(std::move(rebased_read_region), read.cigar());
        }
    } else {
        read.realign(std::move(rebased_read_region), read.cigar());
    }
}

void rebase_adjacent(AlignedRead& read, const GenomicRegion& haplotype_region, const CigarString& haplotype_to_reference)
{
    assert(!haplotype_to_reference.empty());
    if (is_before(haplotype_region, read) && is_insertion(haplotype_to_reference.back())) {
        using Flag = CigarOperation::Flag;
        const CigarString padded_haplotype_cigar {haplotype_to_reference.back(),
                                                  CigarOperation {region_size(read), Flag::sequenceMatch}};
        auto rebased_read_cigar = rebase(read.cigar(), padded_haplotype_cigar);
        auto rebased_read_region = expand_rhs(head_region(read), reference_size(rebased_read_cigar));
        read.realign(std::move(rebased_read_region), std::move(rebased_read_cigar));
    }
}

void rebase(AlignedRead& read, const GenomicRegion& haplotype_region, const CigarString& haplotype_to_reference)
{
    const auto rebase_shift = calculate_rebase_shift(read, haplotype_region, haplotype_to_reference);
    if (overlaps(read, haplotype_region)) {
        rebase_overlapped(read, haplotype_region, haplotype_to_reference, rebase_shift);
    } else if (rebase_shift != 0) {
        auto rebased_read_region = shift(mapped_region(read), rebase_shift);
        if (rebase_shift < 0) {
            rebase_not_overlapped(read, haplotype_region, haplotype_to_reference, std::move(rebased_read_region));
        } else {
            read.realign(std::move(rebased_read_region), read.cigar());
        }
    } else if (are_adjacent(haplotype_region, read)) {
        rebase_adjacent(read, haplotype_region, haplotype_to_reference);
    }
}

} // namespace

void rebase(std::vector<AlignedRead>& reads, const Haplotype& haplotype, boost::optional<ThreadPool&> workers)
{
    const auto haplotype_cigar = haplotype.cigar();
    using octopus::for_each;
    for_each(std::begin(reads), std::end(reads), [&] (auto& read) {
        rebase(read, haplotype.mapped_region(), haplotype_cigar);
    });
}

void realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype,
                          HaplotypeLikelihoodModel model)
{
    realign(reads, haplotype, std::move(model));
    rebase(reads, haplotype);
}

void realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    realign_to_reference(reads, haplotype, make_default_haplotype_likelihood_model());
}

std::vector<AlignedRead>
realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model)
{
    auto result = realign(reads, haplotype, std::move(model));
    rebase(result, haplotype);
    return result;
}

std::vector<AlignedRead> realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    return realign_to_reference(reads, haplotype, make_default_haplotype_likelihood_model());
}

void safe_realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model)
{
    safe_realign(reads, haplotype, std::move(model));
    rebase(reads, haplotype);
}

void safe_realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    safe_realign_to_reference(reads, haplotype, make_default_haplotype_likelihood_model());
}

std::vector<AlignedRead>
safe_realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model)
{
    auto result = safe_realign(reads, haplotype, std::move(model));
    rebase(result, haplotype);
    return result;
}

std::vector<AlignedRead> safe_realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    return safe_realign_to_reference(reads, haplotype, make_default_haplotype_likelihood_model());
}

} // namespace octopus
