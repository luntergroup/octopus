// Copyright (c) 2017 Daniel Cooke
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
#include "core/models/error/error_model_factory.hpp"

namespace octopus {

namespace {

auto deletion_length(const AlignedRead& read) noexcept
{
    return region_size(read) > sequence_size(read) ? sequence_size(read) - region_size(read) : 0;
}

GenomicRegion::Size max_deletion_length(const std::vector<AlignedRead>& reads) noexcept
{
    if (reads.empty()) return 0;
    auto max_itr = std::max_element(std::cbegin(reads), std::cend(reads),
                                    [] (const auto& lhs, const auto& rhs) {
                                        return deletion_length(lhs) < deletion_length(rhs);
                                    });
    return deletion_length(*max_itr);
}

auto safe_expand(const GenomicRegion& region, const GenomicRegion::Size n)
{
    if (region.begin() < n) {
        return GenomicRegion {region.contig_name(), 0, region.end() + 2 * n - region.begin()};
    } else {
        return expand(region, n);
    }
}

} // namespace

Haplotype expand_for_realignment(const Haplotype& haplotype, const std::vector<AlignedRead>& reads)
{
    const auto reads_region = encompassing_region(reads);
    auto required_flank_pad = 2 * HaplotypeLikelihoodModel::pad_requirement();
    if (region_size(haplotype) > sequence_size(haplotype)) {
        required_flank_pad += region_size(haplotype) - sequence_size(haplotype);
    }
    const auto required_haplotype_region = safe_expand(reads_region, required_flank_pad);
    if (contains(haplotype, required_haplotype_region)) {
        return haplotype;
    } else {
        return remap(haplotype, encompassing_region(haplotype, required_haplotype_region));
    }
}

namespace {

auto compute_read_hashes(const std::vector<AlignedRead>& reads)
{
    static constexpr unsigned char mapperKmerSize {6};
    std::vector<KmerPerfectHashes> result {};
    result.reserve(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(result),
                   [=] (const AlignedRead& read) { return compute_kmer_hashes<mapperKmerSize>(read.sequence()); });
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

} // namespace

void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    realign(reads, haplotype, HaplotypeLikelihoodModel {nullptr, make_indel_error_model(), false});
}

void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model)
{
    const auto read_hashes = compute_read_hashes(reads);
    static constexpr unsigned char mapperKmerSize {6};
    auto haplotype_hashes = init_kmer_hash_table<mapperKmerSize>();
    populate_kmer_hash_table<mapperKmerSize>(haplotype.sequence(), haplotype_hashes);
    auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
    model.reset(haplotype);
    for (std::size_t i {0}; i < reads.size(); ++i) {
        auto mapping_positions = map_query_to_target(read_hashes[i], haplotype_hashes, haplotype_mapping_counts);
        reset_mapping_counts(haplotype_mapping_counts);
        realign(reads[i], haplotype, model.align(reads[i], mapping_positions));
    }
}

std::vector<AlignedRead> realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype,
                                 HaplotypeLikelihoodModel model)
{
    auto result = reads;
    realign(result, haplotype, model);
    return result;
}

std::vector<AlignedRead> realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    return realign(reads, haplotype, HaplotypeLikelihoodModel {nullptr, make_indel_error_model(), false});
}

void safe_realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    auto expanded_haplotype = expand_for_realignment(haplotype, reads);
    try {
        realign(reads, expanded_haplotype);
    } catch (const HaplotypeLikelihoodModel::ShortHaplotypeError& e) {
        expanded_haplotype = expand(expanded_haplotype, e.required_extension());
        realign(reads, expanded_haplotype);
    }
}

std::vector<AlignedRead> safe_realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    auto result = reads;
    safe_realign(result, haplotype);
    return result;
}

namespace {

bool is_insertion(CigarOperation::Flag flag) noexcept
{
    return flag == CigarOperation::Flag::insertion;
}

bool is_insertion(const CigarOperation& op) noexcept
{
    return is_insertion(op.flag());
}

bool is_deletion(CigarOperation::Flag flag) noexcept
{
    return flag == CigarOperation::Flag::deletion;
}

bool is_deletion(const CigarOperation& op) noexcept
{
    return is_deletion(op.flag());
}

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
    assert(is_match(haplotype) && is_match(read));
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
        if (is_match(read_op)) {
            for (unsigned n {0}; n < read_op.size();) {
                assert(hap_flag_itr != std::cend(haplotypes_ops));
                if (is_match(*hap_flag_itr)) {
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
        result = copy_sequence(haplotype_to_reference, left_overhang_size(haplotype_region, read_region), sequence_size(haplotype_to_reference));
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
                        result.emplace_back(rhs_pad_size, Flag::sequenceMatch);
                    }
                }
            } else {
                if (lhs_pad_size > 0) result.emplace_back(lhs_pad_size, Flag::sequenceMatch);
                if (is_sequence_match(haplotype_to_reference.back())) {
                    result.insert(result.cend(), haplotype_to_reference.cbegin(), haplotype_to_reference.cend() - 1);
                    result.emplace_back(haplotype_to_reference.back().size() + rhs_pad_size, Flag::sequenceMatch);
                } else {
                    utils::append(haplotype_to_reference, result);
                    if (rhs_pad_size) result.emplace_back(rhs_pad_size, Flag::sequenceMatch);
                }
            }
        } else if (begins_before(read_region, haplotype_region)) {
            assert(ends_before(read_region, haplotype_region));
            const auto lhs_pad_size = left_overhang_size(read_region, haplotype_region);
            if (is_sequence_match(haplotype_to_reference.front())) {
                result.emplace_back(lhs_pad_size + haplotype_to_reference.front().size(), Flag::sequenceMatch);
                result.insert(result.cend(), haplotype_to_reference.cbegin() + 1, haplotype_to_reference.cend());
            } else {
                result.emplace_back(lhs_pad_size, Flag::sequenceMatch);
                utils::append(haplotype_to_reference, result);
            }
        } else {
            assert(begins_before(haplotype_region, read_region) && ends_before(haplotype_region, read_region));
            const auto offset = left_overhang_size(haplotype_region, read_region);
            const auto rhs_pad_size = right_overhang_size(read_region, haplotype_region);
            result = copy_sequence(haplotype_to_reference, offset, sequence_size(haplotype_to_reference));
            assert(!result.empty());
            if (is_sequence_match(result.back())) {
                result.back() = CigarOperation {result.back().size() + rhs_pad_size, Flag::sequenceMatch};
            } else if (rhs_pad_size > 0) {
                result.emplace_back(rhs_pad_size, Flag::sequenceMatch);
            }
        }
    }
    if (sequence_size(result) < reference_size(read_to_haplotype)) {
        assert(!result.empty());
        const auto rhs_pad_size = reference_size(read_to_haplotype) - sequence_size(result);
        if (is_sequence_match(result.back())) {
            result.back() = CigarOperation {result.back().size() + rhs_pad_size, Flag::sequenceMatch};
        } else {
            result.emplace_back(rhs_pad_size, Flag::sequenceMatch);
        }
    }
    return result;
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

void rebase(AlignedRead& read, const GenomicRegion& haplotype_region, const CigarString& haplotype_to_reference)
{
    const auto rebase_sift = calculate_rebase_shift(read, haplotype_region, haplotype_to_reference);
    if (overlaps(read, haplotype_region)) {
        const auto padded_haplotype_cigar = pad_reference(read.mapped_region(), read.cigar(), haplotype_region, haplotype_to_reference);
        auto rebased_cigar = rebase(read.cigar(), padded_haplotype_cigar);
        auto rebased_read_region = expand_rhs(shift(head_region(read), rebase_sift), reference_size(rebased_cigar));
        read.realign(std::move(rebased_read_region), std::move(rebased_cigar));
    } else if (rebase_sift != 0) {
        read.realign(shift(mapped_region(read), rebase_sift), read.cigar());
    }
}

} // namespace

void rebase(std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    const auto haplotype_cigar = haplotype.cigar();
    for (auto& read : reads) {
        rebase(read, haplotype.mapped_region(), haplotype_cigar);
    }
    std::sort(std::begin(reads), std::end(reads));
}

void realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    realign(reads, haplotype);
    rebase(reads, haplotype);
}

std::vector<AlignedRead> realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    auto result = realign(reads, haplotype);
    rebase(result, haplotype);
    return result;
}

void safe_realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    safe_realign(reads, haplotype);
    rebase(reads, haplotype);
}

std::vector<AlignedRead> safe_realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    auto result = safe_realign(reads, haplotype);
    rebase(result, haplotype);
    return result;
}

} // namespace octopus
