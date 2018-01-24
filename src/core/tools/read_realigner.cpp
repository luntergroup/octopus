// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_realigner.hpp"

#include <iterator>
#include <algorithm>
#include <utility>

#include "basics/genomic_region.hpp"
#include "utils/maths.hpp"
#include "utils/kmer_mapper.hpp"
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

AlignedRead realign(const AlignedRead& read, GenomicRegion region, CigarString cigar)
{
    return AlignedRead {read.name(), std::move(region), read.sequence(), read.base_qualities(),
                        std::move(cigar), read.mapping_quality(), read.flags()};
}

AlignedRead realign(const AlignedRead& read, const Haplotype& haplotype,
                    GenomicRegion::Position mapping_position, CigarString alignment)
{
    const auto remapped_read_begin = mapped_begin(haplotype) + mapping_position;
    const auto remapped_read_end = remapped_read_begin + reference_size(alignment);
    return realign(read, GenomicRegion {contig_name(read), remapped_read_begin, remapped_read_end}, std::move(alignment));
}

AlignedRead realign(const AlignedRead& read, const Haplotype& haplotype, HaplotypeLikelihoodModel::Alignment alignment)
{
    return realign(read, haplotype, alignment.mapping_position, std::move(alignment.cigar));
}

} // namespace

std::vector<AlignedRead> realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype,
                                 HaplotypeLikelihoodModel model)
{
    const auto read_hashes = compute_read_hashes(reads);
    static constexpr unsigned char mapperKmerSize {6};
    auto haplotype_hashes = init_kmer_hash_table<mapperKmerSize>();
    populate_kmer_hash_table<mapperKmerSize>(haplotype.sequence(), haplotype_hashes);
    auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
    model.reset(haplotype);
    std::vector<AlignedRead> result {};
    result.reserve(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::cbegin(read_hashes), std::back_inserter(result),
                   [&] (const auto& read, const auto& read_hash) {
                       auto mapping_positions = map_query_to_target(read_hash, haplotype_hashes, haplotype_mapping_counts);
                       reset_mapping_counts(haplotype_mapping_counts);
                       return realign(read, haplotype, model.align(read, mapping_positions));
                   });
    return result;
}

std::vector<AlignedRead> realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    return realign(reads, haplotype, HaplotypeLikelihoodModel {nullptr, make_indel_error_model(), false});
}

std::vector<AlignedRead> safe_realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype)
{
    auto expanded_haplotype = expand_for_realignment(haplotype, reads);
    try {
        return realign(reads, expanded_haplotype);
    } catch (const HaplotypeLikelihoodModel::ShortHaplotypeError& e) {
        expanded_haplotype = expand(expanded_haplotype, e.required_extension());
        return realign(reads, expanded_haplotype);
    }
}

} // namespace octopus
