// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "local_reassembler.hpp"

#include <algorithm>
#include <iterator>
#include <deque>
#include <stdexcept>
#include <thread>
#include <future>
#include <cassert>

#include "tandem/tandem.hpp"

#include "concepts/mappable.hpp"
#include "concepts/mappable_range.hpp"
#include "basics/contig_region.hpp"
#include "basics/cigar_string.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/sequence_utils.hpp"
#include "utils/append.hpp"
#include "utils/global_aligner.hpp"
#include "utils/read_stats.hpp"
#include "utils/free_memory.hpp"
#include "io/reference/reference_genome.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace coretools {

namespace {

void remove_duplicates(std::vector<unsigned>& kmer_sizes)
{
    std::sort(std::begin(kmer_sizes), std::end(kmer_sizes));
    kmer_sizes.erase(std::unique(std::begin(kmer_sizes), std::end(kmer_sizes)), std::end(kmer_sizes));
}

auto generate_fallback_kmer_sizes(std::vector<unsigned>& result,
                                  const std::vector<unsigned>& default_kmer_sizes,
                                  const unsigned num_fallbacks, const unsigned interval_size)
{
    assert(!default_kmer_sizes.empty());
    result.resize(num_fallbacks);
    auto k = default_kmer_sizes.back();
    std::generate_n(std::begin(result), num_fallbacks,
                    [&k, interval_size] () noexcept -> decltype(k) {
                        return k += interval_size;
                    });
}

} // namespace

LocalReassembler::LocalReassembler(const ReferenceGenome& reference, Options options)
: execution_policy_ {ExecutionPolicy::par}
, reference_ {reference}
, default_kmer_sizes_ {std::move(options.kmer_sizes)}
, fallback_kmer_sizes_ {}
, read_buffer_ {}
, max_bin_size_ {options.bin_size}
, max_bin_overlap_ {options.bin_overlap}
, mask_threshold_ {options.mask_threshold}
, min_kmer_observations_ {options.min_kmer_observations}
, max_bubbles_ {options.max_bubbles}
, min_bubble_score_ {options.min_bubble_score}
, max_variant_size_ {options.max_variant_size}
, cycle_tolerance_ {options.cycle_tolerance}
, ignore_strand_bias_ {options.ignore_strand_bias}
{
    if (max_bin_size_ == 0) {
        throw std::runtime_error {"bin size must be greater than zero"};
    }
    if (max_bin_overlap_ >= max_bin_size_) {
        max_bin_overlap_ = max_bin_size_ - 1;
    }
    if (options.fallback_interval_size == 0) {
        throw std::runtime_error {"fallback interval size must be greater than zero"};
    }
    if (default_kmer_sizes_.empty()) return;
    remove_duplicates(default_kmer_sizes_);
    generate_fallback_kmer_sizes(fallback_kmer_sizes_, default_kmer_sizes_,
                                 options.num_fallbacks, options.fallback_interval_size);
}

std::unique_ptr<VariantGenerator> LocalReassembler::do_clone() const
{
    return std::make_unique<LocalReassembler>(*this);
}

LocalReassembler::Bin::Bin(GenomicRegion region)
: region {std::move(region)}
{}

const GenomicRegion& LocalReassembler::Bin::mapped_region() const noexcept
{
    return region;
}

void LocalReassembler::Bin::add(const AlignedRead& read, const std::size_t sample_index)
{
    if (read_region) {
        read_region = encompassing_region(*read_region, contig_region(read));
    } else {
        read_region = contig_region(read);
    }
    if (is_forward_strand(read)) {
        forward_read_sequences.push_back({read.sequence(), read.base_qualities(), sample_index});
    } else {
        reverse_read_sequences.push_back({read.sequence(), read.base_qualities(), sample_index});
    }
}

void LocalReassembler::Bin::add(const AlignedRead& read, const NucleotideSequence& masked_sequence, const std::size_t sample_index)
{
    if (read_region) {
        read_region = encompassing_region(*read_region, contig_region(read));
    } else {
        read_region = contig_region(read);
    }
    if (is_forward_strand(read)) {
        forward_read_sequences.push_back({masked_sequence, read.base_qualities(), sample_index});
    } else {
        reverse_read_sequences.push_back({masked_sequence, read.base_qualities(), sample_index});
    }
}

void LocalReassembler::Bin::clear() noexcept
{
    forward_read_sequences.clear();
    forward_read_sequences.shrink_to_fit();
    reverse_read_sequences.clear();
    reverse_read_sequences.shrink_to_fit();
}

std::size_t LocalReassembler::Bin::size() const noexcept
{
    return forward_read_sequences.size() + reverse_read_sequences.size();
}

bool LocalReassembler::Bin::empty() const noexcept
{
    return forward_read_sequences.empty() && reverse_read_sequences.empty();
}

bool LocalReassembler::do_requires_reads() const noexcept
{
    return true;
}

namespace {

bool has_low_quality_flank(const AlignedRead& read, const AlignedRead::BaseQuality good_quality) noexcept
{
    if (is_soft_clipped(read)) {
        if (is_front_soft_clipped(read) && read.base_qualities().front() < good_quality) {
            return true;
        } else {
            return is_back_soft_clipped(read) && read.base_qualities().back() < good_quality;
        }
    } else {
        return false;
    }
}

bool has_low_quality_match(const AlignedRead& read, const AlignedRead::BaseQuality good_quality) noexcept
{
    if (good_quality == 0) return false;
    auto quality_itr = std::cbegin(read.base_qualities());
    return std::any_of(std::cbegin(read.cigar()), std::cend(read.cigar()),
                       [&] (const auto& op) {
                           if (is_match_or_substitution(op)) {
                               auto result = std::any_of(quality_itr, std::next(quality_itr, op.size()),
                                                         [=] (auto q) { return q < good_quality; });
                               std::advance(quality_itr, op.size());
                               return result;
                           } else if (advances_sequence(op)) {
                               std::advance(quality_itr, op.size());
                           }
                           return false;
                       });
}

bool requires_masking(const AlignedRead& read, const AlignedRead::BaseQuality good_quality) noexcept
{
    return has_low_quality_flank(read, good_quality) || has_low_quality_match(read, good_quality);
}

auto find_first_sequence_op(const std::vector<CigarOperation::Flag>& cigar) noexcept
{
    return std::find_if_not(std::cbegin(cigar), std::cend(cigar),
                            [] (auto op) { return op == CigarOperation::Flag::hardClipped; });
}

template <typename T>
auto make_optional(bool b, T&& value)
{
    if (b) {
        return boost::optional<T> {std::forward<T>(value)};
    } else {
        return boost::optional<T> {};
    }
}

auto transform_low_quality_matches_to_reference(AlignedRead::NucleotideSequence read_sequence,
                                                const AlignedRead::BaseQualityVector& base_qualities,
                                                const AlignedRead::NucleotideSequence& reference_sequence,
                                                const CigarString& cigar,
                                                const AlignedRead::BaseQuality min_quality)
{
    const auto expanded_cigar = decompose(cigar);
    auto ref_itr   = std::cbegin(reference_sequence);
    auto cigar_itr = find_first_sequence_op(expanded_cigar);
    bool has_masked {false};
    std::transform(std::cbegin(read_sequence), std::cend(read_sequence), std::cbegin(base_qualities),
                   std::begin(read_sequence), [&] (const auto read_base, const auto base_quality) {
        using Flag = CigarOperation::Flag;
        // Deletions are excess reference sequence so we need to move the
        // reference iterator to the next non-deleted read base
        while (cigar_itr != std::cend(expanded_cigar) && *cigar_itr == Flag::deletion) {
            ++cigar_itr;
            ++ref_itr;
        }
        const auto op = *cigar_itr++;
        if (is_match_or_substitution(op)) {
            const auto ref_base = *ref_itr++; // Don't forget to increment ref_itr!
            if (base_quality >= min_quality) {
                return read_base;
            } else {
                has_masked = true;
                return ref_base;
            }
        } else {
            if (op != Flag::insertion) ++ref_itr;
            return read_base;
        }
    });
    return make_optional(has_masked, std::move(read_sequence));
}

auto transform_low_quality_matches_to_reference(const AlignedRead& read,
                                                const AlignedRead::BaseQuality min_quality,
                                                const ReferenceGenome& reference)
{
    return transform_low_quality_matches_to_reference(read.sequence(), read.base_qualities(),
                                                      reference.fetch_sequence(mapped_region(read)),
                                                      read.cigar(), min_quality);
}

auto get_removable_flank_sizes(const AlignedRead& read, const AlignedRead::BaseQuality min_quality) noexcept
{
    CigarOperation::Size front_clip, back_clip;
    std::tie(front_clip, back_clip) = get_soft_clipped_sizes(read);
    AlignedRead::NucleotideSequence::size_type front {0}, back {0};
    const auto is_low_quality = [min_quality] (auto q) { return q < min_quality; };
    const auto& base_qualities = read.base_qualities();
    if (front_clip > 0) {
        const auto begin = std::cbegin(base_qualities);
        const auto first_good = std::find_if_not(begin, std::next(begin, front_clip), is_low_quality);
        front = std::distance(begin, first_good);
    }
    if (back_clip > 0) {
        const auto begin = std::crbegin(base_qualities);
        const auto first_good = std::find_if_not(begin, std::next(begin, back_clip), is_low_quality);
        back = std::distance(begin, first_good);
    }
    return std::make_pair(front, back);
}

auto mask(const AlignedRead& read, const AlignedRead::BaseQuality min_quality, const ReferenceGenome& reference)
{
    auto result = transform_low_quality_matches_to_reference(read, min_quality, reference);
    if (result && has_low_quality_flank(read, min_quality)) {
        const auto p = get_removable_flank_sizes(read, min_quality);
        assert(p.first + p.second < sequence_size(read));
        result->erase(std::prev(std::cend(*result), p.second), std::cend(*result));
        result->erase(std::cbegin(*result), std::next(std::cbegin(*result), p.first));
    }
    return result;
}

} // namespace

template <typename Container, typename M>
auto overlapped_bins(Container& bins, const M& mappable)
{
    return bases(overlap_range(std::begin(bins), std::end(bins), mappable, BidirectionallySortedTag {}));
}

void LocalReassembler::do_add_read(const SampleName& sample, const AlignedRead& read)
{
    read_buffer_[sample].insert(read);
}

void LocalReassembler::do_add_reads(const SampleName& sample, ReadVectorIterator first, ReadVectorIterator last)
{
    read_buffer_[sample].insert(first, last);
}

void LocalReassembler::do_add_reads(const SampleName& sample, ReadFlatSetIterator first, ReadFlatSetIterator last)
{
    read_buffer_[sample].insert(first, last);
}

void remove_nonoverlapping(std::vector<Variant>& candidates, std::vector<GenomicRegion>& active_regions)
{
    const auto it = std::remove_if(std::begin(candidates), std::end(candidates),
                                   [&] (const Variant& candidate) {
                                       return !has_overlapped(active_regions, candidate);
                                   });
    candidates.erase(it, std::end(candidates));
}

void remove_duplicates(std::deque<Variant>& variants)
{
    std::sort(std::begin(variants), std::end(variants));
    variants.erase(std::unique(std::begin(variants), std::end(variants)), std::end(variants));
}

void remove_larger_than(std::deque<Variant>& variants, const Variant::MappingDomain::Size max_size)
{
    variants.erase(std::remove_if(std::begin(variants), std::end(variants),
                                  [max_size] (const auto& variant) { return region_size(variant) > max_size; }),
                   std::end(variants));
}

std::vector<Variant> LocalReassembler::do_generate(const RegionSet& regions, OptionalThreadPool workers) const
{
    BinList bins {};
    SequenceBuffer masked_sequence_buffer {};
    for (auto& p : read_buffer_) {
        std::sort(std::begin(p.second), std::end(p.second));
    }
    for (const auto& region : regions) {
        prepare_bins(region, bins);
        std::size_t sample_idx {0};
        for (const auto& p : read_buffer_) {
            for (const auto& read : overlap_range(p.second, region)) {
                auto active_bins = overlapped_bins(bins, read);
                assert(!active_bins.empty());
                if (requires_masking(read, mask_threshold_)) {
                    auto masked_sequence = mask(read, mask_threshold_, reference_);
                    if (masked_sequence) {
                        masked_sequence_buffer.emplace_back(std::move(*masked_sequence));
                        for (auto& bin : active_bins) {
                            bin.add(read, std::cref(masked_sequence_buffer.back()), sample_idx);
                        }
                    }
                } else {
                    for (auto& bin : active_bins) bin.add(read, sample_idx);
                }
            }
            ++sample_idx;
        }
    }
    finalise_bins(bins, regions);
    if (bins.empty()) return {};
    std::deque<Variant> candidates {};
    if (workers && bins.size() > 1 && workers->size() > 1) {
        std::vector<std::future<std::deque<Variant>>> bin_futures {};
        bin_futures.reserve(bins.size());
        for (auto& bin : bins) {
            bin_futures.push_back(workers->try_push([&] () {
                if (debug_log_) {
                    stream(*debug_log_) << "Assembling " << bin.size() << " reads in bin " << mapped_region(bin);
                }
                std::deque<Variant> result {};
                    const auto num_default_failures = try_assemble_with_defaults(bin, result);
                    if (num_default_failures == default_kmer_sizes_.size()) {
                        try_assemble_with_fallbacks(bin, result);
                    }
                    bin.clear();
                    return result;
            }));
        }
        for (auto&& f : bin_futures) {
            utils::append(f.get(), candidates);
        }
    } else {
        for (auto& bin : bins) {
            if (debug_log_) {
                stream(*debug_log_) << "Assembling " << bin.size() << " reads in bin " << mapped_region(bin);
            }
            const auto num_default_failures = try_assemble_with_defaults(bin, candidates);
            if (num_default_failures == default_kmer_sizes_.size()) {
                try_assemble_with_fallbacks(bin, candidates);
            }
            bin.clear();
        }
    }
    remove_duplicates(candidates);
    remove_larger_than(candidates, max_variant_size_);
    return {std::make_move_iterator(std::begin(candidates)), std::make_move_iterator(std::end(candidates))};
}

void LocalReassembler::do_clear() noexcept
{
    free_memory(read_buffer_);
}

std::string LocalReassembler::name() const
{
    return "LocalReassembler";
}

// private methods

template <typename MappableTp>
auto decompose(const MappableTp& mappable, const GenomicRegion::Position n,
               const GenomicRegion::Size overlap = 0)
{
    if (overlap >= n) {
        throw std::runtime_error {"decompose: overlap must be less than n"};
    }
    std::vector<GenomicRegion> result {};
    if (n == 0) return result;
    const auto num_elements = region_size(mappable) / (n - overlap);
    if (num_elements == 0) return result;
    result.reserve(num_elements);
    const auto& contig = contig_name(mappable);
    auto curr = mapped_begin(mappable);
    std::generate_n(std::back_inserter(result), num_elements, [&contig, &curr, n, overlap] () {
        auto tmp = curr;
        curr += (n - overlap);
        return GenomicRegion {contig, tmp, tmp + n};
    });
    return result;
}

void LocalReassembler::prepare_bins(const GenomicRegion& region, BinList& bins) const
{
    assert(bins.empty() || is_after(region, bins.back()));
    if (size(region) > max_bin_size_) {
        auto bin_region = expand_rhs(head_region(region), max_bin_size_);
        while (ends_before(bin_region, region)) {
            bins.push_back(bin_region);
            bin_region = shift(bin_region, max_bin_overlap_);
        }
        if (overlap_size(region, bin_region) > 0) {
            bins.push_back(*overlapped_region(region, bin_region));
        }
    } else {
        bins.push_back(region);
    }
}

bool LocalReassembler::should_assemble_bin(const Bin& bin) const
{
    return !bin.empty();
}

void LocalReassembler::finalise_bins(BinList& bins, const RegionSet& active_regions) const
{
    bins.erase(std::remove_if(std::begin(bins), std::end(bins),
                              [this] (const Bin& bin) { return !should_assemble_bin(bin); }),
               std::end(bins));
    for (auto& bin : bins) {
        if (bin.read_region) {
            if (size(bin.region) < max_bin_size_) {
                bin.region = expand(bin.region, (max_bin_size_ - size(bin.region)) / 2);
            }
            assert(overlaps(bin.region.contig_region(), *bin.read_region));
            bin.region = GenomicRegion {bin.region.contig_name(), *overlapped_region(bin.region.contig_region(), *bin.read_region)};
        }
    }
    // unique in reverse order as we want to keep bigger bins, which
    // are sorted after smaller bins with the same starting point
    bins.erase(std::begin(bins), std::unique(std::rbegin(bins), std::rend(bins),
                                             [] (const Bin& lhs, const Bin& rhs) {
                                                 return begins_equal(lhs, rhs);
                                             }).base());
}

namespace {

template <typename L>
void log_success(L& log, const char* type, const unsigned k)
{
    if (log) stream(*log, 8) << type << " assembler with kmer size " << k << " completed";
}

template <typename L>
void log_partial_success(L& log, const char* type, const unsigned k)
{
    if (log) stream(*log, 8) << type << " assembler with kmer size " << k << " partially completed";
}

template <typename L>
void log_failure(L& log, const char* type, const unsigned k)
{
    if (log) stream(*log, 8) << type << " assembler with kmer size " << k << " failed";
}

} // namespace

unsigned LocalReassembler::try_assemble_with_defaults(const Bin& bin, std::deque<Variant>& result) const
{
    unsigned num_failures {0};
    for (const auto k : default_kmer_sizes_) {
        const auto status = assemble_bin(k, bin, result);
        switch (status) {
            case AssemblerStatus::success:
                log_success(debug_log_, "Default", k);
                break;
            case AssemblerStatus::partial_success:
                log_partial_success(debug_log_, "Default", k);
                ++num_failures;
                break;
            default:
                log_failure(debug_log_, "Default", k);
                ++num_failures;
        }
    }
    return num_failures;
}

void LocalReassembler::try_assemble_with_fallbacks(const Bin& bin, std::deque<Variant>& result) const
{
    auto prev_k = default_kmer_sizes_.back();
    for (const auto k : fallback_kmer_sizes_) {
        const auto status = assemble_bin(k, bin, result);
        switch (status) {
            case AssemblerStatus::success:
                log_success(debug_log_, "Fallback", k);
                if (k - prev_k > 5) {
                    const auto gap = k - prev_k;
                    assemble_bin(k - gap / 2, bin, result);
                    assemble_bin(k + gap / 2, bin, result);
                }
                return;
            case AssemblerStatus::partial_success:
                log_partial_success(debug_log_, "Fallback", k);
                break;
            default:
                log_failure(debug_log_, "Fallback", k);
        }
        prev_k = k;
    }
}

GenomicRegion LocalReassembler::propose_assembler_region(const GenomicRegion& input_region, unsigned kmer_size) const
{
    if (input_region.begin() < kmer_size) {
        const auto& contig = input_region.contig_name();
        if (reference_.get().contig_size(contig) >= kmer_size) {
            return GenomicRegion {contig, 0, input_region.end() + kmer_size};
        } else {
            return reference_.get().contig_region(contig);
        }
    } else {
        auto ideal_proposal = expand(input_region, kmer_size);
        if (reference_.get().contains(ideal_proposal)) {
            return ideal_proposal;
        } else {
            const auto& contig = input_region.contig_name();
            const auto end = reference_.get().contig_size(contig);
            return GenomicRegion {contig, input_region.begin() - kmer_size, end};
        }
    }
}

void LocalReassembler::load(const Bin& bin, Assembler& assembler) const
{
    for (const auto& read : bin.forward_read_sequences) {
        assembler.insert_read(read.sequence, read.base_qualities, Assembler::Direction::forward, read.sample_index);
    }
    for (const auto& read : bin.reverse_read_sequences) {
        assembler.insert_read(read.sequence, read.base_qualities, Assembler::Direction::reverse, read.sample_index);
    }
}

void LocalReassembler::load(const Bin& bin, const std::size_t sample_idx, Assembler& assembler) const
{
    for (const auto& read : bin.forward_read_sequences) {
        if (read.sample_index == sample_idx) {
            assembler.insert_read(read.sequence, read.base_qualities, Assembler::Direction::forward, read.sample_index);
        }
    }
    for (const auto& read : bin.reverse_read_sequences) {
        if (read.sample_index == sample_idx) {
            assembler.insert_read(read.sequence, read.base_qualities, Assembler::Direction::reverse, read.sample_index);
        }
    }
}

LocalReassembler::AssemblerStatus
LocalReassembler::assemble_bin(const unsigned kmer_size, const Bin& bin, std::deque<Variant>& result) const
{
    if (bin.empty()) return AssemblerStatus::success;
    const auto assemble_region = propose_assembler_region(bin.region, kmer_size);
    if (size(assemble_region) < kmer_size) return AssemblerStatus::failed;
    const auto reference_sequence = reference_.get().fetch_sequence(assemble_region);
    if (!utils::is_canonical_dna(reference_sequence)) return AssemblerStatus::failed;
    Assembler::Parameters assembler_params {kmer_size};
    assembler_params.use_strand_bias = !ignore_strand_bias_;
    Assembler assembler {assembler_params, reference_sequence};
    if (assembler.is_unique_reference()) {
        load(bin, assembler);
        auto status = try_assemble_region(assembler, reference_sequence, assemble_region, result);
        const auto num_samples = read_buffer_.size();
        if (num_samples > 1 && status != AssemblerStatus::success) {
            for (std::size_t sample_idx {0}; sample_idx < num_samples; ++sample_idx) {
                Assembler sample_assembler {assembler_params, reference_sequence};
                load(bin, sample_idx, sample_assembler);
                try_assemble_region(sample_assembler, reference_sequence, assemble_region, result);
            }
        }
        return status;
    } else {
        return AssemblerStatus::failed;
    }
}

bool is_inversion(const Assembler::Variant& v) noexcept
{
    return v.ref.size() > 2
           && utils::are_reverse_complements(v.ref, v.alt)
           && !utils::is_homopolymer(v.ref)
           && !std::equal(std::next(std::cbegin(v.ref)), std::prev(std::cend(v.ref)), std::next(std::cbegin(v.alt)));
}

void trim_reference(Assembler::Variant& v)
{
    using std::cbegin; using std::cend; using std::crbegin; using std::crend;
    const auto p1 = std::mismatch(crbegin(v.ref), crend(v.ref), crbegin(v.alt), crend(v.alt));
    v.ref.erase(p1.first.base(), cend(v.ref));
    v.alt.erase(p1.second.base(), cend(v.alt));
    const auto p2 = std::mismatch(cbegin(v.ref), cend(v.ref), cbegin(v.alt), cend(v.alt));
    v.begin_pos += std::distance(cbegin(v.ref), p2.first);
    v.ref.erase(cbegin(v.ref), p2.first);
    v.alt.erase(cbegin(v.alt), p2.second);
}

void trim_reference(std::deque<Assembler::Variant>& variants)
{
    for (auto& v : variants) trim_reference(v);
}

bool is_complex(const Assembler::Variant& v) noexcept
{
    return (v.ref.size() > 1 && !v.alt.empty()) || (v.alt.size() > 1 && !v.ref.empty());
}

bool is_decomposable(const Assembler::Variant& v) noexcept
{
    return is_complex(v) && !is_inversion(v);
}

auto partition_decomposable(std::deque<Assembler::Variant>& variants)
{
    return std::stable_partition(std::begin(variants), std::end(variants),
                                 [] (const auto& candidate) { return !is_decomposable(candidate); });
}

bool is_mnv(const Assembler::Variant& v) noexcept
{
    return v.ref.size() == v.alt.size()
           && (v.ref.size() <= 2
               || std::equal(std::next(std::cbegin(v.ref)), std::prev(std::cend(v.ref)), std::next(std::cbegin(v.alt))));
}

auto split_mnv(Assembler::Variant&& mnv)
{
    assert(mnv.ref.size() > 1 && mnv.alt.size() > 1);
    assert(mnv.ref.front() != mnv.alt.front() && mnv.ref.back() != mnv.alt.back());
    std::vector<Assembler::Variant> result {};
    result.reserve(4);
    // Need to allocate new memory for all but the last SNV
    result.emplace_back(mnv.begin_pos, mnv.ref.front(), mnv.alt.front());
    const auto first_ref_itr       = std::cbegin(mnv.ref);
    const auto penultimate_ref_itr = std::prev(std::cend(mnv.ref));
    const auto first_alt_itr       = std::cbegin(mnv.alt);
    const auto penultimate_alt_itr = std::prev(std::cend(mnv.alt));
    auto p = std::mismatch(std::next(first_ref_itr), penultimate_ref_itr, std::next(first_alt_itr));
    while (p.first != penultimate_ref_itr) {
        assert(p.first < penultimate_ref_itr && p.second < penultimate_alt_itr);
        result.emplace_back(mnv.begin_pos + std::distance(first_ref_itr, p.first), *p.first, *p.second);
        p = std::mismatch(std::next(p.first), penultimate_ref_itr, std::next(p.second));
    }
    // So just need to remove the unwanted sequence from the last one
    const auto last_snv_begin = mnv.begin_pos + mnv.ref.size() - 1;
    mnv.ref.erase(first_ref_itr, penultimate_ref_itr);
    mnv.alt.erase(first_alt_itr, penultimate_alt_itr);
    result.emplace_back(last_snv_begin, std::move(mnv.ref), std::move(mnv.alt));
    return result;
}

struct Repeat : public Mappable<Repeat>
{
    ContigRegion region;
    unsigned period;
    Assembler::NucleotideSequence::const_iterator begin_itr, end_itr;
    const auto& mapped_region() const noexcept { return region; }
    auto begin() const noexcept { return begin_itr; }
    auto end() const noexcept { return end_itr; }
    Repeat() = default;
    Repeat(const tandem::Repeat& repeat, const Assembler::NucleotideSequence& sequence) noexcept
    : region {repeat.pos, repeat.pos + repeat.length}
    , period {repeat.period}
    , begin_itr {std::next(std::cbegin(sequence), repeat.pos)}
    , end_itr {std::next(begin_itr, repeat.length)}
    {}
};

auto find_repeats(Assembler::NucleotideSequence& sequence, const unsigned max_period = 5)
{
    sequence.push_back('$');
    auto repeats = tandem::extract_exact_tandem_repeats(sequence, 1, max_period);
    sequence.pop_back();
    std::vector<Repeat> result(repeats.size());
    std::transform(std::cbegin(repeats), std::cend(repeats), std::begin(result),
                   [&] (auto repeat) { return Repeat {repeat, sequence}; });
    std::sort(std::begin(result), std::end(result));
    return result;
}

struct VariantReference : public Mappable<VariantReference>
{
    using Position = ContigRegion::Position;
    ContigRegion region;
    std::reference_wrapper<Assembler::Variant> variant;
    const auto& mapped_region() const noexcept { return region; }
    VariantReference(Assembler::Variant& v)
    : region {static_cast<Position>(v.begin_pos), static_cast<Position>(v.begin_pos + v.ref.size())}
    , variant {v}
    {}
    const auto& ref() const noexcept { return variant.get().ref; }
    auto& ref() noexcept { return variant.get().ref; }
    const auto& alt() const noexcept { return variant.get().alt; }
    auto& alt() noexcept { return variant.get().alt; }
};

bool matches_rhs(const Repeat& repeat, const Assembler::NucleotideSequence& sequence) noexcept
{
    if (sequence.size() < repeat.period) return false;
    if (sequence.size() == repeat.period) {
        return std::equal(std::cbegin(sequence), std::cend(sequence), std::cbegin(repeat));
    } else if (utils::is_tandem_repeat(sequence, repeat.period)) {
        assert(std::distance(std::cbegin(repeat), std::cend(repeat)) >= 2 * repeat.period);
        const auto repeat_match_end_itr = std::next(std::cbegin(repeat), 2 * repeat.period);
        auto match_itr = std::search(std::cbegin(repeat), repeat_match_end_itr,
                                     std::cbegin(sequence), std::next(std::cbegin(sequence), repeat.period));
        return match_itr != repeat_match_end_itr;
    } else {
        return false;
    }
}

template <typename Range>
auto rotate_right(Range& range, std::size_t n)
{
    return std::rotate(std::rbegin(range), std::next(std::rbegin(range), n), std::rend(range));
}

void complete_partial_ref_repeat(Assembler::Variant& v, const Repeat& repeat)
{
    assert(v.ref.size() >= repeat.period);
    const auto partial_repeat_len = v.ref.size() % repeat.period;
    if (partial_repeat_len > 0) {
        const auto num_remaining_repeat_bases = repeat.period - partial_repeat_len;
        v.ref.reserve(v.ref.size() + num_remaining_repeat_bases);
        v.alt.reserve(v.alt.size() + num_remaining_repeat_bases);
        const auto rest_repeat_begin_itr = std::next(std::cbegin(v.ref), partial_repeat_len);
        const auto rest_repeat_end_itr   = std::next(rest_repeat_begin_itr, num_remaining_repeat_bases);
        v.ref.insert(std::cend(v.ref), rest_repeat_begin_itr, rest_repeat_end_itr);
        v.alt.insert(std::cend(v.alt), rest_repeat_begin_itr, rest_repeat_end_itr);
        rotate_right(v.ref, num_remaining_repeat_bases);
        rotate_right(v.alt, num_remaining_repeat_bases);
        v.begin_pos -= num_remaining_repeat_bases;
    }
}

void complete_partial_alt_repeat(Assembler::Variant& v, const Repeat& repeat)
{
    assert(v.alt.size() >= repeat.period);
    const auto partial_repeat_len = v.alt.size() % repeat.period;
    if (partial_repeat_len > 0) {
        const auto num_remaining_repeat_bases = repeat.period - partial_repeat_len;
        v.ref.reserve(v.ref.size() + num_remaining_repeat_bases);
        v.alt.reserve(v.alt.size() + num_remaining_repeat_bases);
        const auto rest_repeat_begin_itr = std::next(std::cbegin(v.alt), partial_repeat_len);
        const auto rest_repeat_end_itr   = std::next(rest_repeat_begin_itr, num_remaining_repeat_bases);
        v.ref.insert(std::cend(v.ref), rest_repeat_begin_itr, rest_repeat_end_itr);
        v.alt.insert(std::cend(v.alt), rest_repeat_begin_itr, rest_repeat_end_itr);
    }
}

std::vector<Assembler::Variant> try_to_split_repeats(Assembler::Variant& v, const std::vector<Repeat>& ref_repeats)
{
    VariantReference v_ref {v};
    const auto ref_repeat_itr = max_overlapped(ref_repeats, v_ref);
    if (ref_repeat_itr == std::cend(ref_repeats)) return {};
    const auto& ref_repeat = *ref_repeat_itr;
    if (!contains(ref_repeat, v_ref)) return {};
    if (v.ref.size() < 2 * ref_repeat.period) return {};
    const auto ref_repeat_is_lhs = left_overhang_size(ref_repeat, v_ref) > ref_repeat.period;
    if (ref_repeat_is_lhs) {
        if (ends_before(ref_repeat, v_ref)) return {};
        complete_partial_ref_repeat(v, ref_repeat);
    } else {
        auto alt_repeat_ritr = std::make_reverse_iterator(ref_repeat_itr);
        auto alt_repeat_match_ritr = std::crend(ref_repeats);
        for (; alt_repeat_ritr != std::crend(ref_repeats); ++alt_repeat_ritr) {
            if (is_before(*alt_repeat_ritr, ref_repeat)) break;
            if (matches_rhs(*alt_repeat_ritr, v.alt)) alt_repeat_match_ritr = alt_repeat_ritr;
        }
        if (alt_repeat_match_ritr == std::crend(ref_repeats)) return {};
        if (v.alt.size() < alt_repeat_match_ritr->period) return {};
        complete_partial_alt_repeat(v, *alt_repeat_match_ritr);
    }
    Assembler::Variant deletion {v.begin_pos, std::move(v.ref), ""}, insertion {v.begin_pos, "", std::move(v.alt)};
    if (ref_repeat_is_lhs) insertion.begin_pos += deletion.ref.size();
    return {std::move(deletion), std::move(insertion)};
}

auto extract_variants(const Assembler::NucleotideSequence& ref, const Assembler::NucleotideSequence& alt,
                      const CigarString& cigar, std::size_t ref_offset)
{
    std::vector<Assembler::Variant> result {};
    result.reserve(cigar.size());
    auto ref_itr = std::cbegin(ref);
    auto alt_itr = std::cbegin(alt);
    for (const auto& op : cigar) {
        using Flag = CigarOperation::Flag;
        using NucleotideSequence = Assembler::NucleotideSequence;
        switch(op.flag()) {
            case Flag::sequenceMatch:
            {
                ref_offset += op.size();
                ref_itr += op.size();
                alt_itr += op.size();
                break;
            }
            case Flag::substitution:
            {
                const auto next_ref_itr = std::next(ref_itr, op.size());
                std::transform(ref_itr, next_ref_itr, alt_itr, std::back_inserter(result),
                               [&ref_offset] (const auto ref, const auto alt) {
                                   return Assembler::Variant {ref_offset++, ref, alt};
                               });
                ref_itr = next_ref_itr;
                alt_itr += op.size();
                break;
            }
            case Flag::insertion:
            {
                const auto next_alt_itr = std::next(alt_itr, op.size());
                result.emplace_back(ref_offset, "", NucleotideSequence {alt_itr, next_alt_itr});
                alt_itr = next_alt_itr;
                break;
            }
            case Flag::deletion:
            {
                const auto next_ref_itr = std::next(ref_itr, op.size());
                result.emplace_back(ref_offset, NucleotideSequence {ref_itr, next_ref_itr}, "");
                ref_offset += op.size();
                ref_itr = next_ref_itr;
                break;
            }
            default:
                throw std::runtime_error {"LocalReassembler: unexpected cigar op"};
        }
        assert(ref_itr <= std::cend(ref) && alt_itr <= std::cend(alt));
    }
    return result;
}

auto align(const Assembler::Variant& v, const Model& model)
{
    return align(v.ref, v.alt, model).cigar;
}

struct VariantDecompositionConfig
{
    enum class ComplexAction { decompose, mnv, ignore } complex;
};

auto count_variant_types(const CigarString& cigar) noexcept
{
    bool has_snv {false}, has_insertion {false}, has_deletion {false};
    for (const auto& op : cigar) {
        switch (op.flag()) {
            case CigarOperation::Flag::substitution: has_snv = true; break;
            case CigarOperation::Flag::insertion: has_insertion = true; break;
            case CigarOperation::Flag::deletion: has_deletion = true; break;
            default: break;
        }
    }
    return has_snv + has_insertion + has_deletion;
}

bool is_complex_alignment(const CigarString& cigar, const Assembler::Variant& v) noexcept
{
    const auto min_allele_size = std::min(v.ref.size(), v.alt.size());
    const auto num_variant_types = count_variant_types(cigar);
    return (min_allele_size > 5 && cigar.size() >= min_allele_size && num_variant_types > 1)
           || (min_allele_size > 8 && cigar.size() > 2 * min_allele_size / 3 && num_variant_types > 1)
           || (min_allele_size > 20 && cigar.size() > min_allele_size / 2 && num_variant_types > 2)
           || (min_allele_size > 50 && cigar.size() > 2 * min_allele_size / 3 && num_variant_types > 2)
           || (min_allele_size > 100 && cigar.size() > min_allele_size / 3 && num_variant_types > 2);
}

bool is_good_alignment(const CigarString& cigar, const Assembler::Variant& v) noexcept
{
    return !is_complex_alignment(cigar, v);
}

bool is_unaligned(const CigarString& cigar) noexcept
{
    return cigar.size() == 2 && (is_deletion(cigar.front()) && is_insertion(cigar.back()));
}

std::vector<Assembler::Variant>
decompose_with_aligner(Assembler::Variant v, const Model& model, const VariantDecompositionConfig config)
{
    const auto cigar = align(v, model);
    if (is_good_alignment(cigar, v)) {
        auto result = extract_variants(v.ref, v.alt, cigar, v.begin_pos);
        if (is_unaligned(cigar)) {
            // If the sequences don't align then the insertion (variant sequence)
            // could be before or after the deletion (reference).
            assert(result.size() == 2);
            Assembler::Variant insertion {result.back()};
            insertion.begin_pos = result.front().begin_pos;
            result.push_back(std::move(insertion));
        }
        return result;
    } else if (config.complex == VariantDecompositionConfig::ComplexAction::mnv) {
        return {std::move(v)};
    } else {
        return {};
    }
}

auto decompose_with_aligner(Assembler::Variant v, const VariantDecompositionConfig config)
{
    Model model {4, -6, -8, -1};
    return decompose_with_aligner(std::move(v), model, config);
}

bool is_fully_decomposed(const Assembler::Variant v) noexcept
{
    return v.ref.empty() && v.alt.empty();
}

std::vector<Assembler::Variant> 
decompose_complex(Assembler::Variant v, const std::vector<Repeat>& reference_repeats, const VariantDecompositionConfig config)
{
    auto result = try_to_split_repeats(v, reference_repeats);
    if (result.empty() || !is_fully_decomposed(v)) {
        utils::append(decompose_with_aligner(std::move(v), config), result);
    }
    return result;
}

std::vector<Assembler::Variant> 
decompose(Assembler::Variant v, const std::vector<Repeat>& reference_repeats, const VariantDecompositionConfig config)
{
    if (is_mnv(v)) {
        return split_mnv(std::move(v));
    } else {
        return decompose_complex(std::move(v), reference_repeats, config);
    }
}

struct VariantLess
{
    using Variant = Assembler::Variant;
    bool operator()(const Variant& lhs, const Variant& rhs) const noexcept
    {
        if (lhs.begin_pos == rhs.begin_pos) {
            if (lhs.ref.size() == rhs.ref.size()) {
                return lhs.alt < rhs.alt;
            }
            return lhs.ref.size() < rhs.ref.size();
        }
        return lhs.begin_pos < rhs.begin_pos;
    }
};

using VariantIterator = std::deque<Assembler::Variant>::iterator;

auto decompose(VariantIterator first, VariantIterator last, const std::vector<Repeat>& reference_repeats, const VariantDecompositionConfig config)
{
    using std::begin; using std::end; using std::make_move_iterator;
    std::deque<Assembler::Variant> result {};
    std::for_each(make_move_iterator(first), make_move_iterator(last), [&] (auto&& complex) {
        utils::append(decompose(std::move(complex), reference_repeats, config), result);
    });
    std::sort(begin(result), end(result), VariantLess {});
    result.erase(std::unique(begin(result), end(result)), end(result));
    return result;
}

void merge(std::deque<Assembler::Variant>&& decomposed, std::deque<Assembler::Variant>& variants,
           std::deque<Assembler::Variant>::iterator first_complex)
{
    using std::begin; using std::end; using std::make_move_iterator;
    assert(!decomposed.empty());
    assert(begin(variants) <= first_complex && first_complex <= end(variants));
    const auto num_complex = static_cast<std::size_t>(std::distance(first_complex, end(variants)));
    // The variants in [first_complex, end(variants)) where moved from so can now be assigned to
    if (decomposed.size() <= num_complex) {
        variants.erase(std::move(begin(decomposed), end(decomposed), first_complex), end(variants));
    } else {
        const auto last_assignable = std::next(begin(decomposed), num_complex);
        assert(last_assignable <= end(decomposed));
        std::move(begin(decomposed), last_assignable, first_complex);
        first_complex = variants.insert(end(variants),
                                        make_move_iterator(last_assignable),
                                        make_move_iterator(end(decomposed)));
        first_complex -= num_complex;
    }
    assert(begin(variants) <= first_complex && first_complex <= end(variants));
    std::inplace_merge(begin(variants), first_complex, end(variants), VariantLess {});
}

void decompose(std::deque<Assembler::Variant>& variants, const ReferenceGenome::GeneticSequence& reference, const VariantDecompositionConfig config)
{
    auto tmp = reference;
    auto reference_repeats = find_repeats(tmp);
    const auto first_decomposable = partition_decomposable(variants);
    if (first_decomposable != std::end(variants)) {
        auto decomposed = decompose(first_decomposable, std::end(variants), reference_repeats, config);
        if (!decomposed.empty()) {
            merge(std::move(decomposed), variants, first_decomposable);
        } else {
            variants.erase(first_decomposable, std::end(variants));
        }
    }
}

void add_to_mapped_variants(std::deque<Assembler::Variant>&& variants, std::deque<Variant>& result,
                            const GenomicRegion& assemble_region)
{
    for (auto& variant : variants) {
        result.emplace_back(contig_name(assemble_region), assemble_region.begin() + variant.begin_pos,
                            std::move(variant.ref), std::move(variant.alt));
    }
}

void remove_large_deletions(std::deque<Assembler::Variant>& variants, const unsigned max_size)
{
    variants.erase(std::remove_if(std::begin(variants), std::end(variants),
                                  [=] (const auto& variant) { return variant.ref.size() >= max_size && variant.alt.empty(); }),
                   std::end(variants));
}

LocalReassembler::AssemblerStatus
LocalReassembler::try_assemble_region(Assembler& assembler, const NucleotideSequence& reference_sequence,
                                      const GenomicRegion& assemble_region, std::deque<Variant>& result) const
{
    assert(assembler.is_unique_reference());
    assembler.try_recover_dangling_branches();
    assembler.prune(min_kmer_observations_);
    auto status = AssemblerStatus::success;
    if (!assembler.is_acyclic()) {
        if (cycle_tolerance_ == Options::CyclicGraphTolerance::none) {
            const auto num_samples = read_buffer_.size();
            if (num_samples > 1) {
                const auto cyclic_samples = assembler.find_cyclic_samples();
                if (num_samples == cyclic_samples.size()) {
                    return AssemblerStatus::failed;
                } else {
                    assembler.clear(cyclic_samples);
                    assembler.prune(min_kmer_observations_);
                    if (assembler.is_acyclic()) {
                        status = AssemblerStatus::partial_success;
                    } else {
                        return AssemblerStatus::failed;
                    }
                }
            } else {
                return AssemblerStatus::failed;
            }
        } else {
            assembler.remove_nonreference_cycles();
            status = AssemblerStatus::partial_success;
        }
    }
    assembler.cleanup();
    if (assembler.is_empty() || assembler.is_all_reference()) {
        return status;
    }
    const auto min_bubble_score = [&] (const auto ref_head_idx, const auto ref_tail_idx) -> double {
        const auto ref_bubble_region = expand_rhs(shift(head_region(assemble_region), ref_head_idx), ref_tail_idx - ref_head_idx);
        return calculate_min_bubble_score(ref_bubble_region);
    };
    auto variants = assembler.extract_variants(max_bubbles_, min_bubble_score);
    assembler.clear();
    if (!variants.empty()) {
        trim_reference(variants);
        std::sort(std::begin(variants), std::end(variants), VariantLess {});
        variants.erase(std::unique(std::begin(variants), std::end(variants)), std::end(variants));
        VariantDecompositionConfig decomposition_config {};
        if (status == AssemblerStatus::success || cycle_tolerance_ == Options::CyclicGraphTolerance::high) {
            decomposition_config.complex = VariantDecompositionConfig::ComplexAction::mnv;
        } else {
            decomposition_config.complex = VariantDecompositionConfig::ComplexAction::ignore;
        }
        decompose(variants, reference_sequence, decomposition_config);
        if (status == AssemblerStatus::partial_success) {
            // TODO: Some false positive large deletions are being generated for small kmer sizes.
            // Until Assembler is better able to remove these automatically, filter them here.
            if (assembler.kmer_size() <= 10) {
                remove_large_deletions(variants, 100);
            } else if (assembler.kmer_size() <= 15) {
                remove_large_deletions(variants, 150);
            } else if (assembler.kmer_size() <= 20) {
                remove_large_deletions(variants, 200);
            } else if (assembler.kmer_size() <= 30) {
                remove_large_deletions(variants, 250);
            }
        }
        add_to_mapped_variants(std::move(variants), result, assemble_region);
    }
    return status;
}

double LocalReassembler::calculate_min_bubble_score(const GenomicRegion& bubble_region) const
{
    ReadBaseCountMap read_counts {};
    read_counts.reserve(read_buffer_.size());
    for (const auto& p : read_buffer_) {
        read_counts.emplace(p.first, count_base_pairs(p.second, bubble_region));
    }
    return min_bubble_score_(bubble_region, read_counts);
}

double DepthBasedBubbleScoreSetter::operator()(const GenomicRegion& region, const LocalReassembler::ReadBaseCountMap& read_counts) const
{
    auto result = min_score_;
    for (const auto& p : read_counts) {
        const auto mean_depth = p.second / size(region);
        result = std::max(result, min_allele_frequency_ * mean_depth);
    }
    return result;
}

DepthBasedBubbleScoreSetter::DepthBasedBubbleScoreSetter(double min_score, double min_allele_frequency)
: min_score_ {min_score}
, min_allele_frequency_ {min_allele_frequency}
{}

} // namespace coretools
} // namespace octopus
