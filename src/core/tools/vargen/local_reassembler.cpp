// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "local_reassembler.hpp"

#include <algorithm>
#include <iterator>
#include <deque>
#include <stdexcept>
#include <cassert>

#include <config/common.hpp>
#include <basics/cigar_string.hpp>
#include <basics/aligned_read.hpp>
#include <concepts/mappable_range.hpp>
#include <utils/mappable_algorithms.hpp>
#include <utils/sequence_utils.hpp>
#include <utils/append.hpp>
#include <io/reference/reference_genome.hpp>
#include <logging/logging.hpp>

#include "utils/global_aligner.hpp"

namespace octopus { namespace coretools {

LocalReassembler::LocalReassembler(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, default_kmer_sizes_ {std::move(options.kmer_sizes)}
, fallback_kmer_sizes_ {}
, bin_size_ {1000}
, bins_ {}
, mask_threshold_ {options.mask_threshold}
, min_supporting_reads_ {options.min_supporting_reads}
, max_variant_size_ {options.max_variant_size}
{
    using std::begin; using std::end;
    
    if (default_kmer_sizes_.empty()) return;
    
    std::sort(begin(default_kmer_sizes_), end(default_kmer_sizes_));
    
    default_kmer_sizes_.erase(std::unique(begin(default_kmer_sizes_), end(default_kmer_sizes_)),
                              end(default_kmer_sizes_));
    fallback_kmer_sizes_.resize(options.num_fallbacks);
    
    auto k = default_kmer_sizes_.back();
    std::generate_n(begin(fallback_kmer_sizes_), options.num_fallbacks,
                    [&] () {
                        k += options.fallback_interval_size;
                        return k;
                    });
}

std::unique_ptr<VariantGenerator> LocalReassembler::do_clone() const
{
    return std::make_unique<LocalReassembler>(*this);
}

LocalReassembler::Bin::Bin(GenomicRegion region)
:
region {std::move(region)}
{}

const GenomicRegion& LocalReassembler::Bin::mapped_region() const noexcept
{
    return region;
}

void LocalReassembler::Bin::insert(const AlignedRead& read)
{
    read_sequences.emplace_back(read.sequence());
}

void LocalReassembler::Bin::insert(const NucleotideSequence& sequence)
{
    read_sequences.emplace_back(sequence);
}

void LocalReassembler::Bin::clear() noexcept
{
    read_sequences.clear();
    read_sequences.shrink_to_fit();
}

bool LocalReassembler::Bin::empty() const noexcept
{
    return read_sequences.empty();
}

bool LocalReassembler::do_requires_reads() const noexcept
{
    return true;
}

namespace {

bool has_low_quality_match(const AlignedRead& read, const AlignedRead::BaseQuality good_quality) noexcept
{
    if (good_quality == 0) return false;
    auto quality_itr = std::cbegin(read.qualities());
    return std::any_of(std::cbegin(read.cigar()), std::cend(read.cigar()),
                       [&] (const auto& op) {
                           if (is_match(op)) {
                               auto result = std::any_of(quality_itr, std::next(quality_itr, op.size()),
                                                         [=] (auto q) { return q < good_quality; });
                               std::advance(quality_itr, op.size());
                               return result;
                           } else if (op.advances_sequence()) {
                               std::advance(quality_itr, op.size());
                           }
                           return false;
                       });
}

using ExpandedCigarString = std::vector<CigarOperation::Flag>;

auto expand_cigar(const CigarString& cigar, const std::size_t size_hint = 0)
{
    ExpandedCigarString result {};
    result.reserve(size_hint);
    
    for (const auto& op : cigar) {
        utils::append(result, op.size(), op.flag());
    }
    
    return result;
}

auto expand_cigar(const AlignedRead& read)
{
    return expand_cigar(read.cigar(), sequence_size(read));
}

auto find_first_sequence_op(const ExpandedCigarString& cigar) noexcept
{
    return std::find_if_not(std::cbegin(cigar), std::cend(cigar),
                            [] (auto op) { return op == CigarOperation::Flag::hardClipped; });
}

bool is_match(const CigarOperation::Flag op) noexcept
{
    switch (op) {
        case CigarOperation::Flag::alignmentMatch:
        case CigarOperation::Flag::sequenceMatch:
        case CigarOperation::Flag::substitution: return true;
        default: return false;
    }
}

auto transform_low_quality_matches_to_reference(AlignedRead::NucleotideSequence read_sequence,
                                                const AlignedRead::BaseQualityVector& base_qualities,
                                                const AlignedRead::NucleotideSequence& reference_sequence,
                                                const ExpandedCigarString& cigar,
                                                const AlignedRead::BaseQuality min_quality)
{
    auto ref_itr = std::cbegin(reference_sequence);
    auto cigar_itr = find_first_sequence_op(cigar);
    std::transform(std::cbegin(read_sequence), std::cend(read_sequence), std::cbegin(base_qualities),
                   std::begin(read_sequence), [&] (const auto read_base, const auto base_quality) {
                       using Flag = CigarOperation::Flag;
                       const auto op = *cigar_itr++;
                       if (!is_match(op)) {
                           if (op != Flag::insertion) ++ref_itr;
                           // Deletions are excess reference sequence so we need to move the
                           // reference iterator to the next nondeleted read base
                           while (cigar_itr != std::cend(cigar) && *cigar_itr == Flag::deletion) {
                               ++cigar_itr;
                               ++ref_itr;
                           }
                           return read_base;
                       }
                       const auto ref_base = *ref_itr++; // Don't forget to increment ref_itr!
                       return base_quality >= min_quality ? read_base : ref_base;
                   });
    return read_sequence;
}

AlignedRead::NucleotideSequence
transform_low_quality_matches_to_reference(const AlignedRead& read,
                                           const AlignedRead::BaseQuality min_quality,
                                           const ReferenceGenome& reference)
{
    return transform_low_quality_matches_to_reference(read.sequence(), read.qualities(),
                                                      reference.fetch_sequence(mapped_region(read)),
                                                      expand_cigar(read), min_quality);
}

} // namespace

template <typename C, typename R>
auto overlapped_bins(C& bins, const R& read)
{
    return bases(overlap_range(std::begin(bins), std::end(bins), read, BidirectionallySortedTag {}));
}

void LocalReassembler::do_add_read(const AlignedRead& read)
{
    prepare_bins_to_insert(read);
    
    auto active_bins = overlapped_bins(bins_, read);
    
    assert(!active_bins.empty());
    
    if (!has_low_quality_match(read, mask_threshold_)) {
        for (auto& bin : active_bins) {
            bin.insert(read);
        }
    } else {
        auto masked_sequence = transform_low_quality_matches_to_reference(read, mask_threshold_,
                                                                          reference_);
        masked_sequence_buffer_.emplace_back(std::move(masked_sequence));
        for (auto& bin : active_bins) {
            bin.insert(std::cref(masked_sequence_buffer_.back()));
        }
    }
}

void LocalReassembler::do_add_reads(VectorIterator first, VectorIterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { do_add_read(read); });
}

void LocalReassembler::do_add_reads(FlatSetIterator first, FlatSetIterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { do_add_read(read); });
}

template <typename Container>
void remove_nonoverlapping(Container& candidates, const GenomicRegion& region)
{
    const auto it = std::remove_if(std::begin(candidates), std::end(candidates),
                                   [&region] (const Variant& candidate) {
                                       return !overlaps(candidate, region);
                                   });
    candidates.erase(it, std::end(candidates));
}

template <typename L>
void log_success(L& log, const char* type, const unsigned k)
{
    if (log) stream(*log, 8) << type << " assembler with kmer size " << k << " completed";
}

template <typename L>
void log_failure(L& log, const char* type, const unsigned k)
{
    if (log) stream(*log, 8) << type << " assembler with kmer size " << k << " failed";
}

auto extract_unique(std::deque<Variant>&& variants)
{
    std::vector<Variant> result {std::make_move_iterator(std::begin(variants)),
        std::make_move_iterator(std::end(variants))};
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

void remove_oversized(std::vector<Variant>& variants, const Variant::MappingDomain::Size max_size)
{
    variants.erase(std::remove_if(std::begin(variants), std::end(variants),
                                  [max_size] (const auto& variant) {
                                        return region_size(variant) > max_size;
                                  }),
                   std::end(variants));
}

std::vector<Variant> LocalReassembler::do_generate_variants(const GenomicRegion& region)
{
    if (bins_.empty()) return {};
    
    auto active_bins = overlapped_bins(bins_, region);
    
    std::deque<Variant> candidates {};
    
    for (Bin& bin : active_bins) {
        if (bin.read_sequences.empty()) continue;
        
        if (debug_log_) {
            stream(*debug_log_) << "Assembling " << bin.read_sequences.size()
                                << " reads in bin " << mapped_region(bin);
        }
        
        unsigned num_defaults_unsuccessful {0};
        
        for (const auto k : default_kmer_sizes_) {
            const auto success = assemble_bin(k, bin, candidates);
            if (success) {
                log_success(debug_log_, "Default", k);
            } else {
                log_failure(debug_log_, "Default", k);
                ++num_defaults_unsuccessful;
            }
        }
        
        if (num_defaults_unsuccessful == default_kmer_sizes_.size()) {
            for (const auto k : fallback_kmer_sizes_) {
                const auto success = assemble_bin(k, bin, candidates);
                if (success) {
                    log_success(debug_log_, "Fallback", k);
                    break;
                } else {
                    log_failure(debug_log_, "Fallback", k);
                }
            }
        }
        
        bin.clear();
    }
    
    auto result = extract_unique(std::move(candidates));
    remove_oversized(result, max_variant_size_);
    remove_nonoverlapping(result, region); // as we expanded original region
    return result;
}

void LocalReassembler::do_clear() noexcept
{
    bins_.clear();
    bins_.shrink_to_fit();
    masked_sequence_buffer_.clear();
    masked_sequence_buffer_.shrink_to_fit();
}

std::string LocalReassembler::name() const
{
    return "Local reassembly";
}

void LocalReassembler::prepare_bins_to_insert(const AlignedRead& read)
{
    const auto& read_region = mapped_region(read);
    
    if (bins_.empty()) {
        if (region_size(read_region) > bin_size_) {
            for (auto subregion : decompose(read_region, bin_size_)) {
                bins_.emplace_back(std::move(subregion));
            }
            bins_.emplace_back(shift(mapped_region(bins_.back()), bin_size_));
        } else if (region_size(read_region) == bin_size_) {
            bins_.emplace_back(read_region);
        } else {
            bins_.emplace_back(expand_rhs(head_region(read_region), bin_size_));
        }
    } else if (!contains(encompassing_region(bins_.front(), bins_.back()), read_region)) {
        while (begins_before(read_region, bins_.front())) {
            bins_.emplace_front(shift(mapped_region(bins_.front()), -bin_size_));
        }
        while (ends_before(bins_.back(), read_region)) {
            bins_.emplace_back(shift(mapped_region(bins_.back()), bin_size_));
        }
    }
    
    assert(contains(encompassing_region(bins_.front(), bins_.back()), read_region));
}

GenomicRegion LocalReassembler::propose_assembler_region(const GenomicRegion& input_region,
                                                         unsigned kmer_size) const
{
    return expand(input_region, kmer_size);
}

void trim_reference(Assembler::Variant& v)
{
    const auto p = std::mismatch(std::begin(v.ref), std::end(v.ref),
                                 std::begin(v.alt), std::end(v.alt));

    v.begin_pos += std::distance(std::begin(v.ref), p.first);
    
    v.ref.erase(std::begin(v.ref), p.first);
    v.alt.erase(std::begin(v.alt), p.second);
    
    const auto p2 = std::mismatch(std::rbegin(v.ref), std::rend(v.ref),
                                  std::rbegin(v.alt), std::rend(v.alt));
    
    v.ref.erase(p2.first.base(), std::end(v.ref));
    v.alt.erase(p2.second.base(), std::end(v.alt));
}

template <typename Container>
void trim_reference(Container& variants)
{
    for (auto& v : variants) trim_reference(v);
}

bool is_complex(const Assembler::Variant& v)
{
    return (v.ref.size() > 1 && !v.alt.empty()) || (v.alt.size() > 1 && !v.ref.empty());
}

bool is_mnp(const Assembler::Variant& v)
{
    return v.ref.size() > 1 && v.ref.size() == v.alt.size();
}

auto split_mnp(Assembler::Variant&& v)
{
    using std::begin; using std::end; using std::next; using std::prev; using std::distance;
    
    std::vector<Assembler::Variant> result {};
    result.reserve(4);
    result.emplace_back(v.begin_pos, v.ref.front(), v.alt.front());
    
    auto p = std::mismatch(next(begin(v.ref)), prev(end(v.ref)), next(begin(v.alt)));
    
    while (p.first != prev(end(v.ref))) {
        result.emplace_back(v.begin_pos + distance(begin(v.ref), p.first), *p.first, *p.second);
        p = std::mismatch(next(p.first), prev(end(v.ref)), next(p.second));
    }
    
    const auto pos = v.begin_pos + v.ref.size() - 1;
    
    v.ref.erase(begin(v.ref), prev(end(v.ref)));
    v.alt.erase(begin(v.alt), prev(end(v.alt)));
    
    result.emplace_back(pos, std::move(v.ref), std::move(v.alt));
    
    return result;
}

auto extract_variants(const Assembler::NucleotideSequence& ref, const Assembler::NucleotideSequence& alt,
                      const CigarString& cigar, std::size_t ref_offset)
{
    std::vector<Assembler::Variant> result {};
    result.reserve(cigar.size());
    
    auto ref_it = std::cbegin(ref);
    auto alt_it = std::cbegin(alt);
    
    for (const auto& op : cigar) {
        using Flag = CigarOperation::Flag;
        using NucleotideSequence = Assembler::NucleotideSequence;
        
        switch(op.flag()) {
            case Flag::sequenceMatch:
            {
                ref_offset += op.size();
                ref_it += op.size();
                alt_it += op.size();
                break;
            }
            case Flag::substitution:
            {
                std::transform(ref_it, std::next(ref_it, op.size()), alt_it, std::back_inserter(result),
                               [&ref_offset] (const auto ref, const auto alt) {
                                   return Assembler::Variant {ref_offset++, ref, alt};
                               });
                ref_it += op.size();
                alt_it += op.size();
                break;
            }
            case Flag::insertion:
            {
                result.emplace_back(ref_offset, "", NucleotideSequence {alt_it, std::next(alt_it, op.size())});
                alt_it += op.size();
                break;
            }
            case Flag::deletion:
            {
                result.emplace_back(ref_offset, NucleotideSequence {ref_it, std::next(ref_it, op.size())}, "");
                ref_offset += op.size();
                ref_it += op.size();
                break;
            }
            default:
                throw std::runtime_error {"LocalReassembler: unexpected cigar op"};
        }
    }
    
    return result;
}

auto align(const Assembler::Variant& v)
{
    return parse_cigar(align(v.ref, v.alt).cigar);
}

auto decompose_complex_indel(Assembler::Variant&& v)
{
    const auto cigar = align(v);
    return extract_variants(v.ref, v.alt, cigar, v.begin_pos);
}

auto decompose_complex(Assembler::Variant&& v)
{
    return is_mnp(v) ? split_mnp(std::move(v)) : decompose_complex_indel(std::move(v));
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
            return lhs.ref.size() < rhs.alt.size();
        }
        return lhs.begin_pos < rhs.begin_pos;
    }
};

template <typename Container>
auto partition_complex(Container& variants)
{
    return std::stable_partition(std::begin(variants), std::end(variants),
                                 [] (const auto& candidate) {
                                     return !is_complex(candidate);
                                 });
}

template <typename InputIt>
auto decompose_complex(InputIt first_complex, InputIt last_complex)
{
    using std::begin; using std::end; using std::make_move_iterator;
    std::deque<Assembler::Variant> result {};
    std::for_each(make_move_iterator(first_complex), make_move_iterator(last_complex),
                  [&result] (auto&& complex) {
                      utils::append(decompose_complex(std::move(complex)), result);
                  });
    std::sort(begin(result), end(result), VariantLess {});
    result.erase(std::unique(begin(result), end(result)), end(result));
    return result;
}

template <typename Container1, typename Container2, typename Iterator>
void merge(Container1&& decomposed, Container2& variants, Iterator first_complex)
{
    using std::begin; using std::end; using std::make_move_iterator;
    const auto num_complex = std::distance(first_complex, end(variants));
    if (decomposed.size() < num_complex) {
        variants.erase(std::move(begin(decomposed), end(decomposed), first_complex), end(variants));
        std::inplace_merge(begin(variants), first_complex, end(variants), VariantLess {});
    } else {
        // The variants in [first_complex, end(variants)) where moved from so can now be assigned to
        const auto last_assignable = std::next(begin(decomposed), num_complex);
        std::move(begin(decomposed), last_assignable, first_complex);
        first_complex = variants.insert(end(variants),
                                        make_move_iterator(last_assignable),
                                        make_move_iterator(end(decomposed)));
        first_complex -= num_complex;
        std::inplace_merge(begin(variants), first_complex, end(variants), VariantLess {});
    }
}

template <typename Container>
void decompose_complex(Container& variants)
{
    const auto first_complex = partition_complex(variants);
    merge(decompose_complex(first_complex, std::end(variants)), variants, first_complex);
}

template <typename C1, typename C2>
void add_to_mapped_variants(C1& result, C2&& variants, const GenomicRegion& region)
{
    for (auto& variant : variants) {
        result.emplace_back(contig_name(region), region.begin() + variant.begin_pos,
                            std::move(variant.ref), std::move(variant.alt));
    }
}

bool LocalReassembler::assemble_bin(const unsigned kmer_size, const Bin& bin,
                                    std::deque<Variant>& result) const
{
    if (bin.empty()) return true;
    
    const auto assembler_region   = propose_assembler_region(bin.region, kmer_size);
    const auto reference_sequence = reference_.get().fetch_sequence(assembler_region);
    
    if (utils::has_ns(reference_sequence)) return false;
    
    Assembler assembler {kmer_size, reference_sequence};
    
    for (const auto& sequence : bin.read_sequences) {
        assembler.insert_read(sequence);
    }
    
    return try_assemble_region(assembler, reference_sequence, assembler_region, result);
}

bool LocalReassembler::try_assemble_region(Assembler& assembler,
                                           const NucleotideSequence& reference_sequence,
                                           const GenomicRegion& reference_region,
                                           std::deque<Variant>& result) const
{
    if (!assembler.prune(min_supporting_reads_)) {
        return false;
    }
    
    auto variants = assembler.extract_variants();
    
    assembler.clear();
    
    if (variants.empty()) {
        return true;
    }
    
    trim_reference(variants);
    decompose_complex(variants);
    add_to_mapped_variants(result, std::move(variants), reference_region);
    
    return true;
}

} // namespace coretools
} // namespace octopus
