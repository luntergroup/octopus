//
//  local_reassembler.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "local_reassembler.hpp"

#include <algorithm>
#include <iterator>
#include <deque>
#include <stdexcept>
#include <cassert>

#include "common.hpp"
#include "reference_genome.hpp"
#include "aligned_read.hpp"
#include "mappable_ranges.hpp"
#include "mappable_algorithms.hpp"
#include "sequence_utils.hpp"
#include "global_aligner.hpp"
#include "cigar_string.hpp"
#include "logging.hpp"

#include "timers.hpp" // BENCHMARK

namespace octopus { namespace coretools {

LocalReassembler::LocalReassembler(const ReferenceGenome& reference, Options options)
:
reference_ {reference},
default_kmer_sizes_ {std::move(options.kmer_sizes)},
fallback_kmer_sizes_ {},
bin_size_ {1000},
bins_ {},
mask_threshold_ {options.mask_threshold},
min_supporting_reads_ {options.min_supporting_reads},
max_variant_size_ {options.max_variant_size}
{
    using std::begin; using std::end;
    
    if (default_kmer_sizes_.empty()) {
        return;
    }
    
    std::sort(begin(default_kmer_sizes_), end(default_kmer_sizes_));
    
    default_kmer_sizes_.erase(std::unique(begin(default_kmer_sizes_), end(default_kmer_sizes_)),
                              end(default_kmer_sizes_));
    
    constexpr unsigned num_fallbacks {6};
    constexpr unsigned fallback_interval_size {10};
    
    fallback_kmer_sizes_.resize(num_fallbacks);
    
    auto k = default_kmer_sizes_.back();
    std::generate_n(begin(fallback_kmer_sizes_), num_fallbacks,
                    [&k, fallback_interval_size] () {
                        k += fallback_interval_size;
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

bool has_low_quality_match(const AlignedRead& read, const AlignedRead::BaseQuality good_quality)
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

namespace
{
    auto expand_cigar(const AlignedRead& read)
    {
        std::string result {};
        result.reserve(sequence_size(read));
        
        for (const auto& op : read.cigar()) {
            result.append(op.size(), op.flag());
        }
        
        return result;
    }
    
    bool is_match(const CigarOperation::Flag op)
    {
        switch (op) {
            case CigarOperation::ALIGNMENT_MATCH:
            case CigarOperation::SEQUENCE_MATCH:
            case CigarOperation::SUBSTITUTION: return true;
            default: return false;
        }
    }
}

AlignedRead::NucleotideSequence
transform_low_quality_matches_to_reference(const AlignedRead& read,
                                           const AlignedRead::BaseQuality min_quality,
                                           const ReferenceGenome& reference)
{
    auto result = read.sequence();
    
    const auto cigar = expand_cigar(read);
    
    auto cigar_op_itr = std::find_if_not(std::cbegin(cigar), std::cend(cigar),
                                         [] (auto op) { return op == CigarOperation::HARD_CLIPPED; });
    
    const auto ref_sequence = reference.fetch_sequence(mapped_region(read));
    
    auto ref_itr = std::cbegin(ref_sequence);
    
    std::transform(std::cbegin(result), std::cend(result), std::cbegin(read.qualities()),
                   std::begin(result), [&] (const char base, const auto quality) {
                       const auto cigar_op = *cigar_op_itr++;
                       
                       if (!is_match(cigar_op)) {
                           if (cigar_op != CigarOperation::INSERTION) ++ref_itr;
                           while (*cigar_op_itr == CigarOperation::DELETION) {
                               ++cigar_op_itr;
                               ++ref_itr;
                           }
                           return base;
                       }
                       
                       auto ref_base = *ref_itr++;
                       
                       return quality >= min_quality ? base : ref_base;
                   });
    
    return result;
}

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

std::vector<Variant>
LocalReassembler::do_generate_variants(const GenomicRegion& region)
{
    if (bins_.empty()) return {};
    
    auto active_bins = overlapped_bins(bins_, region);
    
    std::deque<Variant> candidates {};
    
    for (Bin& bin : active_bins) {
        if (bin.read_sequences.empty()) continue;
        
        if (debug_log_) {
            stream(*debug_log_) << "Assembling " << bin.read_sequences.size() << " reads in bin " << mapped_region(bin);
        }
        
        unsigned num_defaults_unsuccessful {0};
        
        for (const auto k : default_kmer_sizes_) {
            try {
                const auto success = assemble_bin(k, bin, candidates);
                
                if (success) {
                    log_success(debug_log_, "Default", k);
                } else {
                    log_failure(debug_log_, "Default", k);
                    ++num_defaults_unsuccessful;
                }
            } catch (std::exception& e) {
                break;
            }
        }
        
        if (num_defaults_unsuccessful == default_kmer_sizes_.size()) {
            for (const auto k : fallback_kmer_sizes_) {
                try {
                    const auto success = assemble_bin(k, bin, candidates);
                    
                    if (success) {
                        log_success(debug_log_, "Fallback", k);
                        break;
                    } else {
                        log_failure(debug_log_, "Fallback", k);
                    }
                } catch (std::exception& e) {
                    break;
                }
            }
        }
        
        bin.clear();
    }
    
    std::vector<Variant> result {std::make_move_iterator(std::begin(candidates)),
                                 std::make_move_iterator(std::end(candidates))};
    
    std::sort(std::begin(result), std::end(result));
    
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    
    result.erase(std::remove_if(std::begin(result), std::end(result),
                                [this] (const auto& variant) {
                                    return region_size(variant) > max_variant_size_;
                                }), std::end(result));
    
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

GenomicRegion
LocalReassembler::propose_assembler_region(const GenomicRegion& input_region,
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

std::vector<Assembler::Variant> split_complex(Assembler::Variant&& v)
{
    using std::begin; using std::end; using std::next; using std::prev; using std::distance;
    
    std::vector<Assembler::Variant> result {};
    
    if (is_mnp(v)) {
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
    
    auto cigar = parse_cigar(align(v.ref, v.alt).cigar);
    
    auto ref_it = std::cbegin(v.ref);
    auto alt_it = std::cbegin(v.alt);
    
    auto pos = v.begin_pos;
    
    result.reserve(cigar.size());
    
    for (const auto& op : cigar) {
        switch(op.flag()) {
            case CigarOperation::SEQUENCE_MATCH:
            {
                pos += op.size();
                ref_it += op.size();
                alt_it += op.size();
                break;
            }
            case CigarOperation::SUBSTITUTION:
            {
                std::transform(ref_it, std::next(ref_it, op.size()), alt_it,
                               std::back_inserter(result),
                               [&pos] (const auto ref, const auto alt) {
                                   return Assembler::Variant {pos++, ref, alt};
                               });
                ref_it += op.size();
                alt_it += op.size();
                break;
            }
            case CigarOperation::INSERTION:
            {
                result.emplace_back(pos, "", Assembler::NucleotideSequence {alt_it, std::next(alt_it, op.size())});
                alt_it += op.size();
                break;
            }
            case CigarOperation::DELETION:
            {
                result.emplace_back(pos, Assembler::NucleotideSequence {ref_it, std::next(ref_it, op.size())}, "");
                pos += op.size();
                ref_it += op.size();
                break;
            }
        }
    }
    
    return result;
}

template <typename Container>
void split_complex(Container& candidates)
{
    using std::begin; using std::end; using std::make_move_iterator;
    
    const auto it = std::stable_partition(begin(candidates), end(candidates),
                                          [] (const auto& candidate) {
                                              return !is_complex(candidate);
                                          });
    
    std::deque<Assembler::Variant> snps {};
    
    std::for_each(make_move_iterator(it), make_move_iterator(end(candidates)),
                  [&snps] (auto&& mnv) {
                      auto split = split_complex(std::move(mnv));
                      snps.insert(end(snps),
                                  make_move_iterator(begin(split)),
                                  make_move_iterator(end(split)));
                  });
    
    const auto variant_less = [] (const auto& lhs, const auto& rhs) {
        return lhs.begin_pos < rhs.begin_pos || lhs.alt < rhs.alt;
    };
    
    std::sort(begin(snps), end(snps), variant_less);
    
    snps.erase(std::unique(begin(snps), end(snps)), end(snps));
    
    const auto num_mnvs = std::distance(it, end(candidates));
    
    if (snps.size() < num_mnvs) {
        candidates.erase(std::move(begin(snps), end(snps), it), end(candidates));
        std::inplace_merge(begin(candidates), it, end(candidates), variant_less);
    } else {
        const auto it2 = std::next(begin(snps), num_mnvs);
        
        std::move(begin(snps), it2, it);
        
        const auto it3 = candidates.insert(end(candidates),
                                           make_move_iterator(it2),
                                           make_move_iterator(end(snps)));
        
        std::inplace_merge(begin(candidates), it3, end(candidates), variant_less);
    }
}

template <typename C1, typename C2>
void add_to_mapped_variants(C1& result, C2&& variants, const GenomicRegion& region)
{
    for (auto& variant : variants) {
        result.emplace_back(contig_name(region),
                            region.begin() + variant.begin_pos,
                            std::move(variant.ref),
                            std::move(variant.alt));
    }
}

bool LocalReassembler::assemble_bin(const unsigned kmer_size, const Bin& bin,
                                                      std::deque<Variant>& result) const
{
    if (bin.empty()) return true;
    
    const auto assembler_region = propose_assembler_region(bin.region, kmer_size);
    
    const auto reference_sequence = reference_.get().fetch_sequence(assembler_region);
    
    if (utils::has_ns(reference_sequence)) {
        throw std::runtime_error {"Bad reference"};
    }
    
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
    split_complex(variants);
    add_to_mapped_variants(result, std::move(variants), reference_region);
    
    return true;
}

} // namespace coretools
} // namespace octopus
