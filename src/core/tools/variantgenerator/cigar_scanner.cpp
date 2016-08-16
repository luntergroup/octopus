// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cigar_scanner.hpp"

#include <iterator>
#include <algorithm>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include <config/common.hpp>
#include <basics/aligned_read.hpp>
#include <basics/cigar_string.hpp>
#include <io/reference/reference_genome.hpp>
#include <concepts/mappable_range.hpp>
#include <utils/mappable_algorithms.hpp>
#include <logging/logging.hpp>

#include <iostream> // DEBUG

namespace octopus { namespace coretools {

std::unique_ptr<VariantGenerator> CigarScanner::do_clone() const
{
    return std::make_unique<CigarScanner>(*this);
}

CigarScanner::CigarScanner(const ReferenceGenome& reference, Options options)
:
reference_ {reference},
options_ {options},
match_ {[] (const Variant& lhs, const Variant& rhs) -> bool {
    return (is_insertion(lhs) && is_same_region(lhs, rhs)) || lhs == rhs;
}},
candidates_ {},
max_seen_candidate_size_ {}
{}

bool CigarScanner::do_requires_reads() const noexcept
{
    return true;
}

namespace {

template <typename Sequence>
bool is_good_sequence(const Sequence& sequence) noexcept
{
    return std::none_of(std::cbegin(sequence), std::cend(sequence),
                        [] (auto base) { return base == 'N'; });
}

template <typename Sequence, typename P, typename S>
Sequence splice(const Sequence& sequence, const P pos, const S size)
{
    const auto it = std::next(std::cbegin(sequence), pos);
    return Sequence {it, std::next(it, size)};
}

template <typename Q, typename P, typename S>
bool is_all_good_quality(const Q& qualities, const P pos, const S size,
                         const typename Q::value_type min_quality)
{
    const auto it = std::next(std::cbegin(qualities), pos);
    return std::all_of(it, std::next(it, size),
                       [min_quality] (const auto quality) {
                           return quality >= min_quality;
                       });
}

template <typename Q, typename P, typename S>
auto count_bad_qualities(const Q& qualities, const P pos, const S size,
                         const typename Q::value_type min_quality)
{
    const auto it = std::next(std::cbegin(qualities), pos);
    return std::count_if(it, std::next(it, size),
                         [min_quality] (const auto quality) {
                             return quality < min_quality;
                         });
}

} // namespace

void CigarScanner::do_add_read(const AlignedRead& read)
{
    using std::cbegin; using std::next; using std::move;
    
    const auto& read_contig   = contig_name(read);
    const auto& read_sequence = read.sequence();
    
    auto sequence_itr  = cbegin(read_sequence);
    auto qualities_itr = cbegin(read.qualities());
    
    auto ref_index = mapped_begin(read);
    std::size_t read_index {0};
    GenomicRegion region;
    
    for (const auto& cigar_operation : read.cigar()) {
        const auto op_size = cigar_operation.size();
        
        using Flag = CigarOperation::Flag;
        
        switch (cigar_operation.flag()) {
            case Flag::AlignmentMatch:
                add_snvs_in_match_range(GenomicRegion {read_contig, ref_index, ref_index + op_size},
                                        next(sequence_itr, read_index),
                                        next(sequence_itr, read_index + op_size),
                                        next(qualities_itr, read_index));
                
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            case Flag::SequenceMatch:
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            case Flag::Substitution:
            {
                region = GenomicRegion {read_contig, ref_index, ref_index + op_size};
                
                if (op_size <= options_.max_variant_size) {
                    auto removed_sequence = reference_.get().fetch_sequence(region);
                    auto added_sequence   = splice(read_sequence, read_index, op_size);
                    
                    if (is_good_sequence(removed_sequence) && is_good_sequence(added_sequence)) {
                        add_candidate(move(region), move(removed_sequence), move(added_sequence));
                    }
                }
                
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            }
            case Flag::Insertion:
            {
                if (count_bad_qualities(read.qualities(), read_index, op_size, options_.min_base_quality)
                        <= options_.max_poor_quality_insertion_bases) {
                    auto added_sequence = splice(read_sequence, read_index, op_size);
                    
                    if (is_good_sequence(added_sequence)) {
                        add_candidate(GenomicRegion {read_contig, ref_index, ref_index},
                                      "", move(added_sequence));
                    }
                }
                
                read_index += op_size;
                
                break;
            }
            case Flag::Deletion:
            {
                if (*next(qualities_itr, read_index) > 0 && *next(qualities_itr, read_index + 1) > 0) {
                    region = GenomicRegion {read_contig, ref_index, ref_index + op_size};
                    
                    auto removed_sequence = reference_.get().fetch_sequence(region);
                    
                    add_candidate(move(region), move(removed_sequence), "");
                }
                
                ref_index += op_size;
                
                break;
            }
            case Flag::SoftClipped:
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            case Flag::HardClipped:
                break;
            case Flag::Padding:
            case Flag::Skipped:
                ref_index += op_size;
                
                break;
        }
    }
}

void CigarScanner::do_add_reads(VectorIterator first, VectorIterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { do_add_read(read); });
}

void CigarScanner::do_add_reads(FlatSetIterator first, FlatSetIterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { do_add_read(read); });
}

template <typename ForwardIt, typename Container>
auto copy_overlapped_indels(ForwardIt first, const ForwardIt last, const Container& candidates)
{
    std::deque<Variant> result {};
    
    if (first == last || candidates.empty()) return result;
    
    const auto max_indel_size = size(largest_region(first, last));
    
    for (const auto& candidate : candidates) {
        const auto overlapped = overlap_range(first, last, candidate, max_indel_size);
        
        if (!overlapped.empty()) {
            std::copy_if(std::cbegin(overlapped), std::cend(overlapped), std::back_inserter(result),
                         [] (const auto& v) { return is_indel(v); });
            
            first = std::cbegin(overlapped).base();
            
            if (first == last) break;
        }
    }
    
    return result;
}

template <typename Container1, typename Container2>
auto copy_overlapped_indels(const Container1& all_candidates, const Container2& selected_candidates)
{
    using std::cbegin; using std::cend; using std::begin; using std::end;
    
    std::vector<Variant> indels {};
    
    static const auto is_indel = [] (const auto& v) { return octopus::is_indel(v); };
    
    indels.reserve(std::count_if(cbegin(all_candidates), cend(all_candidates), is_indel));
    
    std::copy_if(cbegin(all_candidates), cend(all_candidates), std::back_inserter(indels), is_indel);
    
    indels.erase(std::unique(begin(indels), end(indels)), end(indels));
    
    return copy_overlapped_indels(cbegin(indels), cend(indels), selected_candidates);
}

std::vector<Variant> CigarScanner::do_generate_variants(const GenomicRegion& region)
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::distance;
    
    std::sort(begin(candidates_), end(candidates_));
    
    auto overlapped = overlap_range(candidates_, region, max_seen_candidate_size_);
    
    std::vector<Variant> result {};
    
    if (options_.min_support < 2) {
        result.insert(end(result), cbegin(overlapped), cend(overlapped));
        result.erase(std::unique(begin(result), end(result)), end(result));
    } else {
        result.reserve(size(overlapped, BidirectionallySortedTag {})); // the maximum
        
        while (true) {
            // TODO: Clang 3.8 has a bug which wont let me use the std::function match_ member
            // here.
            const auto TEMP_FIX = [] (const Variant& lhs, const Variant& rhs) -> bool {
                return (is_insertion(lhs) && is_same_region(lhs, rhs)) || lhs == rhs;
            };
            
            const auto it = std::adjacent_find(cbegin(overlapped), cend(overlapped), TEMP_FIX);
            
            if (it == cend(overlapped)) break;
            
            const Variant& duplicate {*it};
            
            const auto it2 = std::find_if_not(std::next(it), cend(overlapped),
                                              [this, &duplicate] (const auto& variant) {
                                                  return match_(variant, duplicate);
                                              });
            
            const auto duplicate_count = distance(it, it2);
            
            if (duplicate_count >= options_.min_support) {
                if (duplicate == *std::prev(it2)) {
                    result.push_back(duplicate);
                } else {
                    std::unique_copy(it, it2, std::back_inserter(result));
                }
            }
            
            if (it2 == cend(overlapped)) break;
            
            overlapped.advance_begin(distance(cbegin(overlapped), it) + duplicate_count);
        }
        
        if (options_.always_include_overlapping_indels) {
            auto overlapped_indels = copy_overlapped_indels(candidates_, result);
            
            std::sort(begin(overlapped_indels), end(overlapped_indels));
            
            result.reserve(result.size() + overlapped_indels.size());
            
            const auto it = end(result);
            
            std::unique_copy(begin(overlapped_indels), end(overlapped_indels),
                             std::back_inserter(result));
            
            std::inplace_merge(begin(result), it, end(result));
            
            result.erase(std::unique(begin(result), end(result)), end(result));
        }
        
        result.shrink_to_fit();
    }
    
    return result;
}

void CigarScanner::do_clear() noexcept
{
    candidates_.clear();
    candidates_.shrink_to_fit();
}

std::string CigarScanner::name() const
{
    return "raw CIGAR";
}

// private methods

void CigarScanner::add_snvs_in_match_range(const GenomicRegion& region, const SequenceIterator first_base,
                                          const SequenceIterator last_base, const QualitiesIterator first_quality)
{
    using boost::make_zip_iterator; using std::for_each; using std::cbegin; using std::cend;
    
    using Tuple = boost::tuple<char, char, AlignedRead::BaseQuality>;
    
    const NucleotideSequence ref_segment {reference_.get().fetch_sequence(region)};
    
    const auto& contig = region.contig_name();
    
    const auto last_quality = std::next(first_quality, std::distance(first_base, last_base));
    
    auto ref_index = mapped_begin(region);
    
    for_each(make_zip_iterator(boost::make_tuple(cbegin(ref_segment), first_base, first_quality)),
             make_zip_iterator(boost::make_tuple(cend(ref_segment), last_base, last_quality)),
             [this, &contig, &ref_index] (const Tuple& t) {
                 const char ref_base  {t.get<0>()}, read_base {t.get<1>()};
                 
                 if (ref_base != read_base && ref_base != 'N' && read_base != 'N'
                     && t.get<2>() >= options_.min_base_quality) {
                     add_candidate(GenomicRegion {contig, ref_index, ref_index + 1}, ref_base, read_base);
                 }
                 
                 ++ref_index;
             });
}

} // namespace coretools
} // namespace octopus
