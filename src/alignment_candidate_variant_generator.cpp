//
//  alignment_candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "alignment_candidate_variant_generator.hpp"

#include <iterator>
#include <algorithm>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "common.hpp"
#include "reference_genome.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"
#include "cigar_string.hpp"
#include "mappable_ranges.hpp"
#include "mappable_algorithms.hpp"
#include "logging.hpp"

#include <iostream> // DEBUG

namespace Octopus
{

// public methods

AlignmentCandidateVariantGenerator::AlignmentCandidateVariantGenerator(const ReferenceGenome& reference,
                                                                       QualityType min_base_quality,
                                                                       unsigned min_support,
                                                                       SizeType max_variant_size)
:
reference_ {reference},
min_base_quality_ {min_base_quality},
max_poor_quality_insertion_bases_ {1},
min_support_ {min_support},
max_variant_size_ {max_variant_size},
candidates_ {},
max_seen_candidate_size_ {}
{}

bool AlignmentCandidateVariantGenerator::requires_reads() const noexcept
{
    return true;
}

namespace
{
    template <typename SequenceType>
    bool is_good_sequence(const SequenceType& sequence) noexcept
    {
        return std::none_of(std::cbegin(sequence), std::cend(sequence),
                            [] (auto base) { return base == 'N'; });
    }
    
    template <typename SequenceType, typename T>
    SequenceType splice(const SequenceType& sequence, const T pos, const T size)
    {
        const auto it = std::next(std::cbegin(sequence), pos);
        return SequenceType {it, std::next(it, size)};
    }
    
    template <typename Q, typename T>
    bool is_all_good_quality(const Q& qualities, const T pos, const T size,
                             const typename Q::value_type min_quality)
    {
        const auto it = std::next(std::cbegin(qualities), pos);
        return std::all_of(it, std::next(it, size),
                           [min_quality] (const auto quality) {
                               return quality >= min_quality;
                           });
    }
    
    template <typename Q, typename T>
    auto count_bad_qualities(const Q& qualities, const T pos, const T size,
                             const typename Q::value_type min_quality)
    {
        const auto it = std::next(std::cbegin(qualities), pos);
        return std::count_if(it, std::next(it, size),
                             [min_quality] (const auto quality) {
                                 return quality < min_quality;
                           });
    }
} // namespace

void AlignmentCandidateVariantGenerator::add_read(const AlignedRead& read)
{
    const auto& read_contig   = contig_name(read);
    const auto& read_sequence = read.get_sequence();
    
    auto sequence_itr  = std::cbegin(read_sequence);
    auto qualities_itr = std::cbegin(read.get_qualities());
    
    auto ref_index = region_begin(read);
    AlignedRead::SizeType read_index {0};
    GenomicRegion region;
    
    for (const auto& cigar_operation : read.get_cigar_string()) {
        const auto op_size = cigar_operation.get_size();
        
        switch (cigar_operation.get_flag()) {
            case CigarOperation::ALIGNMENT_MATCH:
                add_snvs_in_match_range(GenomicRegion {read_contig, ref_index, ref_index + op_size},
                                        sequence_itr + read_index, sequence_itr + read_index + op_size,
                                        qualities_itr + read_index);
                
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            case CigarOperation::SEQUENCE_MATCH:
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            case CigarOperation::SUBSTITUTION:
            {
                region = GenomicRegion {read_contig, ref_index, ref_index + op_size};
                
                if (op_size <= max_variant_size_) {
                    auto removed_sequence = reference_.get().get_sequence(region);
                    auto added_sequence   = splice(read_sequence, read_index, op_size);
                    
                    if (is_good_sequence(removed_sequence) && is_good_sequence(added_sequence)) {
                        add_candidate(std::move(region), std::move(removed_sequence),
                                      std::move(added_sequence));
                    }
                }
                
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            }
            case CigarOperation::INSERTION:
            {
                if (count_bad_qualities(read.get_qualities(), read_index, op_size, min_base_quality_)
                        <= max_poor_quality_insertion_bases_) {
                    auto added_sequence = splice(read_sequence, read_index, op_size);
                    
                    if (is_good_sequence(added_sequence)) {
                        add_candidate(GenomicRegion {read_contig, ref_index, ref_index},
                                      "", std::move(added_sequence));
                    }
                }
                
                read_index += op_size;
                
                break;
            }
            case CigarOperation::DELETION:
            {
                region = GenomicRegion {read_contig, ref_index, ref_index + op_size};
                
                if (count_bad_qualities(read.get_qualities(), read_index, op_size,
                                        min_base_quality_) < op_size) {
                    auto removed_sequence = reference_.get().get_sequence(region);
                    add_candidate(std::move(region), std::move(removed_sequence), "");
                }
                
                ref_index += op_size;
                
                break;
            }
            case CigarOperation::SKIPPED:
                ref_index += op_size;
                
                break;
            case CigarOperation::SOFT_CLIPPED:
                read_index += op_size;
                ref_index  += op_size;
                
                break;
        }
    }
}

void AlignmentCandidateVariantGenerator::add_reads(std::vector<AlignedRead>::const_iterator first,
                                                   std::vector<AlignedRead>::const_iterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { add_read(read); });
    candidates_.shrink_to_fit();
}

void AlignmentCandidateVariantGenerator::add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                                                   MappableFlatMultiSet<AlignedRead>::const_iterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { add_read(read); });
    candidates_.shrink_to_fit();
}

std::vector<Variant> AlignmentCandidateVariantGenerator::generate_candidates(const GenomicRegion& region)
{
    using std::begin; using std::end; using std::distance;
    
    std::sort(begin(candidates_), end(candidates_));
    
    auto overlapped = overlap_range(candidates_, region, max_seen_candidate_size_);
    
    std::vector<Variant> result {};
    
    if (min_support_ < 2) {
        result.insert(end(result), begin(overlapped), end(overlapped));
        result.erase(std::unique(begin(result), end(result)), end(result));
    } else {
        result.reserve(size(overlapped, BidirectionallySortedTag {})); // the maximum
        
        while (true) {
            const auto it = std::adjacent_find(begin(overlapped), end(overlapped));
            
            if (it == end(overlapped)) break;
            
            const Variant& duplicate {*it};
            
            const auto it2 = std::find_if_not(std::next(it), end(overlapped),
                                              [&] (const auto& variant) {
                                                  return variant == duplicate;
                                              });
            
            const auto duplicate_count = distance(it, it2);
            
            if (duplicate_count >= min_support_) {
                result.emplace_back(duplicate);
            }
            
            if (it2 == end(overlapped)) break;
            
            overlapped.advance_begin(distance(begin(overlapped), it) + duplicate_count);
        }
        
        result.shrink_to_fit();
    }
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        debug::print_generated_candidates(stream(log), result, "raw CIGAR strings");
    }
    
    return result;
}

void AlignmentCandidateVariantGenerator::clear()
{
    candidates_.clear();
}

// private methods

void AlignmentCandidateVariantGenerator::
add_snvs_in_match_range(const GenomicRegion& region, const SequenceIterator first_base,
                        const SequenceIterator last_base, const QualitiesIterator first_quality)
{
    using boost::make_zip_iterator; using std::for_each; using std::cbegin; using std::cend;
    
    using Tuple = boost::tuple<char, char, QualityType>;
    
    const SequenceType ref_segment {reference_.get().get_sequence(region)};
    
    const auto& contig = region.get_contig_name();
    
    const auto last_quality = std::next(first_quality, std::distance(first_base, last_base));
    
    auto ref_index = region_begin(region);
    
    for_each(make_zip_iterator(boost::make_tuple(cbegin(ref_segment), first_base, first_quality)),
             make_zip_iterator(boost::make_tuple(cend(ref_segment), last_base, last_quality)),
             [this, &contig, &ref_index] (const Tuple& t) {
                 const char ref_base  {t.get<0>()}, read_base {t.get<1>()};
                 
                 if (ref_base != read_base && ref_base != 'N' && read_base != 'N'
                     && t.get<2>() >= min_base_quality_) {
                     add_candidate(GenomicRegion {contig, ref_index, ref_index + 1}, ref_base, read_base);
                 }
                 
                 ++ref_index;
             });
}
} // namespace Octopus
