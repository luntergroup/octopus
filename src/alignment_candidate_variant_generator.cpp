//
//  alignment_candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "alignment_candidate_variant_generator.hpp"

#include <iterator>  // std::cbegin etc, std::distance
#include <algorithm> // std::for_each, std::sort, std::unique, std::adjacent_find, std::find_if_not
#include <boost/range/combine.hpp>

#include "reference_genome.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"
#include "cigar_string.hpp"
#include "mappable_algorithms.hpp"

#include <iostream> // DEBUG

namespace Octopus {

template <typename SequenceType>
bool is_good_sequence(const SequenceType& sequence) noexcept
{
    return std::none_of(std::cbegin(sequence), std::cend(sequence), [] (auto base) { return base == 'N'; });
}

// public methods

AlignmentCandidateVariantGenerator::AlignmentCandidateVariantGenerator(ReferenceGenome& reference,
                                                                       QualityType min_base_quality,
                                                                       unsigned min_supporting_reads,
                                                                       SizeType max_variant_size)
:
reference_ {reference},
min_base_quality_ {min_base_quality},
min_supporting_reads_ {min_supporting_reads},
max_variant_size_ {max_variant_size},
candidates_ {},
are_candidates_sorted_ {true},
max_seen_candidate_size_ {}
{}

void AlignmentCandidateVariantGenerator::add_read(const AlignedRead& read)
{
    const auto& contig_name   = get_contig_name(read);
    const auto& read_sequence = read.get_sequence();
    
    auto sequence_itr  = std::cbegin(read_sequence);
    auto qualities_itr = std::cbegin(read.get_qualities());
    
    auto ref_index = get_begin(read);
    AlignedRead::SizeType read_index {};
    GenomicRegion region {};
    
    for (const auto& cigar_operation : read.get_cigar_string()) {
        auto op_size = cigar_operation.get_size();
        
        switch (cigar_operation.get_flag()) {
            case CigarOperation::ALIGNMENT_MATCH:
                add_snvs_in_match_range(GenomicRegion {contig_name, ref_index, ref_index + op_size},
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
                region = GenomicRegion {contig_name, ref_index, ref_index + op_size};
                
                auto removed_sequence = reference_.get_sequence(region);
                auto added_sequence   = read_sequence.substr(read_index, op_size);
                
                if (is_good_sequence(removed_sequence) && is_good_sequence(added_sequence)) {
                    add_candidate(region, std::move(removed_sequence), std::move(added_sequence));
                }
                
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            }
            case CigarOperation::INSERTION:
            {
                auto added_sequence = read_sequence.substr(read_index, op_size);
                
                if (is_good_sequence(added_sequence)) {
                    add_candidate(GenomicRegion {contig_name, ref_index, ref_index},
                                  "", std::move(added_sequence));
                }
                
                read_index += op_size;
                
                break;
            }
            case CigarOperation::DELETION:
            {
                region = GenomicRegion {contig_name, ref_index, ref_index + op_size};
                
                auto removed_sequence = reference_.get_sequence(region);
                
                if (is_good_sequence(removed_sequence)) {
                    add_candidate(region, std::move(removed_sequence), "");
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

void AlignmentCandidateVariantGenerator::add_reads(MappableSet<AlignedRead>::const_iterator first,
                                                   MappableSet<AlignedRead>::const_iterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { add_read(read); });
    candidates_.shrink_to_fit();
}

std::vector<Variant> AlignmentCandidateVariantGenerator::get_candidates(const GenomicRegion& region)
{
    using std::begin; using std::end; using std::cbegin; using std::cend;
    
    if (!are_candidates_sorted_) {
        std::sort(begin(candidates_), end(candidates_));
        
        if (min_supporting_reads_ == 1) {
            candidates_.erase(std::unique(begin(candidates_), end(candidates_)), end(candidates_));
        }
        
        are_candidates_sorted_ = true;
    }
    
    auto overlapped = overlap_range(cbegin(candidates_), cend(candidates_), region, max_seen_candidate_size_);
    
    if (min_supporting_reads_ == 1) {
        return std::vector<Variant> {begin(overlapped), end(overlapped)};
    } else {
        std::vector<Variant> result {};
        result.reserve(bases(overlapped).size());
        
        while (!overlapped.empty()) {
            auto it = std::adjacent_find(overlapped.begin(), overlapped.end());
            
            if (it == overlapped.end()) break;
            
            const Variant& duplicate = *it;
            
            auto it2 = std::find_if_not(std::next(it), overlapped.end(),
                                        [&duplicate] (const auto& variant) { return variant == duplicate; });
            
            auto duplicate_count = std::distance(it, it2);
            
            if (duplicate_count >= min_supporting_reads_) {
                result.emplace_back(duplicate);
            }
            
            overlapped.advance_begin(duplicate_count);
        }
        
        result.shrink_to_fit();
        
        return result;
    }
}

void AlignmentCandidateVariantGenerator::clear()
{
    candidates_.clear();
}

// private methods    

void AlignmentCandidateVariantGenerator::
add_snvs_in_match_range(const GenomicRegion& region, SequenceIterator first_base,
                        SequenceIterator last_base, QualitiesIterator first_quality)
{
    using boost::make_zip_iterator;
    
    using Tuple = boost::tuple<SequenceType::value_type, SequenceType::value_type, QualityType>;
    
    auto ref_segment = reference_.get_sequence(region);
    auto ref_index   = get_begin(region);
    
    std::for_each(make_zip_iterator(boost::make_tuple(std::cbegin(ref_segment), first_base, first_quality)),
                  make_zip_iterator(boost::make_tuple(std::cend(ref_segment), last_base, first_quality
                                                   + std::distance(first_base, last_base))),
        [this, &region, &ref_index] (const Tuple& t) {
            auto ref_base  = t.get<0>();
            auto read_base = t.get<1>();
            
            if (ref_base != read_base && ref_base != 'N' && read_base != 'N' && t.get<2>() >= min_base_quality_) {
                add_candidate(GenomicRegion {region.get_contig_name(), ref_index, ref_index + 1},
                              ref_base, read_base);
            }
            
            ++ref_index;
    });
}

} // namespace Octopus
