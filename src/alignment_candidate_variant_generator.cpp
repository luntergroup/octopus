//
//  alignment_candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "alignment_candidate_variant_generator.h"

#include <algorithm> // std::for_each, std::sort, std::unique, std::lower_bound, std::upper_bound
#include <iterator>  // std::cbegin etc, std::distance
#include <boost/range/combine.hpp>

#include "reference_genome.h"
#include "aligned_read.h"
#include "variant.h"
#include "cigar_string.h"
#include "region_algorithms.h"

AlignmentCandidateVariantGenerator::AlignmentCandidateVariantGenerator(ReferenceGenome& the_reference,
                                                                       QualityType min_base_quality)
:
the_reference_ {the_reference},
candidates_ {},
min_base_quality_ {min_base_quality},
are_candidates_sorted_ {true}
{}

void AlignmentCandidateVariantGenerator::add_read(const AlignedRead& a_read)
{
    const auto& contig_name       = get_contig_name(a_read);
    const auto& the_read_sequence = a_read.get_sequence();
    const auto& the_qualities     = a_read.get_qualities();
    
    auto sequence_begin  = std::cbegin(the_read_sequence);
    auto qualities_begin = std::cbegin(the_qualities);
    auto ref_index       = get_begin(a_read);
    AlignedRead::SizeType read_index {};
    CigarOperation::SizeType op_size {};
    GenomicRegion a_region {};
    
    for (const auto& cigar_operation : a_read.get_cigar_string()) {
        op_size = cigar_operation.get_size();
        
        switch (cigar_operation.get_flag()) {
            case CigarOperation::ALIGNMENT_MATCH:
                get_variants_in_match_range(GenomicRegion {contig_name, ref_index, ref_index + op_size},
                                            sequence_begin + read_index,
                                            sequence_begin + read_index + op_size,
                                            qualities_begin + read_index);
                
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            case CigarOperation::SEQUENCE_MATCH:
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            case CigarOperation::SUBSTITUTION:
            {
                a_region = GenomicRegion {contig_name, ref_index, ref_index + op_size};
                auto removed_sequence = the_reference_.get_sequence(a_region);
                auto added_sequence   = the_read_sequence.substr(read_index, op_size);
                
                if (is_good_sequence(removed_sequence) && is_good_sequence(added_sequence)) {
                    add_variant(a_region, std::move(removed_sequence), std::move(added_sequence));
                }
                
                read_index += op_size;
                ref_index  += op_size;
                
                break;
            }
            case CigarOperation::INSERTION:
            {
                auto added_sequence = the_read_sequence.substr(read_index, op_size);
                
                if (is_good_sequence(added_sequence)) {
                    add_variant(GenomicRegion {contig_name, ref_index, ref_index},
                                "", std::move(added_sequence));
                }
                
                read_index += op_size;
                
                break;
            }
            case CigarOperation::DELETION:
            {
                a_region = GenomicRegion {contig_name, ref_index, ref_index + op_size};
                auto removed_sequence = the_reference_.get_sequence(a_region);
                
                if (is_good_sequence(removed_sequence)) {
                    add_variant(a_region, std::move(removed_sequence), "");
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

void AlignmentCandidateVariantGenerator::add_reads(ReadIterator first, ReadIterator last)
{
    candidates_.reserve(estimate_num_variants(std::distance(first, last)));
    std::for_each(first, last, [this] (const auto& a_read ) { add_read(a_read); });
    candidates_.shrink_to_fit();
}

std::vector<Variant> AlignmentCandidateVariantGenerator::get_candidates(const GenomicRegion& a_region)
{
    if (!are_candidates_sorted_) {
        std::sort(std::begin(candidates_), std::end(candidates_));
        auto it = std::unique(std::begin(candidates_), std::end(candidates_));
        candidates_.resize(std::distance(std::begin(candidates_), it));
        are_candidates_sorted_ = true;
    }
    
    auto range = overlap_range(std::cbegin(candidates_), std::cend(candidates_), a_region);
    
    return std::vector<Variant> {range.begin(), range.end()};
}

void AlignmentCandidateVariantGenerator::reserve(std::size_t n)
{
    candidates_.reserve(n);
}

void AlignmentCandidateVariantGenerator::clear()
{
    candidates_.clear();
}

void AlignmentCandidateVariantGenerator::
get_variants_in_match_range(const GenomicRegion& the_region, SequenceIterator first_base,
                            SequenceIterator last_base, QualitiesIterator first_quality)
{
    auto ref_segment = the_reference_.get_sequence(the_region);
    auto ref_index   = get_begin(the_region);
    
    SequenceType::value_type ref_base {'N'}, read_base {'N'};
    
    std::for_each(
        boost::make_zip_iterator(boost::make_tuple(std::cbegin(ref_segment), first_base, first_quality)),
        boost::make_zip_iterator(boost::make_tuple(std::cend(ref_segment), last_base, first_quality
                                                   + std::distance(first_base, last_base))),
        [this, &the_region, &ref_index, &ref_base, &read_base]
                  (const boost::tuple<SequenceType::value_type, SequenceType::value_type, QualityType>& p) {
            ref_base  = p.get<0>();
            read_base = p.get<1>();
            
            if (ref_base != read_base && ref_base != 'N' && read_base != 'N' && p.get<2>() >= min_base_quality_) {
                add_variant(GenomicRegion {the_region.get_contig_name(), ref_index, ref_index + 1},
                            ref_base, read_base);
            }
            
            ++ref_index;
    });
}

bool AlignmentCandidateVariantGenerator::is_good_sequence(const SequenceType& sequence) const noexcept
{
    return std::none_of(std::cbegin(sequence), std::cend(sequence), [] (auto base) {
        return base == 'N';
    });
}

std::size_t AlignmentCandidateVariantGenerator::estimate_num_variants(std::size_t num_reads) const noexcept
{
    return num_reads;
}
