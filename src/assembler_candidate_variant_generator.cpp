//
//  assembler_candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "assembler_candidate_variant_generator.hpp"

#include <algorithm>
#include <iterator>
#include <cassert>

#include "common.hpp"
#include "reference_genome.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"
#include "sequence_utils.hpp"
#include "logging.hpp"

namespace Octopus
{
AssemblerCandidateVariantGenerator::AssemblerCandidateVariantGenerator(const ReferenceGenome& reference,
                                                                       unsigned kmer_size,
                                                                       QualityType min_base_quality,
                                                                       unsigned min_supporting_reads,
                                                                       SizeType max_variant_size)
:
reference_ {reference},
initial_kmer_sizes_ {kmer_size},
fallback_kmer_sizes_ {},
assembler_ {initial_kmer_sizes_.front()},
region_assembled_ {},
min_base_quality_ {min_base_quality},
min_supporting_reads_ {min_supporting_reads},
max_variant_size_ {max_variant_size}
{}

bool AssemblerCandidateVariantGenerator::requires_reads() const noexcept
{
    return true;
}

bool are_all_bases_good_quality(const AlignedRead& read,
                                const AlignedRead::QualityType min_quality)
{
    return std::all_of(std::cbegin(read.get_qualities()),
                       std::cend(read.get_qualities()),
                       [min_quality] (const char quality) {
                           return quality >= min_quality;
                       });
}

AlignedRead::SequenceType
transform_low_base_qualities_to_n(const AlignedRead& read, const AlignedRead::QualityType min_quality)
{
    auto result = read.get_sequence();
    
    std::transform(std::cbegin(result), std::cend(result),
                   std::cbegin(read.get_qualities()),
                   std::begin(result),
                   [min_quality] (const char base, const auto quality) {
                       return (quality >= min_quality) ? base : 'N';
                   });
    
    return result;
}

void AssemblerCandidateVariantGenerator::add_read(const AlignedRead& read)
{
    if (are_all_bases_good_quality(read, min_base_quality_)) {
        assembler_.insert_read(read.get_sequence());
    } else {
        auto masked_sequence = transform_low_base_qualities_to_n(read, min_base_quality_);
         assembler_.insert_read(masked_sequence);
    }
    
    if (region_assembled_) {
        region_assembled_ = encompassing_region(read, *region_assembled_);
    } else {
        region_assembled_ = mapped_region(read);
    }
}

void AssemblerCandidateVariantGenerator::add_reads(std::vector<AlignedRead>::const_iterator first,
                                                   std::vector<AlignedRead>::const_iterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { add_read(read); });
}

void AssemblerCandidateVariantGenerator::add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                                                   MappableFlatMultiSet<AlignedRead>::const_iterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { add_read(read); });
}

void trim_variant(Assembler::Variant& v)
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
void trim_variants(Container& variants)
{
    for (auto& v : variants) {
        trim_variant(v);
    }
}

template <typename C1, typename C2>
void add_to_mapped_variants(C1& result, C2&& variants, const GenomicRegion& region)
{
    result.reserve(result.size() + variants.size());
    
    const auto it = std::end(result);
    
    for (auto& variant : variants) {
        result.emplace_back(contig_name(region),
                            region.get_begin() + static_cast<GenomicRegion::SizeType>(variant.begin_pos),
                            std::move(variant.ref),
                            std::move(variant.alt)
                            );
    }
    
    std::sort(it, std::end(result));
    
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
}

template <typename C, typename S>
void assert_all_reference_consistent(const C& variants, const S& ref_sequence)
{
    assert(std::all_of(std::cbegin(variants), std::cend(variants),
                       [&ref_sequence] (const auto& v) {
                           return v.begin_pos + v.ref.size() < ref_sequence.size()
                                    && std::equal(std::cbegin(v.ref), std::cend(v.ref),
                                                  std::next(std::begin(ref_sequence),
                                                            v.begin_pos));
                       }));
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

std::vector<Variant> AssemblerCandidateVariantGenerator::generate_candidates(const GenomicRegion& region)
{
    std::vector<Variant> result {};
    
    if (!region_assembled_) {
        return result;
    }
    
    const auto kmer_size = initial_kmer_sizes_.front();
    
    const auto reference_region = expand(*region_assembled_, kmer_size);
    
    const auto ref_sequence = reference_.get().get_sequence(reference_region);
    
    if (has_ns(ref_sequence)) {
        return result;
    }
    
    //std::cout << "Reference region = " << reference_region << std::endl;
    //std::cout << "Reference = " << ref_sequence << std::endl;
    
    assembler_.insert_reference(ref_sequence);
    
    if (assembler_.is_all_reference()) {
        return result;
    }
    
    assembler_.remove_trivial_nonreference_cycles();
    
    if (!assembler_.prune(min_supporting_reads_)) {
        Logging::WarningLogger log {};
        log << "Assembler could not generate candidates due to cyclic graph";
        return result;
    }
    
//    std::cout << "Final graph:" << std::endl;
//    debug::print(assembler_);
    
    auto variants = assembler_.extract_variants();
    
    assembler_.clear();
    
    trim_variants(variants);
    
//    add_to_mapped_variants(result, std::move(variants), reference_region);
//    debug::print_generated_candidates(result, "local re-assembly");
//    exit(0);
    
    assert_all_reference_consistent(variants, ref_sequence);
    
    add_to_mapped_variants(result, std::move(variants), reference_region);
    
    remove_nonoverlapping(result, region); // as we expanded original region
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        debug::print_generated_candidates(stream(log), result, "local re-assembly");
    }
    
    return result;
}

void AssemblerCandidateVariantGenerator::clear()
{
    assembler_.clear();
}

} // namespace Octopus
