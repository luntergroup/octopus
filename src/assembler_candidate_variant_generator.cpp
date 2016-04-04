//
//  assembler_candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "assembler_candidate_variant_generator.hpp"

#include <algorithm>

#include "common.hpp"
#include "reference_genome.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"
#include "logging.hpp"

namespace Octopus
{
AssemblerCandidateVariantGenerator::AssemblerCandidateVariantGenerator(const ReferenceGenome& reference,
                                                                       unsigned kmer_size,
                                                                       SizeType max_variant_size)
:
reference_ {reference},
kmer_size_ {31},
assembler_ {kmer_size_},
region_assembled_ {},
max_variant_size_ {max_variant_size}
{}

bool AssemblerCandidateVariantGenerator::requires_reads() const noexcept
{
    return true;
}

void AssemblerCandidateVariantGenerator::add_read(const AlignedRead& read)
{
    assembler_.insert_read(read.get_sequence());
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

std::vector<Variant> AssemblerCandidateVariantGenerator::generate_candidates(const GenomicRegion& region)
{
    std::vector<Variant> result {};
    
    if (!region_assembled_) {
        return result;
    }
    
    const auto reference_region = expand(*region_assembled_, kmer_size_);
    
    const auto ref_sequence = reference_.get().get_sequence(region);
    
    std::cout << "Reference = " << ref_sequence << std::endl;
    
    assembler_.insert_reference(ref_sequence);
    
    if (assembler_.is_all_reference()) {
        return result;
    }
    
    assembler_.remove_trivial_nonreference_cycles();
    assembler_.prune(2);
    
    std::cout << "Final graph:" << std::endl;
    debug::print_edges(assembler_);
    
    if (!assembler_.is_acyclic()) {
        Logging::WarningLogger log {};
        log << "Assembler could not generate candidates due to cyclic graph";
        return result;
    }
    
    std::cout << "no cycles. Extracting variants..." << std::endl;
    
    auto variants = assembler_.extract_variants();
    
    for (auto v : variants) {
        std::cout << v.begin_pos << " " << v.ref << " " << v.alt << std::endl;
    }
    
    return result;
}

void AssemblerCandidateVariantGenerator::clear()
{
    assembler_.clear();
}

} // namespace Octopus
