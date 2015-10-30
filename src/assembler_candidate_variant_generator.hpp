//
//  assembler_candidate_variant_generator.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__assembler_candidate_variant_generator__
#define __Octopus__assembler_candidate_variant_generator__

#include <vector>
#include <cstddef> // std::size_t

#include "i_candidate_variant_generator.hpp"
#include "variant_assembler.hpp"

class ReferenceGenome;
class AlignedRead;
class GenomicRegion;
class Variant;

namespace Octopus {

class AssemblerCandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    using SizeType = GenomicRegion::SizeType;
    
    AssemblerCandidateVariantGenerator() = delete;
    explicit AssemblerCandidateVariantGenerator(ReferenceGenome& reference, unsigned kmer_size,
                                                SizeType max_variant_size = 100);
    ~AssemblerCandidateVariantGenerator() override = default;
    
    AssemblerCandidateVariantGenerator(const AssemblerCandidateVariantGenerator&)            = default;
    AssemblerCandidateVariantGenerator& operator=(const AssemblerCandidateVariantGenerator&) = default;
    AssemblerCandidateVariantGenerator(AssemblerCandidateVariantGenerator&&)                 = default;
    AssemblerCandidateVariantGenerator& operator=(AssemblerCandidateVariantGenerator&&)      = default;
    
    void add_read(const AlignedRead& read) override;
    void add_reads(std::vector<AlignedRead>::const_iterator first, std::vector<AlignedRead>::const_iterator last) override;
    void add_reads(MappableSet<AlignedRead>::const_iterator first, MappableSet<AlignedRead>::const_iterator last) override;
    std::vector<Variant> get_candidates(const GenomicRegion& region) override;
    void clear() override;
    
private:
    ReferenceGenome& reference_;
    VariantAssembler assembler_;
    SizeType max_variant_size_;
};

} // namespace Octopus

#endif /* defined(__Octopus__assembler_candidate_variant_generator__) */
