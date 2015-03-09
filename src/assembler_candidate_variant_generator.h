//
//  assembler_candidate_variant_generator.h
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__assembler_candidate_variant_generator__
#define __Octopus__assembler_candidate_variant_generator__

#include <vector>
#include <cstddef> // std::size_t

#include "i_variant_candidate_generator.h"
#include "variant_assembler.h"

class ReferenceGenome;
class AlignedRead;
class GenomicRegion;
class Variant;

class AssemblerCandidateVariantGenerator : public IVariantCandidateGenerator
{
public:
    AssemblerCandidateVariantGenerator() = delete;
    explicit AssemblerCandidateVariantGenerator(ReferenceGenome& the_reference, unsigned kmer_size);
    ~AssemblerCandidateVariantGenerator() override = default;
    
    AssemblerCandidateVariantGenerator(const AssemblerCandidateVariantGenerator&)            = default;
    AssemblerCandidateVariantGenerator& operator=(const AssemblerCandidateVariantGenerator&) = default;
    AssemblerCandidateVariantGenerator(AssemblerCandidateVariantGenerator&&)                 = default;
    AssemblerCandidateVariantGenerator& operator=(AssemblerCandidateVariantGenerator&&)      = default;
    
    void add_read(const AlignedRead& a_read) override;
    std::vector<Variant> get_candidates(const GenomicRegion& a_region) override;
    void reserve(std::size_t n) override;
    void clear() override;
    
private:
    ReferenceGenome& the_reference_;
    VariantAssembler the_variant_assembler_;
};

#endif /* defined(__Octopus__assembler_candidate_variant_generator__) */
