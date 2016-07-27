//
//  candidate_variant_generator.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__candidate_variant_generator__
#define __Octopus__candidate_variant_generator__

#include <vector>
#include <memory>
#include <cstddef>

#include "i_candidate_variant_generator.hpp"

class GenomicRegion;
class AlignedRead;
class Variant;

namespace Octopus {

class CandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    CandidateVariantGenerator() = default;
    
    CandidateVariantGenerator(const CandidateVariantGenerator&)            = delete;
    CandidateVariantGenerator& operator=(const CandidateVariantGenerator&) = delete;
    CandidateVariantGenerator(CandidateVariantGenerator&&)                 = default;
    CandidateVariantGenerator& operator=(CandidateVariantGenerator&&)      = default;
    
    ~CandidateVariantGenerator() override = default;
    
    void register_generator(std::unique_ptr<ICandidateVariantGenerator> generator);
    
    bool requires_reads() const noexcept override;
    
    void add_read(const AlignedRead& read) override;
    void add_reads(std::vector<AlignedRead>::const_iterator first,
                   std::vector<AlignedRead>::const_iterator last) override;
    void add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                   MappableFlatMultiSet<AlignedRead>::const_iterator last) override;
    
    std::vector<Variant> generate_candidates(const GenomicRegion& region) override;
    
    void reserve(std::size_t n) override;
    void clear() override;
    
private:
    std::vector<std::unique_ptr<ICandidateVariantGenerator>> generators_;
};
} // namespace Octopus

#endif /* defined(__Octopus__variant_candidate_generator__) */
