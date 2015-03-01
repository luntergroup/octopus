//
//  prior_variant_candidates.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__prior_variant_candidates__
#define __Octopus__prior_variant_candidates__

#include "i_variant_candidate_generator.h"

class VariantFile;
class GenomicRegion;

class PriorVariantCandidates : public IVariantCandidateGenerator
{
public:
    PriorVariantCandidates() = delete;
    explicit PriorVariantCandidates(VariantFile& a_variant_source);
    ~PriorVariantCandidates() override = default;
    
    PriorVariantCandidates(const PriorVariantCandidates&)            = default;
    PriorVariantCandidates& operator=(const PriorVariantCandidates&) = default;
    PriorVariantCandidates(PriorVariantCandidates&&)                 = default;
    PriorVariantCandidates& operator=(PriorVariantCandidates&&)      = default;
    
    void add_read(const AlignedRead& a_read) override;
    std::set<Variant> get_candidates(const GenomicRegion& a_region) override;
    void clear() override;
    
private:
    VariantFile& a_variant_file_;
};

inline void PriorVariantCandidates::add_read(const AlignedRead& a_read) {}
inline void PriorVariantCandidates::clear() {}

#endif /* defined(__Octopus__prior_variant_candidates__) */
