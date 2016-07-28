//
//  reference_call.cpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "reference_call.hpp"

#include "mappable.hpp"

namespace octopus
{
    const GenomicRegion& ReferenceCall::mapped_region() const noexcept
    {
        return ::mapped_region(reference_);
    }
    
    const Allele& ReferenceCall::reference() const noexcept
    {
        return reference_;
    }
    
    void ReferenceCall::replace_called_alleles(const char old_base, const char replacement_base)
    {
        const auto& ref_sequence = reference_.sequence();
        
        auto it = std::find(std::cbegin(ref_sequence), std::cend(ref_sequence), old_base);
        
        if (it != std::cend(ref_sequence)) {
            Allele::NucleotideSequence new_sequence {};
            new_sequence.reserve(ref_sequence.size());
            
            new_sequence.insert(std::end(new_sequence), std::cbegin(ref_sequence), it);
            
            std::replace_copy(it, std::cend(ref_sequence), std::back_inserter(new_sequence),
                              old_base, replacement_base);
            
            reference_ = Allele {reference_.mapped_region(), std::move(new_sequence)};
        }
    }
    
    void ReferenceCall::replace(const Allele& old, Allele replacement)
    {
        if (reference_ == old) reference_ = std::move(replacement);
    }
    
    void ReferenceCall::replace_uncalled_genotype_alleles(const Allele& replacement, const char ignore)
    {
        
    }
    
    void ReferenceCall::decorate(VcfRecord::Builder& record) const
    {
        
    }
} // namespace octopus
