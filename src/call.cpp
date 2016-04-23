//
//  call.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "call.hpp"

#include <algorithm>

namespace Octopus
{
    double Call::get_quality() const noexcept
    {
        return quality_;
    }
    
    const Call::GenotypeCall& Call::get_genotype_call(const SampleIdType& sample) const
    {
        return genotype_calls_.at(sample);
    }
    
    bool Call::is_phased(const SampleIdType& sample) const
    {
        return static_cast<bool>(genotype_calls_.at(sample).phase);
    }
    
    bool Call::all_phased() const noexcept
    {
        return std::all_of(std::cbegin(genotype_calls_), std::cend(genotype_calls_),
                           [] (const auto& p) {
                               return static_cast<bool>(p.second.phase);
                           });
    }
    
    void Call::set_phase(const SampleIdType& sample, PhaseCall phase)
    {
        genotype_calls_.at(sample).phase = std::move(phase);
    }
} // namespace Octopus
