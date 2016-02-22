//
//  phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 20/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "phaser.hpp"

#include <cassert>

namespace Octopus
{
    Phaser::Phaser(const double min_phase_score)
    :
    min_phase_score_ {min_phase_score}
    {}
    
    boost::optional<Phaser::PhaseSet>
    Phaser::try_phase(const std::vector<Haplotype>& haplotypes,
                      const GenotypePosteriorMap& genotype_posteriors)
    {
        return boost::none;
    }
    
    Phaser::PhaseSet
    Phaser::force_phase(const std::vector<Haplotype>& haplotypes,
                        const GenotypePosteriorMap& genotype_posteriors)
    {
        assert(!haplotypes.empty());
        assert(!genotype_posteriors.empty1());
        assert(!genotype_posteriors.empty2());
        
        const auto& haplotype_region = haplotypes.front().get_region();
        
        PhaseSet result {haplotype_region};
        
        result.phase_regions.reserve(genotype_posteriors.size1());
        
        for (const auto& p : genotype_posteriors) {
            result.phase_regions[p.first].emplace_back(haplotype_region, 1);
        }
        
        return result;
    }
    
} // namespace Ocotpus
