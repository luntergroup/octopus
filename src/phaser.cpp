//
//  phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 20/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "phaser.hpp"

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
        return PhaseSet {haplotypes.front().get_region()};
    }
    
} // namespace Ocotpus
