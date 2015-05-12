//
//  variant_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 29/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_caller.h"
#include "bayesian_genotype_model.h"

#include <algorithm>

using BayesianGenotypeModel::probability_allele_in_sample;

VariantCalls<HaplotypePhaser::SampleIdType, HaplotypePhaser::RealType>
call_variants(const std::vector<HaplotypePhaser::SampleIdType>& the_samples,
              const HaplotypePhaser::PhasedRegions& the_phased_regions,
              const std::vector<Variant>& the_candidates, const VariationalBayesGenotypeModel& the_model)
{
    VariantCalls<HaplotypePhaser::SampleIdType, HaplotypePhaser::RealType> result {};
    
    for (const auto& candidate : the_candidates) {
        Allele reference_allele {candidate.get_reference_allele()};
        Allele alternative_allele {candidate.get_alternative_allele()};
        
        auto it = std::find_if(std::cbegin(the_phased_regions), std::cend(the_phased_regions),
                               [&alternative_allele] (const auto& phased_region) {
                                   return overlaps(alternative_allele, phased_region.the_region);
                               });
        for (const auto& sample : the_samples ) {
            auto p = probability_allele_in_sample(alternative_allele, it->the_haplotypes,
                                                  it->the_latent_posteriors.genotype_posteriors.at(sample),
                                                  it->the_genotypes);
            
            if (p > 0.9) {
                result[sample][candidate] = p;
            }
        }
    }
    
    return result;
}
