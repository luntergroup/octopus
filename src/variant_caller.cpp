//
//  variant_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 29/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_caller.h"

#include <algorithm> // std::find_if

#include "bayesian_genotype_model.h"
#include "mappable_algorithms.h"

namespace Octopus
{

namespace VariantCaller
{
    using BayesianGenotypeModel::probability_allele_in_sample;

    HaplotypePhaser::PhasedRegions::const_iterator
    find_phased_region(const HaplotypePhaser::PhasedRegions& the_phased_regions, const Allele& an_allele)
    {
        return std::find_if(std::cbegin(the_phased_regions), std::cend(the_phased_regions),
                            [&an_allele] (const auto& phased_region) {
                                return overlaps(an_allele, phased_region.the_region);
                            });
    }
    
    AllelePosteriors<Octopus::SampleIdType, Octopus::ProbabilityType>
    get_allele_posteriors(const std::vector<Octopus::SampleIdType>& the_samples,
                          const HaplotypePhaser::PhasedRegions& the_phased_regions,
                          const std::vector<Allele>& the_alleles)
    {
        AllelePosteriors<Octopus::SampleIdType, Octopus::ProbabilityType> result {};
        
        for (const auto& allele : the_alleles) {
            auto it = find_phased_region(the_phased_regions, allele);
            
            for (const auto& sample : the_samples ) {
                result[sample][allele] = probability_allele_in_sample(allele, it->the_haplotypes,
                                                                      it->the_latent_posteriors.genotype_probabilities.at(sample),
                                                                      it->the_genotypes);
            }
        }
        
        return result;
    }
    
    AllelePosteriors<Octopus::SampleIdType, Octopus::ProbabilityType>
    get_allele_posteriors(const std::vector<Octopus::SampleIdType>& the_samples,
                          const HaplotypePhaser::PhasedRegions& the_phased_regions,
                          const std::vector<Variant>& the_candidates)
    {
        std::vector<Allele> the_alleles {};
        
        for (const Variant& a_variant : the_candidates) {
            the_alleles.emplace_back(a_variant.get_reference_allele());
            the_alleles.emplace_back(a_variant.get_alternative_allele());
        }
        
        return get_allele_posteriors(the_samples, the_phased_regions, the_alleles);
    }
    
    VariantCalls<Octopus::SampleIdType, Octopus::ProbabilityType>
    call_variants(const std::vector<Octopus::SampleIdType>& the_samples,
                  const HaplotypePhaser::PhasedRegions& the_phased_regions,
                  const std::vector<Variant>& the_candidates)
    {
        VariantCalls<Octopus::SampleIdType, Octopus::ProbabilityType> result {};
        
        for (const auto& candidate : the_candidates) {
            Allele reference_allele {candidate.get_reference_allele()};
            Allele alternative_allele {candidate.get_alternative_allele()};
            
            auto it = find_phased_region(the_phased_regions, alternative_allele);
            
            for (const auto& sample : the_samples ) {
                auto p = probability_allele_in_sample(alternative_allele, it->the_haplotypes,
                                                      it->the_latent_posteriors.genotype_probabilities.at(sample),
                                                      it->the_genotypes);
                
                if (p > 0.9) {
                    result[sample][candidate] = p;
                }
            }
        }
        
        return result;
    }
    
    } // end namespace VariantCaller
    
} // end namespace Octopus
