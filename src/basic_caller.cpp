//
//  basic_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "basic_caller.h"

#include <unordered_map>

#include "genomic_region.h"
#include "read_manager.h"
#include "variant.h"
#include "haplotype.h"
#include "genotype.h"
#include "haplotype_tree.h"
#include "search_regions.h"
#include "vcf_record.h"

#include "genotype_model.h"
#include "population_genotype_model.h"

#include <iostream> // TEST

BasicVariantCaller::BasicVariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                                       ReadFilter read_filter, ReadTransform read_transform,
                                       CandidateVariantGenerator& candidate_generator)
:
VariantCaller {reference, read_manager, read_filter, read_transform, candidate_generator}
{}

GenomicRegion BasicVariantCaller::get_init_region(const GenomicRegion& region)
{
    return region;
}

GenomicRegion BasicVariantCaller::get_next_region(const GenomicRegion& current_region)
{
    return GenomicRegion {"TEST", 0, 0};
}

std::unordered_map<Haplotype, double>
get_haplotype_posteriors(const Octopus::GenotypeModel::SampleGenotypeProbabilities& genotype_posteriors)
{
    std::unordered_map<Haplotype, double> result {};
    
    for (const auto& genotype_posterior: genotype_posteriors) {
        for (const Haplotype& haplotype : genotype_posterior.first) {
            result[haplotype] += genotype_posterior.second;
        }
    }
    
    return result;
}



std::vector<VcfRecord> BasicVariantCaller::call_variants(const GenomicRegion& region,
                                                         const std::vector<Variant>& candidates,
                                                         const ReadMap& reads)
{
    std::vector<VcfRecord> result {};
    
    Octopus::HaplotypeTree tree {reference_};
    extend_tree(candidates, tree);
    auto haplotypes = tree.get_haplotypes(region);
    
    std::cout << "there are " << haplotypes.size() << " haplotypes" << std::endl;
    
    auto genotype_model = std::make_unique<Octopus::PopulationGenotypeModel>(1, 2);
    
    auto genotype_posteriors = genotype_model->evaluate(haplotypes, reads);
    
    for (const auto sample_genotype_posteriors : genotype_posteriors) {
        auto haplotype_posteriors = get_haplotype_posteriors(sample_genotype_posteriors.second);
        
        auto max = *std::max_element(std::cbegin(haplotype_posteriors), std::cend(haplotype_posteriors), [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
        
        std::cout << sample_genotype_posteriors.first << std::endl;
        std::cout << max.first << " " << max.second << std::endl;
    }
    
    return result;
}
