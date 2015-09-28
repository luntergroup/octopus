//
//  cancer_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 16/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_caller.hpp"

#include <unordered_map>
#include <numeric>

#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "haplotype_tree.hpp"
#include "search_regions.hpp"
#include "vcf_record.hpp"

#include "genotype_model.hpp"
#include "cancer_genotype_model.hpp"

#include <iostream> // TEST

CancerVariantCaller::CancerVariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                                         ReadFilter read_filter, ReadTransform read_transform,
                                         CandidateVariantGenerator& candidate_generator)
:
VariantCaller {reference, read_manager, read_filter, read_transform, candidate_generator}
{}

GenomicRegion CancerVariantCaller::get_init_region(const GenomicRegion& region)
{
    return region;
}

GenomicRegion CancerVariantCaller::get_next_region(const GenomicRegion& current_region)
{
    return GenomicRegion {"TEST", 0, 0};
}

std::vector<VcfRecord> CancerVariantCaller::call_variants(const GenomicRegion& region,
                                                          const std::vector<Variant>& candidates,
                                                          const ReadMap& reads)
{
    std::vector<VcfRecord> result {};
    
    Octopus::HaplotypeTree tree {reference_};
    extend_tree(candidates, tree);
    
    auto haplotypes = tree.get_haplotypes(region);
    
    std::cout << "there are " << haplotypes.size() << " haplotypes" << std::endl;
    
    auto genotype_model = std::make_unique<Octopus::CancerGenotypeModel>(1, "1");
    
    auto genotype_posteriors = genotype_model->evaluate(haplotypes, reads);
    
    return result;
}
