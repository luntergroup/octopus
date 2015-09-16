//
//  basic_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "basic_caller.h"

#include <unordered_map>
#include <numeric>

#include "genomic_region.h"
#include "read_manager.h"
#include "allele.h"
#include "variant.h"
#include "haplotype.h"
#include "genotype.h"
#include "haplotype_tree.h"
#include "search_regions.h"
#include "vcf_record.h"

#include "mappable_algorithms.h"
#include "variant_utils.h"
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
get_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                         const Octopus::GenotypeModel::SampleGenotypeProbabilities& genotype_posteriors)
{
    std::unordered_map<Haplotype, double> result {};
    result.reserve(haplotypes.size());
    
    for (const auto& haplotype : haplotypes) {
        for (const auto& genotype_posterior: genotype_posteriors) {
            if (genotype_posterior.first.contains(haplotype)) {
                result[haplotype] += genotype_posterior.second;
            }
        }
    }
    
    return result;
}

double marginalise(const Allele& allele, const std::unordered_map<Haplotype, double>& haplotype_posteriors)
{
    double result {};
    
    for (const auto& haplotype_posterior : haplotype_posteriors) {
        if (haplotype_posterior.first.contains(allele)) {
            result += haplotype_posterior.second;
        }
    }
    
    return result;
}

std::unordered_map<Allele, double>
get_allele_posteriors(const std::unordered_map<Haplotype, double>& haplotype_posteriors,
                      const std::vector<Variant>& variants)
{
    auto alleles = decompose(variants);
    
    std::unordered_map<Allele, double> result {};
    result.reserve(alleles.size());
    
    for (const auto& allele : alleles) {
        result.emplace(allele, marginalise(allele, haplotype_posteriors));
    }
    
    return result;
}

double marginalise(const Allele& allele, const Octopus::GenotypeModel::SampleGenotypeProbabilities& genotype_posteriors)
{
    double result {};
    
    for (const auto& genotype_posterior : genotype_posteriors) {
        if (std::any_of(std::cbegin(genotype_posterior.first), std::cend(genotype_posterior.first),
                        [&allele] (const auto& haplotype) {
                            return haplotype.contains(allele);
                        })) {
                            result += genotype_posterior.second;
                        }
    }
    
    return result;
}

std::unordered_map<Allele, double>
get_allele_posteriors(const Octopus::GenotypeModel::SampleGenotypeProbabilities& genotype_posteriors,
                      const std::vector<Allele>& alleles)
{
    std::unordered_map<Allele, double> result {};
    result.reserve(alleles.size());
    
    for (const auto& allele : alleles) {
        result.emplace(allele, marginalise(allele, genotype_posteriors));
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
        auto haplotype_posteriors = get_haplotype_posteriors(haplotypes, sample_genotype_posteriors.second);
        
        std::cout << "haplotype posteriors" << std::endl;
        for (auto& h : haplotype_posteriors) {
            h.first.print_explicit_alleles();
            std::cout << h.second << std::endl;
        }
        
        auto alleles = decompose(candidates);
        
        //auto allele_posteriors = get_allele_posteriors(haplotype_posteriors, candidates);
        auto allele_posteriors = get_allele_posteriors(sample_genotype_posteriors.second, alleles);
        
        auto segments = segment(alleles);
        
        for (const auto& segment : segments) {
            
        }
        
        std::cout << "allele posteriors" << std::endl;
        for (auto& a : allele_posteriors) {
            std::cout << a.first << " " << a.second << std::endl;
        }
    }
    
    return result;
}
