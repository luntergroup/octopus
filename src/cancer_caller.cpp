//
//  cancer_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 16/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_caller.hpp"

#include <unordered_map>
#include <algorithm> // std::max_element, std::any_of
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

#include "mappable_algorithms.hpp"
#include "variant_utils.hpp"

#include "cancer_genotype.hpp"
#include "genotype_model.hpp"
#include "cancer_genotype_model.hpp"

#include <iostream> // TEST
#include "haplotype_prior_model.hpp"

namespace Octopus
{

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
    
    std::ostream& operator<<(std::ostream& os, const std::array<double, 3>& arr)
    {
        os << arr[0] << " " << arr[1] << " " << arr[2];
        return os;
    }
    
    std::ostream& operator<<(std::ostream& os, const std::unordered_map<Octopus::SampleIdType, std::array<double, 3>>& m)
    {
        for (const auto& p : m) os << p.first << ": " << p.second << "\n";
        return os;
    }
    
    static std::vector<GenomicRegion> get_segment_regions(const std::vector<std::vector<Allele>>& segments)
    {
        std::vector<GenomicRegion> result {};
        result.reserve(segments.size());
        for (const auto& segment : segments) result.push_back(segment.front().get_region());
        return result;
    }
    
    auto find_map_genotype(const GenotypeModel::Cancer::GenotypeProbabilities& genotype_posteriors)
    {
        return *std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                 [] (const auto& lhs, const auto& rhs) {
                                     return lhs.second < rhs.second;
                                 });
    }
    
    bool is_cancer(const GenotypeModel::Cancer::Latents& latents)
    {
        if (find_map_genotype(latents.genotype_posteriors).second > 0.95) {
            for (const auto& sample_genotype_weights : latents.genotype_weights) {
                if (sample_genotype_weights.second[2] > 0.2) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    bool is_homozygous_reference(const CancerGenotype<Haplotype>& genotype, const Haplotype& reference)
    {
        return genotype.get_normal_genotype().count(reference) == genotype.ploidy();
    }
    
    bool has_cancer_sample(const GenotypeModel::Cancer::GenotypeWeights& genotype_weights, double min_weight = 0.1)
    {
        return std::any_of(std::cbegin(genotype_weights), std::cend(genotype_weights),
                           [min_weight] (const auto& p) { return p.second[2] >= min_weight; });
    }
    
    
    
    std::vector<VcfRecord> CancerVariantCaller::call_variants(const GenomicRegion& region,
                                                              const std::vector<Variant>& candidates,
                                                              const ReadMap& reads)
    {
        std::vector<VcfRecord> result {};
        
        Octopus::HaplotypeTree tree {reference_};
        extend_tree(candidates, tree);
        
        auto haplotypes = tree.get_haplotypes(region);
        
        //std::cout << "there are " << haplotypes.size() << " haplotypes" << std::endl;
        
        SampleIdType normal_sample {"HG00102"};
        
        GenotypeModel::Cancer genotype_model {normal_sample};
        
        auto latents = genotype_model.evaluate(haplotypes, reads, reference_);
        
        auto map_genotype = find_map_genotype(latents.genotype_posteriors);
        
        if (map_genotype.second >= min_posterior_) {
            Haplotype reference_haplotype {reference_, region};
            
            auto has_variant = !is_homozygous_reference(map_genotype.first, reference_haplotype);
            auto has_cancer  = has_cancer_sample(latents.genotype_weights);
            
            if (has_variant || has_cancer) {
                if (has_variant && has_cancer) {
                    std::cout << "somatic mutation somewhere" << std::endl;
                } else if (has_variant && !has_cancer) {
                    std::cout << "probable normal variation" << std::endl;
                } else {
                    std::cout << "wierdness" << std::endl;
                }
            }
            
            auto alleles  = decompose(candidates);
            auto segments = segment(alleles);
            auto regions  = get_segment_regions(segments);
            
            
        }
        
        print_alleles(map_genotype.first);
        std::cout << map_genotype.second << std::endl;
        std::cout << latents.genotype_weights << std::endl;
        std::cout << std::boolalpha << is_cancer(latents) << std::endl;
        
        return result;
    }

} // namespace Octopus
