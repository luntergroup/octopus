//
//  cancer_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_genotype_model.hpp"

#include <algorithm>
#include <iterator>
#include <functional>
#include <cassert>
#include <iostream>

#include "maths.hpp"
#include "logging.hpp"

namespace Octopus
{
    namespace GenotypeModel
    {
    // public methods
    
    Cancer::Cancer(std::vector<SampleIdType> samples, const SampleIdType& normal_sample,
                   Priors priors)
    :
    Cancer {std::move(samples), normal_sample, std::move(priors), AlgorithmParameters {}}
    {}
    
    Cancer::Cancer(std::vector<SampleIdType> samples, const SampleIdType& normal_sample,
                   Priors priors, AlgorithmParameters parameters)
    :
    samples_ {std::move(samples)},
    normal_sample_ {normal_sample},
    priors_ {std::move(priors)},
    parameters_ {parameters}
    {
        const auto it = std::find(std::begin(samples_), std::end(samples_), normal_sample_);
        assert(it != std::end(samples_));
        std::iter_swap(std::begin(samples_), it);
    }
    
//        Cancer::Latents
//        run_variational_bayes(const std::vector<CancerGenotype<Haplotype>>& genotypes)
//        {
//            
//        }
    
    Cancer::Latents
    Cancer::infer_latents(const std::vector<Haplotype>& haplotypes,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
    {
        assert(!haplotypes.empty());
        
        static_assert(K > 1, "K <= 1");
        
        auto genotypes = generate_all_cancer_genotypes(haplotypes, K - 1);
        
        assert(!genotypes.empty());
        
        if (DEBUG_MODE) {
            Logging::DebugLogger log {};
            stream(log) << "There are " << genotypes.size() << " initial cancer genotypes";
        }
        
        exit(0);
        return Latents {};
    }
    
    namespace debug
    {
        
    } // namespace debug
    
    } // namespace GenotypeModel
} // namespace Octopus
