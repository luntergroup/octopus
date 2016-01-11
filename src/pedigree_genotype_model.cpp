//
//  pedigree_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "pedigree_genotype_model.hpp"

#include <algorithm>

#include "haplotype_likelihood_cache.hpp"
#include "read_model.hpp"
#include "maths.hpp"
#include "read_utils.hpp"

namespace Octopus
{
    namespace GenotypeModel
    {
        using GenotypeLogLikelihood        = double;
        using SampleGenotypeLogLikelihoods = std::vector<GenotypeLogLikelihood>;
        struct GenotypeLogLikelihoods
        {
            GenotypeLogLikelihoods(SampleGenotypeLogLikelihoods mother, SampleGenotypeLogLikelihoods father,
                                   SampleGenotypeLogLikelihoods child)
            : mother {mother}, father {father}, child {child} {}
            SampleGenotypeLogLikelihoods mother, father, child;
        };
        
        using GenotypePosterior        = double;
        using SampleGenotypePosteriors = std::vector<GenotypePosterior>;
        struct GenotypePosteriors
        {
            GenotypePosteriors(SampleGenotypePosteriors mother, SampleGenotypePosteriors father,
                               SampleGenotypePosteriors child)
            : mother {mother}, father {father}, child {child} {}
            SampleGenotypePosteriors mother, father, child;
        };
        
        GenotypeLogLikelihoods
        compute_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                         const ReadMap& reads, ReadModel& read_model,
                                         const SampleIdType& mother, const SampleIdType& father,
                                         const SampleIdType& child)
        {
            using std::cbegin; using std::cend; using std::begin; using std::transform;
            
            SampleGenotypeLogLikelihoods mother_log_likelihoods(genotypes.size());
            
            const auto& mother_reads = reads.at(mother);
            
            transform(cbegin(genotypes), cend(genotypes), begin(mother_log_likelihoods),
                      [&mother_reads, &read_model] (const auto& genotype) {
                          return read_model.log_probability(cbegin(mother_reads), cend(mother_reads),
                                                            genotype);
                      });
            
            SampleGenotypeLogLikelihoods father_log_likelihoods(genotypes.size());
            
            const auto& father_reads = reads.at(father);
            
            transform(cbegin(genotypes), cend(genotypes), begin(father_log_likelihoods),
                      [&father_reads, &read_model] (const auto& genotype) {
                          return read_model.log_probability(cbegin(father_reads), cend(father_reads),
                                                            genotype);
                      });
            
            SampleGenotypeLogLikelihoods child_log_likelihoods(genotypes.size());
            
            const auto& child_reads = reads.at(child);
            
            transform(cbegin(genotypes), cend(genotypes), begin(child_log_likelihoods),
                      [&child_reads, &read_model] (const auto& genotype) {
                          return read_model.log_probability(cbegin(child_reads), cend(child_reads),
                                                            genotype);
                      });
            
            return GenotypeLogLikelihoods {
                std::move(mother_log_likelihoods),
                std::move(father_log_likelihoods),
                std::move(child_log_likelihoods)
            };
        }
        
//        GenotypePosteriors
//        init_genotype_posteriors(const GenotypeLogMarginals& genotype_log_marginals,
//                                 const GenotypeLogLikelihoods& genotype_log_likilhoods)
//        {
//            
//        }
        
        Pedigree::Latents
        Pedigree::evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads,
                           const ReferenceGenome& reference)
        {
            if (haplotypes.size() == 1) {
                // TODO: catch this case to avoid computing
            }
            
            HaplotypeLikelihoodCache haplotype_likelihoods {reads, haplotypes};
            
            ReadModel read_model {ploidy_, haplotype_likelihoods};
            
            auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
            
            std::cout << "there are " << genotypes.size() << " candidate genotypes" << std::endl;
            
//            const auto genotype_log_likilhoods = compute_genotype_log_likelihoods(genotypes, reads,
//                                                                                  read_model, mother_,
//                                                                                  father_, child_);
//            
//            
            
            Pedigree::Latents result {};
            return result;
        }
    } // namespace GenotypeModel
} // namespace Octopus
