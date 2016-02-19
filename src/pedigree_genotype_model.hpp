//
//  pedigree_genotype_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 23/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef pedigree_genotype_model_hpp
#define pedigree_genotype_model_hpp

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>

#include "common.hpp"
#include "genotype.hpp"
#include "reference_genome.hpp"
#include "haplotype_prior_model.hpp"

class AlignedRead;
class Haplotype;

namespace Octopus
{
    namespace GenotypeModel
    {
        class Pedigree
        {
        public:
            struct Latents
            {
                // TODO
            };
            
            Pedigree(unsigned ploidy, unsigned max_em_iterations = 100, double em_epsilon = 0.001);
            
            Latents infer_latents(const std::vector<Haplotype>& haplotypes, const ReadMap& reads,
                                  const ReferenceGenome& reference);
            
        private:
            const unsigned ploidy_;
            const unsigned max_em_iterations_;
            const double em_epsilon_;
        };
        
    } // namespace GenotypeModel
} // namespace Octopus

#endif /* pedigree_genotype_model_hpp */
