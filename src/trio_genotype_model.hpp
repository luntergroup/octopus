//
//  trio_genotype_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 29/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef trio_genotype_model_hpp
#define trio_genotype_model_hpp

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
        class Trio
        {
        public:
            struct Latents
            {
                // TODO
            };
            
            Trio(unsigned ploidy, unsigned max_em_iterations = 100, double em_epsilon = 0.001);
            
            Latents infer_latents(const std::vector<Haplotype>& haplotypes, const ReadMap& reads,
                                  const ReferenceGenome& reference);
            
        private:
            HaplotypePriorModel haplotype_prior_model_;
            
            const unsigned ploidy_;
            const unsigned max_em_iterations_;
            const double em_epsilon_;
        };
        
    } // namespace GenotypeModel
} // namespace Octopus

#endif /* trio_genotype_model_hpp */
