//
//  fixed_ploidy_genotype_likelihood_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__fixed_ploidy_genotype_likelihood_model__
#define __Octopus__fixed_ploidy_genotype_likelihood_model__

#include <cstddef>
#include <numeric>
#include <functional>
#include <iterator>
#include <vector>

#include "haplotype.hpp"
#include "genotype.hpp"
#include "aligned_read.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "probability_matrix.hpp"

namespace Octopus
{
namespace GenotypeModel
{
    class FixedPloidyGenotypeLikelihoodModel
    {
    public:
        FixedPloidyGenotypeLikelihoodModel()  = delete;
        explicit FixedPloidyGenotypeLikelihoodModel(unsigned ploidy, const HaplotypeLikelihoodCache& haplotype_likelihoods);
        ~FixedPloidyGenotypeLikelihoodModel() = default;
        
        FixedPloidyGenotypeLikelihoodModel(const FixedPloidyGenotypeLikelihoodModel&)            = default;
        FixedPloidyGenotypeLikelihoodModel& operator=(const FixedPloidyGenotypeLikelihoodModel&) = default;
        FixedPloidyGenotypeLikelihoodModel(FixedPloidyGenotypeLikelihoodModel&&)                 = default;
        FixedPloidyGenotypeLikelihoodModel& operator=(FixedPloidyGenotypeLikelihoodModel&&)      = default;
        
        // ln p(read | haplotype)
        double log_probability(const AlignedRead& read, const Haplotype& haplotype) const;
        
        // ln p(read | genotype)
        double log_probability(const AlignedRead& read, const Genotype<Haplotype>& genotype) const;
        
        // ln p(reads | genotype)
        template <typename ForwardIt>
        double log_probability(ForwardIt first_read, ForwardIt last_read,
                               const Genotype<Haplotype>& genotype) const;
        
    private:
        std::reference_wrapper<const HaplotypeLikelihoodCache> haplotype_likelihoods_;
        
        unsigned ploidy_;
        double ln_ploidy_;
        
        // These are just for optimisation
        double log_probability_haploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const;
        double log_probability_diploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const;
        double log_probability_triploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const;
        double log_probability_polyploid(const AlignedRead& read, const Genotype<Haplotype>& genotype) const;
    };
    
    // ln p(reads | genotype) = sum (read in reads} ln p(read | genotype)
    template <typename ForwardIt>
    double FixedPloidyGenotypeLikelihoodModel::log_probability(ForwardIt first_read, ForwardIt last_read,
                                                               const Genotype<Haplotype>& genotype) const
    {
        // These cases are just for optimisation
        switch (genotype.ploidy()) {
            case 1:
                return std::accumulate(first_read, last_read, 0.0,
                                       [this, &genotype] (const auto curr, const auto& read) {
                                           return curr + log_probability_haploid(read, genotype);
                                       });
            case 2:
                return std::accumulate(first_read, last_read, 0.0,
                                       [this, &genotype] (const auto curr, const auto& read) {
                                           return curr + log_probability_diploid(read, genotype);
                                       });
            case 3:
                return std::accumulate(first_read, last_read, 0.0,
                                       [this, &genotype] (const auto curr, const auto& read) {
                                           return curr + log_probability_triploid(read, genotype);
                                       });
            default:
                return std::accumulate(first_read, last_read, 0.0,
                                       [this, &genotype] (const auto curr, const auto& read) {
                                           return curr + log_probability_polyploid(read, genotype);
                                       });
        }
    }
    
    template <typename Container>
    auto log_probability(const Container& reads, const Genotype<Haplotype>& genotype,
                         const FixedPloidyGenotypeLikelihoodModel& read_model)
    {
        return read_model.log_probability(std::cbegin(reads), std::cend(reads), genotype);
    }
    
    ProbabilityMatrix<Genotype<Haplotype>>
    compute_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                     const ReadMap& reads,
                                     const FixedPloidyGenotypeLikelihoodModel& read_model);
    
    namespace debug
    {
        void print_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                            const FixedPloidyGenotypeLikelihoodModel& read_model,
                                            size_t n = 5);
        
        void print_read_genotype_liklihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                            const ReadMap& reads,
                                            const FixedPloidyGenotypeLikelihoodModel& read_model,
                                            size_t n = 3);
    } // namespace debug
} // namespace GenotypeModel
} // namespace Octopus

#endif /* defined(__Octopus__fixed_ploidy_genotype_likelihood_model__) */
