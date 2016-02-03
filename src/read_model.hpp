//
//  read_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_model__
#define __Octopus__read_model__

#include <cstddef>
#include <numeric>
#include <functional>

#include "haplotype.hpp"
#include "genotype.hpp"
#include "aligned_read.hpp"
#include "haplotype_likelihood_cache.hpp"

namespace Octopus
{

class ReadModel
{
public:
    ReadModel()  = delete;
    explicit ReadModel(unsigned ploidy, const HaplotypeLikelihoodCache& haplotype_likelihoods);
    ~ReadModel() = default;
    
    ReadModel(const ReadModel&)            = default;
    ReadModel& operator=(const ReadModel&) = default;
    ReadModel(ReadModel&&)                 = default;
    ReadModel& operator=(ReadModel&&)      = default;
    
    // ln p(read | haplotype)
    double log_probability(const AlignedRead& read, const Haplotype& haplotype);
    
    // ln p(read | genotype)
    double log_probability(const AlignedRead& read, const Genotype<Haplotype>& genotype);
    
    // ln p(reads | genotype)
    template <typename ForwardIterator>
    double log_probability(ForwardIterator first_read, ForwardIterator last_read, const Genotype<Haplotype>& genotype);
    
private:
    std::reference_wrapper<const HaplotypeLikelihoodCache> haplotype_likelihoods_;
    
    const unsigned ploidy_;
    const double ln_ploidy_;
    
    // These are just for optimisation
    double log_probability_haploid(const AlignedRead& read, const Genotype<Haplotype>& genotype);
    double log_probability_diploid(const AlignedRead& read, const Genotype<Haplotype>& genotype);
    double log_probability_triploid(const AlignedRead& read, const Genotype<Haplotype>& genotype);
    double log_probability_polyploid(const AlignedRead& read, const Genotype<Haplotype>& genotype);
};

// ln p(reads | genotype) = sum (read in reads} ln p(read | genotype)
template <typename ForwardIterator>
double ReadModel::log_probability(ForwardIterator first_read, ForwardIterator last_read,
                                  const Genotype<Haplotype>& genotype)
{
    return std::accumulate(first_read, last_read, 0.0,
                           [this, &genotype] (auto curr, const auto& read) {
                               return curr + log_probability(read, genotype);
                           });
}

} // namespace Octopus

#endif /* defined(__Octopus__read_model__) */
