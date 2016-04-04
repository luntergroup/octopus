//
//  coalescent_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "coalescent_model.hpp"



namespace Octopus
{
    CoalescentModel::CoalescentModel(Haplotype reference_haplotype,
                                     double snp_heterozygosity,
                                     double indel_heterozygosity)
    :
    reference_haplotypes_ {std::move(reference_haplotype)},
    snp_heterozygosity_ {snp_heterozygosity},
    indel_heterozygosity_ {indel_heterozygosity}
    {}
    
    CoalescentModel::CoalescentModel(std::vector<Haplotype> reference_haplotypes,
                                     double snp_heterozygosity,
                                     double indel_heterozygosity)
    :
    reference_haplotypes_ {std::move(reference_haplotypes)},
    snp_heterozygosity_ {snp_heterozygosity},
    indel_heterozygosity_ {indel_heterozygosity}
    {
        difference_cache_.reserve(128);
        difference_cache_.emplace(std::cref(reference_haplotypes_.front()), std::vector<Variant> {});
    }
    
    template <typename M>
    unsigned calculate_num_segregating_sites(const std::vector<std::reference_wrapper<const Haplotype>>& haplotypes,
                                             M& difference_cache)
    {
        std::set<Variant> differences {};
        
        const auto& reference = haplotypes.front();
        
        std::for_each(std::next(std::cbegin(haplotypes)), std::cend(haplotypes),
                      [&] (const auto& haplotype) {
                          const auto it = difference_cache.find(haplotype);
                          
                          if (it != std::cend(difference_cache)) {
                              differences.insert(std::cbegin(it->second), std::cend(it->second));
                          } else {
                              auto curr = haplotype.get().difference(reference.get());
                              differences.insert(std::cbegin(curr), std::cend(curr));
                              difference_cache.emplace(std::cref(haplotype), std::move(curr));
                          }
                      });
        
        return static_cast<unsigned>(differences.size());
    }
    
//    double CoalescentModel::evaluate(const std::vector<Haplotype>& haplotypes) const
//    {
//        std::vector<std::reference_wrapper<const Haplotype>> augmented_haplotypes {};
//        augmented_haplotypes.reserve(reference_haplotypes_.size() + haplotypes.size());
//        
//        std::transform(std::cbegin(reference_haplotypes_), std::cend(reference_haplotypes_),
//                       std::back_inserter(augmented_haplotypes),
//                       [] (const auto& haplotype) { return std::cref(haplotype); });
//        std::transform(std::cbegin(haplotypes), std::cend(haplotypes),
//                       std::back_inserter(augmented_haplotypes),
//                       [] (const auto& haplotype) { return std::cref(haplotype); });
//        
//        const auto k = calculate_num_segregating_sites(augmented_haplotypes, difference_cache_);
//        
//        const auto n = static_cast<unsigned>(augmented_haplotypes.size());
//        
//        double result {0};
//        
//        const auto theta = snp_heterozygosity_;
//        
//        for (unsigned i {2}; i <= n; ++i) {
//            result += std::pow(-1, i) * boost::math::binomial_coefficient<double>(n - 1, i - 1)
//                        * ((i - 1) / (theta + i - 1)) * std::pow(theta / (theta + i - 1), k);
//        }
//        
//        return result;
//    }
} // namespace Octopus
