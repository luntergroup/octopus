//
//  coalescent_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef coalescent_model_hpp
#define coalescent_model_hpp

#include <vector>
#include <unordered_map>
#include <functional>
#include <utility>
#include <functional>
#include <iterator>
#include <cmath>
#include <algorithm>
#include <set>

#include <boost/math/special_functions/binomial.hpp>

#include "haplotype.hpp"
#include "variant.hpp"

namespace Octopus
{
    class CoalescentModel
    {
    public:
        CoalescentModel() = default;
        
        explicit CoalescentModel(Haplotype reference_haplotype,
                                 double snp_heterozygosity = 0.001,
                                 double indel_heterozygosity = 0.001);
        explicit CoalescentModel(std::vector<Haplotype> reference_haplotypes,
                                 double snp_heterozygosity = 0.001,
                                 double indel_heterozygosity = 0.001);
        
        ~CoalescentModel() = default;
        
        CoalescentModel(const CoalescentModel&)            = default;
        CoalescentModel& operator=(const CoalescentModel&) = default;
        CoalescentModel(CoalescentModel&&)                 = default;
        CoalescentModel& operator=(CoalescentModel&&)      = default;
        
        template <typename H>
        double evaluate(const H& haplotypes) const;
        
    private:
        std::vector<Haplotype> reference_haplotypes_;
        
        double snp_heterozygosity_, indel_heterozygosity_;
        
        using HaplotypeReference = std::reference_wrapper<const Haplotype>;
        
        mutable std::unordered_map<HaplotypeReference, std::vector<Variant>> difference_cache_;
    };
    
    template <typename H1, typename H2, typename M>
    unsigned calculate_num_segregating_sites(const H1& reference_haplotypes,
                                             const H2& candidate_haplotypes,
                                             M& difference_cache)
    {
        std::set<Variant> differences {};
        
        const auto& reference = reference_haplotypes.front();
        
        for (const auto& haplotype : candidate_haplotypes) {
            const auto it = difference_cache.find(haplotype);
            
            if (it != std::cend(difference_cache)) {
                differences.insert(std::cbegin(it->second), std::cend(it->second));
            } else {
                auto curr = haplotype.difference(reference);
                differences.insert(std::cbegin(curr), std::cend(curr));
                difference_cache.emplace(std::cref(haplotype), std::move(curr));
            }
        }
        
        return static_cast<unsigned>(differences.size());
    }
    
    template <typename H>
    double CoalescentModel::evaluate(const H& haplotypes) const
    {
        const auto k = calculate_num_segregating_sites(reference_haplotypes_, haplotypes,
                                                       difference_cache_);
        
        const auto n = static_cast<unsigned>(reference_haplotypes_.size() + haplotypes.ploidy());
        
        double result {0};
        
        const auto theta = snp_heterozygosity_;
        
        for (unsigned i {2}; i <= n; ++i) {
            result += std::pow(-1, i) * boost::math::binomial_coefficient<double>(n - 1, i - 1)
                        * ((i - 1) / (theta + i - 1)) * std::pow(theta / (theta + i - 1), k);
        }
        
        return result;
    }
} // namespace Octopus

#endif /* coalescent_model_hpp */
