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
#include <cstddef>
#include <tuple>

#include <boost/math/special_functions/binomial.hpp>

#include "haplotype.hpp"
#include "variant.hpp"
#include "maths.hpp"

namespace Octopus
{
    class CoalescentModel
    {
    public:
        CoalescentModel() = delete;
        
        explicit CoalescentModel(Haplotype reference,
                                 double snp_heterozygosity = 0.001,
                                 double indel_heterozygosity = 0.001);
        
        ~CoalescentModel() = default;
        
        CoalescentModel(const CoalescentModel&)            = default;
        CoalescentModel& operator=(const CoalescentModel&) = default;
        CoalescentModel(CoalescentModel&&)                 = default;
        CoalescentModel& operator=(CoalescentModel&&)      = default;
        
        void set_reference(Haplotype reference);
        
        template <typename Container> double evaluate(const Container& haplotypes) const;
        
    private:
        Haplotype reference_;
        
        double snp_heterozygosity_, indel_heterozygosity_;
        
        using HaplotypeReference = std::reference_wrapper<const Haplotype>;
        
        mutable std::unordered_map<HaplotypeReference, std::vector<Variant>> difference_cache_;
        
        template <typename Container>
        std::pair<unsigned, unsigned> count_segregating_sites(const Container& haplotypes) const;
    };
    
    namespace detail
    {
        template <typename Container>
        auto size(const Container& haplotypes)
        {
            // Use this because Genotype template does not have a size member method (uses
            // ploidy instead).
            return std::distance(std::cbegin(haplotypes), std::cend(haplotypes));
        }
    } // namespace detail
    
    template <typename Container>
    double CoalescentModel::evaluate(const Container& haplotypes) const
    {
        unsigned k_snp, k_indel;
        
        std::tie(k_snp, k_indel) = count_segregating_sites(haplotypes);
        
        const auto n = static_cast<unsigned>(detail::size(haplotypes) + 1); // + 1 for reference
        
        double result {0};
        
        const auto theta = snp_heterozygosity_;
        
        const auto k = k_snp + k_indel; // TODO: correct calculation for different indel heterozygosity
        
        for (unsigned i {2}; i <= n; ++i) {
            result += std::pow(-1, i) * boost::math::binomial_coefficient<double>(n - 1, i - 1)
                        * ((i - 1) / (theta + i - 1)) * std::pow(theta / (theta + i - 1), k);
        }
        
        return std::log(result);
    }
    
    // private methods
    
    template <typename Container>
    std::pair<unsigned, unsigned>
    CoalescentModel::count_segregating_sites(const Container& haplotypes) const
    {
        const auto n = detail::size(haplotypes);
        
        std::vector<std::reference_wrapper<const Variant>> sites {};
        sites.reserve(2 * n);
        
        for (const Haplotype& haplotype : haplotypes) {
            const auto it = difference_cache_.find(haplotype);
            
            if (it != std::cend(difference_cache_)) {
                sites.insert(std::end(sites), std::cbegin(it->second), std::cend(it->second));
            } else {
                const auto p = difference_cache_.emplace(std::piecewise_construct,
                                                         std::forward_as_tuple(haplotype),
                                                         std::forward_as_tuple(haplotype.difference(reference_)));
                
                sites.insert(std::end(sites), std::cbegin(p.first->second), std::cend(p.first->second));
            }
        }
        
        auto it1 = std::end(sites);
        
        if (n > 1) {
            std::sort(std::begin(sites), std::end(sites));
            
            it1 = std::unique(std::begin(sites), std::end(sites));
        }
        
        const auto it2 = std::partition(std::begin(sites), it1,
                                        [] (const auto& variant) {
                                            return is_indel(variant);
                                        });
        
        return std::make_pair(static_cast<unsigned>(std::distance(it2, it1)),
                              static_cast<unsigned>(std::distance(std::begin(sites), it2)));
    }
    
    // non-member methods
    
    template <typename Container>
    std::vector<double> calculate_log_priors(const Container& genotypes, const CoalescentModel& model)
    {
        std::vector<double> result(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [&model] (const auto& genotype) {
                           return model.evaluate(genotype);
                       });
        
        Maths::normalise_logs(result);
        
        return result;
    }
} // namespace Octopus

#endif /* coalescent_model_hpp */
