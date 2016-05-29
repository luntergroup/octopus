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
#include <cassert>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/functional/hash.hpp>

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
                                 double indel_heterozygosity = 0.001,
                                 unsigned max_haplotypes = 1024);
        
        ~CoalescentModel() = default;
        
        CoalescentModel(const CoalescentModel&)            = default;
        CoalescentModel& operator=(const CoalescentModel&) = default;
        CoalescentModel(CoalescentModel&&)                 = default;
        CoalescentModel& operator=(CoalescentModel&&)      = default;
        
        void set_reference(Haplotype reference);
        
        template <typename Container> double evaluate(const Container& haplotypes) const;
        
    private:
        using SiteCountTuple = std::tuple<unsigned, unsigned, unsigned>;
        
        struct SiteCountTupleHash
        {
            inline auto operator()(const SiteCountTuple& t) const noexcept
            {
                return boost::hash_value(t);
            }
        };
        
        Haplotype reference_;
        
        double snp_heterozygosity_, indel_heterozygosity_;
        
        using HaplotypeReference = std::reference_wrapper<const Haplotype>;
        
        mutable std::vector<std::reference_wrapper<const Variant>> site_buffer1_, site_buffer2_;
        
        using VariantReference = std::reference_wrapper<const Variant>;
        
        mutable std::unordered_map<HaplotypeReference, std::vector<Variant>> difference_cache_;
        
        mutable std::unordered_map<SiteCountTuple, double, SiteCountTupleHash> result_cache_;
        
        template <typename Container>
        SiteCountTuple count_segregating_sites(const Container& haplotypes) const;
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
        const auto t = count_segregating_sites(haplotypes);
        
        const auto it = result_cache_.find(t);
        
        if (it != std::cend(result_cache_)) return it->second;
        
        unsigned k_snp, k_indel, n;
        
        std::tie(k_snp, k_indel, n) = t;
        
        double result {0};
        
        const auto theta = snp_heterozygosity_;
        
        const auto k = k_snp + k_indel; // TODO: correct calculation for different indel heterozygosity
        
        for (unsigned i {2}; i <= n; ++i) {
            result += std::pow(-1, i) * boost::math::binomial_coefficient<double>(n - 1, i - 1)
                        * ((i - 1) / (theta + i - 1)) * std::pow(theta / (theta + i - 1), k);
        }
        
        result = std::log(result);
        
        result_cache_.emplace(t, result);
        
        return result;
    }
    
    // private methods
    
    template <typename Container>
    CoalescentModel::SiteCountTuple
    CoalescentModel::count_segregating_sites(const Container& haplotypes) const
    {
        assert(site_buffer2_.empty());
        
        site_buffer1_.clear();
        
        for (const Haplotype& haplotype : haplotypes) {
            auto it = difference_cache_.find(haplotype);
            
            if (it == cend(difference_cache_)) {
                it = difference_cache_.emplace(std::piecewise_construct,
                                               std::forward_as_tuple(haplotype),
                                               std::forward_as_tuple(haplotype.difference(reference_))).first;
            }
            
            std::set_union(std::begin(site_buffer1_), std::end(site_buffer1_),
                           std::cbegin(it->second), std::cend(it->second),
                           std::back_inserter(site_buffer2_));
            
            std::swap(site_buffer1_, site_buffer2_);
            
            site_buffer2_.clear();
        }
        
        const auto num_indels = std::count_if(std::cbegin(site_buffer1_), std::cend(site_buffer1_),
                                              [] (const auto& v) { return is_indel(v); });
        
        
        return std::make_tuple(site_buffer1_.size() - num_indels, num_indels,
                               static_cast<unsigned>(detail::size(haplotypes) + 1));
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
