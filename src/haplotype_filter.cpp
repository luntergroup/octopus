//
//  haplotype_filter.cpp
//  Octopus
//
//  Created by Daniel Cooke on 02/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "haplotype_filter.hpp"

#include <unordered_map>
#include <deque>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cassert>

#include "haplotype.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "maths.hpp"

namespace Octopus
{
    double max_read_likelihood(const Haplotype& haplotype,
                               const std::vector<SampleIdType>& samples,
                               const HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        auto result = std::numeric_limits<double>::lowest();
        
        for (const auto& sample : samples) {
            for (const double likelihood : haplotype_likelihoods.log_likelihoods(sample, haplotype)) {
                if (likelihood > result) result = likelihood;
                if (Maths::almost_zero(likelihood)) break;
            }
        }
        
        return result;
    }
    
    double sum_likelihoods(const Haplotype& haplotype,
                           const std::vector<SampleIdType>& samples,
                           const HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        double result {0};
        
        for (const auto& sample : samples) {
            const auto& sample_likelihoods = haplotype_likelihoods.log_likelihoods(sample, haplotype);
            result += std::accumulate(std::cbegin(sample_likelihoods), std::cend(sample_likelihoods), 0.0);
        }
        
        return result;
    }
    
    struct LikelihoodGreater
    {
        explicit LikelihoodGreater(const std::unordered_map<Haplotype, double>& liklihoods)
        : liklihoods_ {liklihoods} {}
        
        bool operator()(const Haplotype& lhs, const double rhs) const
        {
            return liklihoods_.at(lhs) > rhs;
        }
        
        bool operator()(const double& lhs, const Haplotype& rhs) const
        {
            return lhs > liklihoods_.at(rhs);
        }
        
        bool operator()(const Haplotype& lhs, const Haplotype& rhs) const
        {
            return liklihoods_.at(lhs) > liklihoods_.at(rhs);
        }
        
    private:
        const std::unordered_map<Haplotype, double>& liklihoods_;
    };
    
    template <typename RandomAccessIt>
    auto
    filter_by_likelihood_sum(const std::vector<SampleIdType>& samples,
                             const RandomAccessIt first, const RandomAccessIt last,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods,
                             const std::size_t n,
                             std::vector<Haplotype>& result)
    {
        const auto num_haplotypes = static_cast<std::size_t>(std::distance(first, last));
        
        std::unordered_map<Haplotype, double> liklihood_sums {num_haplotypes};
        
        std::for_each(first, last, [&] (const auto& haplotype) {
            liklihood_sums.emplace(haplotype, sum_likelihoods(haplotype, samples,
                                                              haplotype_likelihoods));
        });
        
        const auto nth = std::next(first, n);
        
        std::nth_element(first, nth, last, LikelihoodGreater {liklihood_sums});
        
        return nth;
    }
    
    std::vector<Haplotype>
    filter_by_maximum_likelihood(std::vector<Haplotype>& haplotypes,
                                 const std::vector<SampleIdType>& samples,
                                 const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                 const std::size_t n)
    {
        std::vector<Haplotype> result {};
        
        if (haplotypes.size() <= n) {
            std::sort(std::begin(haplotypes), std::end(haplotypes));
            return result;
        }
        
        std::unordered_map<Haplotype, double> max_liklihoods {haplotypes.size()};
        
        for (const auto& haplotype : haplotypes) {
            max_liklihoods.emplace(haplotype, max_read_likelihood(haplotype, samples, haplotype_likelihoods));
        }
        
        auto nth = std::next(std::begin(haplotypes), n);
        
        std::nth_element(std::begin(haplotypes), nth, std::end(haplotypes),
                         LikelihoodGreater {max_liklihoods});
        
        const auto nth_max_liklihood = max_liklihoods.at(*nth);
        
        const auto it = std::find_if(std::make_reverse_iterator(std::prev(nth)), std::rend(haplotypes),
                                     [nth_max_liklihood, &max_liklihoods] (const Haplotype& haplotype) {
                                         return max_liklihoods.at(haplotype) == nth_max_liklihood;
                                     });
        
        if (it != std::rend(haplotypes)) {
            std::sort(std::begin(haplotypes), nth, LikelihoodGreater {max_liklihoods});
            std::sort(nth, std::end(haplotypes), LikelihoodGreater {max_liklihoods});
            
            const auto er = std::equal_range(std::begin(haplotypes), std::end(haplotypes),
                                             nth_max_liklihood,
                                             LikelihoodGreater {max_liklihoods});
            
            result.assign(std::make_move_iterator(er.second),
                          std::make_move_iterator(std::end(haplotypes)));
            
            haplotypes.erase(er.second, std::end(haplotypes));
            
            // now deal with haplotypes in er
            nth = filter_by_likelihood_sum(samples, er.first, er.second, haplotype_likelihoods,
                                           n - std::distance(std::begin(haplotypes), er.first),
                                           result);
        }
        
        std::sort(std::begin(haplotypes), nth);
        std::sort(nth, std::end(haplotypes));
        
        std::deque<Haplotype> duplicates {};
        
        std::set_intersection(std::begin(haplotypes), nth, nth, std::end(haplotypes),
                              std::back_inserter(duplicates));
        
        result.insert(std::end(result),
                      std::make_move_iterator(nth),
                      std::make_move_iterator(std::end(haplotypes)));
        
        haplotypes.erase(nth, std::end(haplotypes));
        
        for (const auto& duplicate : duplicates) {
            const auto er = std::equal_range(std::begin(haplotypes), std::end(haplotypes), duplicate);
            haplotypes.erase(er.first, er.second);
        }
        
        return result;
    }
    
    std::vector<Haplotype>
    filter_to_n_haplotypes(std::vector<Haplotype>& haplotypes, const std::vector<SampleIdType>& samples,
                           const HaplotypeLikelihoodCache& haplotype_likelihoods, const std::size_t n)
    {
        return filter_by_maximum_likelihood(haplotypes, samples, haplotype_likelihoods, n);
    }
} // namespace Octopus