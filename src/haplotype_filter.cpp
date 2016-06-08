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
#include <type_traits>
#include <functional>
#include <limits>
#include <cassert>

#include "common.hpp"
#include "haplotype.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "maths.hpp"
#include "logging.hpp"

namespace Octopus
{
namespace
{
template <typename T>
struct FilterGreater
{
    FilterGreater(const std::unordered_map<Haplotype, T>& values) : values_ {values} {}
    
    bool operator()(const Haplotype& lhs, const T rhs) const
    {
        return values_.at(lhs) > rhs;
    }
    
    bool operator()(const T& lhs, const Haplotype& rhs) const
    {
        return lhs > values_.at(rhs);
    }
    
    bool operator()(const Haplotype& lhs, const Haplotype& rhs) const
    {
        return values_.at(lhs) > values_.at(rhs);
    }
    
private:
    const std::unordered_map<Haplotype, T>& values_;
};

template <typename F>
std::size_t try_filter(std::vector<Haplotype>& haplotypes,
                       const std::vector<SampleIdType>& samples,
                       const HaplotypeLikelihoodCache& haplotype_likelihoods,
                       const std::size_t n,
                       std::vector<Haplotype>& result,
                       F filter)
{
    using T = std::result_of_t<F(const Haplotype&, decltype(samples), decltype(haplotype_likelihoods))>;
    
    std::unordered_map<Haplotype, T> filter_likelihoods {haplotypes.size()};
    
    for (const auto& haplotype : haplotypes) {
        filter_likelihoods.emplace(haplotype, filter(haplotype, samples, haplotype_likelihoods));
    }
    
    const auto first = std::begin(haplotypes);
    const auto last  = std::end(haplotypes);
    
    const auto nth = std::next(first, n);
    
    const FilterGreater<T> cmp {filter_likelihoods};
    
    std::nth_element(first, nth, last, cmp);
    
    const auto nth_max_liklihood = filter_likelihoods.at(*nth);
    
    const auto rlast = std::make_reverse_iterator(first);
    
    const auto it = std::find_if(std::make_reverse_iterator(std::prev(nth)), rlast,
                                 [nth_max_liklihood, &filter_likelihoods] (const Haplotype& haplotype) {
                                     return filter_likelihoods.at(haplotype) == nth_max_liklihood;
                                 });
    
    if (it != rlast) {
        std::sort(first, nth, cmp);
        std::sort(nth, last, cmp);
        
        const auto er = std::equal_range(first, last, nth_max_liklihood, cmp);
        
        result.insert(std::end(result), std::make_move_iterator(er.second), std::make_move_iterator(last));
        
        const auto num_removed = std::distance(er.second, last);
        
        haplotypes.erase(er.second, last);
        
        return num_removed;
    }
    
    result.insert(std::end(result), std::make_move_iterator(nth), std::make_move_iterator(last));
    
    const auto num_removed = std::distance(nth, last);
    
    haplotypes.erase(nth, last);
    
    return num_removed;
}

template <typename F>
void force_filter(std::vector<Haplotype>& haplotypes,
                  const std::vector<SampleIdType>& samples,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods,
                  const std::size_t n,
                  std::vector<Haplotype>& result,
                  F filter)
{
    using T = std::result_of_t<F(const Haplotype&, decltype(samples), decltype(haplotype_likelihoods))>;
    
    std::unordered_map<Haplotype, T> filter_likelihoods {haplotypes.size()};
    
    for (const auto& haplotype : haplotypes) {
        filter_likelihoods.emplace(haplotype, filter(haplotype, samples, haplotype_likelihoods));
    }
    
    const auto first = std::begin(haplotypes);
    const auto last  = std::end(haplotypes);
    
    const auto nth = std::next(first, n);
    
    const FilterGreater<T> cmp {filter_likelihoods};
    
    std::nth_element(first, nth, last, cmp);
    
    result.insert(std::end(result), std::make_move_iterator(nth), std::make_move_iterator(last));
    
    haplotypes.erase(nth, last);
}

// filters

struct MaxLikelihood
{
    auto operator()(const Haplotype& haplotype, const std::vector<SampleIdType>& samples,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods) const
    {
        auto result = std::numeric_limits<double>::lowest();
        
        for (const auto& sample : samples) {
            for (const auto likelihood : haplotype_likelihoods.log_likelihoods(sample, haplotype)) {
                if (likelihood > result) result = likelihood;
                if (Maths::almost_zero(likelihood)) break;
            }
        }
        
        return result;
    }
};

struct LikelihoodZeroCount
{
    auto operator()(const Haplotype& haplotype, const std::vector<SampleIdType>& samples,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods) const
    {
        std::size_t result {0};
        
        for (const auto& sample : samples) {
            const auto& sample_likelihoods = haplotype_likelihoods.log_likelihoods(sample, haplotype);
            result += std::count_if(std::cbegin(sample_likelihoods), std::cend(sample_likelihoods),
                                    [] (const auto& likelihood) { return likelihood == 0; });
        }
        
        return result;
    }
};

class AssignmentCount
{
    std::unordered_map<Haplotype, float> assignments_;
    
public:
    AssignmentCount() = delete;
    
    AssignmentCount(const std::vector<Haplotype>& haplotypes,
                    const std::vector<SampleIdType>& samples,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        assignments_.reserve(haplotypes.size());
        
        for (const auto& haplotype : haplotypes) {
            assignments_.emplace(haplotype, 0.0);
        }
        
        std::vector<std::reference_wrapper<const Haplotype>> top {};
        
        for (const auto& sample : samples) {
            const auto n = haplotype_likelihoods.num_likelihoods(sample);
            
            for (std::size_t i {0}; i < n; ++i) {
                auto cur_max = std::numeric_limits<HaplotypeLikelihoodCache::Likelihoods::value_type>::lowest();
                
                for (const auto& haplotype : haplotypes) {
                    const auto p = haplotype_likelihoods.log_likelihoods(sample, haplotype)[i];
                    
                    if (Maths::almost_equal(p, cur_max)) {
                        top.emplace_back(haplotype);
                    } else if (p > cur_max) {
                        top.assign({haplotype});
                        cur_max = p;
                    }
                }
                
                const auto top_score = 1.0f / top.size();
                
                for (const auto& haplotype : top) {
                    assignments_.at(haplotype) += top_score;
                }
                
                top.clear();
            }
        }
    }
    
    auto operator()(const Haplotype& haplotype, const std::vector<SampleIdType>& samples,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods) const
    {
        return assignments_.at(haplotype);
    }
};

struct LikelihoodSum
{
    auto operator()(const Haplotype& haplotype, const std::vector<SampleIdType>& samples,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods) const
    {
        double result {0};
        
        for (const auto& sample : samples) {
            const auto& sample_likelihoods = haplotype_likelihoods.log_likelihoods(sample, haplotype);
            result += std::accumulate(std::cbegin(sample_likelihoods), std::cend(sample_likelihoods), 0.0);
        }
        
        return result;
    }
};
} // namespace

// main method

std::vector<Haplotype>
filter_to_n_haplotypes(std::vector<Haplotype>& haplotypes,
                       const std::vector<SampleIdType>& samples,
                       const HaplotypeLikelihoodCache& haplotype_likelihoods,
                       const std::size_t n)
{
    std::vector<Haplotype> result {};
    
    if (haplotypes.size() <= n) {
        return result;
    }
    
    auto num_to_filter = haplotypes.size() - n;
    
    result.reserve(haplotypes.size() - n);
    
    Logging::DebugLogger debug_log {};
    Logging::TraceLogger trace_log {};
    
    if (DEBUG_MODE) {
        stream(debug_log) << "Filtering " << num_to_filter << " of "
                            << haplotypes.size() << " haplotypes";
    }
    
    num_to_filter -= try_filter(haplotypes, samples, haplotype_likelihoods, n, result,
                                MaxLikelihood {});
    
    if (DEBUG_MODE) {
        stream(debug_log) << "There are " << haplotypes.size()
                            << " remaining haplotypes after maximum likelihood filtering";
    }
    
    if (num_to_filter == 0) {
        return result;
    }
    
    num_to_filter -= try_filter(haplotypes, samples, haplotype_likelihoods, n, result,
                                LikelihoodZeroCount {});
    
    if (DEBUG_MODE) {
        stream(debug_log) << "There are " << haplotypes.size()
                            << " remaining haplotypes after likelihood zero count filtering";
    }
    
    if (num_to_filter == 0) {
        return result;
    }
    
    num_to_filter -= try_filter(haplotypes, samples, haplotype_likelihoods, n, result,
                                AssignmentCount {haplotypes, samples, haplotype_likelihoods});
    
    if (DEBUG_MODE) {
        stream(debug_log) << "There are " << haplotypes.size()
                            << " remaining haplotypes after assignment count filtering";
    }
    
    if (num_to_filter == 0) {
        return result;
    }
    
    num_to_filter -= try_filter(haplotypes, samples, haplotype_likelihoods, n, result,
                                LikelihoodSum {});
    
    if (DEBUG_MODE) {
        stream(debug_log) << "There are " << haplotypes.size()
                            << " remaining haplotypes after likelihood sum filtering";
    }
    
    if (num_to_filter == 0) {
        return result;
    }
    
    if (DEBUG_MODE) {
        stream(debug_log) << "Force filtering " << haplotypes.size() << " haplotypes with assignment count filtering";
    }
    
    force_filter(haplotypes, samples, haplotype_likelihoods, n, result,
                 AssignmentCount {haplotypes, samples, haplotype_likelihoods});
    
    return result;
}
} // namespace Octopus