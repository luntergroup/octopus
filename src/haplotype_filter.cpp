//
//  haplotype_filter.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_filter.hpp"

#include <algorithm>
#include <iterator>
#include <utility>
#include <unordered_map>

#include "read_utils.hpp"

namespace Octopus
{
void filter_haplotypes(std::vector<Haplotype>& haplotypes, const ReadMap& reads, size_t n,
                       const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    if (haplotypes.size() <= n) return;
    
    std::unordered_map<std::reference_wrapper<const Haplotype>, double> max_liklihoods {haplotypes.size()};
    
    for (const auto& haplotype : haplotypes) {
        double max_read_liklihood {0};
        
        for (const auto& sample_reads : reads) {
            for (const auto& read : sample_reads.second) {
                const auto cur_read_liklihood = haplotype_likelihoods.log_probability(read, haplotype);
                if (cur_read_liklihood > max_read_liklihood) max_read_liklihood = cur_read_liklihood;
            }
        }
        
        max_liklihoods.emplace(haplotype, max_read_liklihood);
    }
    
    const auto nth = std::next(std::begin(haplotypes), n);
    
    std::nth_element(std::begin(haplotypes), nth, std::end(haplotypes),
                     [&] (const auto& lhs, const auto& rhs) {
                         return max_liklihoods.at(lhs) > max_liklihoods.at(rhs);
                     });
    
    haplotypes.erase(nth, std::end(haplotypes));
}

std::vector<Haplotype>
filter_haplotypes(const std::vector<Haplotype>& haplotypes, const ReadMap& reads, const size_t n,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    if (haplotypes.size() <= n) {
        return std::vector<Haplotype> {std::cbegin(haplotypes), std::cend(haplotypes)};
    }
    
    std::vector<std::pair<std::reference_wrapper<const Haplotype>, double>> haplotype_scores {};
    haplotype_scores.reserve(haplotypes.size());
    
    for (const auto& haplotype : haplotypes) {
        double max_read_liklihood {0};
        
        for (const auto& sample_reads : reads) {
            for (const auto& read : sample_reads.second) {
                const auto cur_read_liklihood = haplotype_likelihoods.log_probability(read, haplotype);
                if (cur_read_liklihood > max_read_liklihood) max_read_liklihood = cur_read_liklihood;
            }
        }
        
        haplotype_scores.emplace_back(haplotype, max_read_liklihood);
    }
    
    const auto nth = std::next(std::begin(haplotype_scores), n);
    
    std::nth_element(std::begin(haplotype_scores), nth, std::end(haplotype_scores),
                     [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
    
    std::vector<Haplotype> result {};
    result.reserve(n);
    
    std::transform(std::begin(haplotype_scores), nth, std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    
    return result;
}

} // namespace Octopus
