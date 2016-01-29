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
#include <functional>

#include "haplotype_likelihood_cache.hpp"
#include "read_utils.hpp"

namespace Octopus
{
std::vector<Haplotype> filter_haplotypes(const std::vector<Haplotype>& haplotypes,
                                         const ReadMap& reads, size_t n,
                                         HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    if (haplotypes.size() <= n) return haplotypes;
    
    std::vector<std::pair<std::reference_wrapper<const Haplotype>, double>> haplotype_scores {};
    haplotype_scores.reserve(haplotypes.size());
    
    for (const auto& haplotype : haplotypes) {
        double max {};
        
        for (const auto& sample_reads : reads) {
            for (const auto& read : sample_reads.second) {
                auto curr = haplotype_likelihoods.log_probability(read, haplotype);
                if (curr > max) max = curr;
            }
        }
        
        haplotype_scores.emplace_back(haplotype, max);
    }
    
    std::sort(std::begin(haplotype_scores), std::end(haplotype_scores),
              [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
    
    std::vector<Haplotype> result {};
    result.reserve(n);
    
    std::transform(std::cbegin(haplotype_scores), std::next(std::cbegin(haplotype_scores), n),
                   std::back_inserter(result), [] (const auto& p) { return p.first; });
    
    return result;
}

} // namespace Octopus
