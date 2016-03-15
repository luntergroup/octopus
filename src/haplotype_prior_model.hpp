//
//  haplotype_prior_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_haplotype_prior_model_hpp
#define Octopus_haplotype_prior_model_hpp

#include <vector>
#include <unordered_map>
#include <functional>
#include <iterator>
#include <cstddef>
#include <algorithm>
#include <cassert>

#include "haplotype.hpp"

namespace Octopus
{
class HaplotypePriorModel
{
public:
    using HaplotypePriorMap  = std::unordered_map<std::reference_wrapper<const Haplotype>, double>;
    
    virtual ~HaplotypePriorModel() = default;
    
    double evaluate(const Haplotype& to, const Haplotype& from) const;
    
    HaplotypePriorMap
    compute_maximum_entropy_haplotype_set(std::vector<Haplotype>& haplotypes) const;
    
private:
    // ln p(to | from)
    virtual double do_evaluate(const Haplotype& to, const Haplotype& from) const = 0;
    
    virtual HaplotypePriorMap
    do_compute_maximum_entropy_haplotype_set(std::vector<Haplotype>& haplotypes) const = 0;
};

namespace debug
{
    template <typename S>
    void print_haplotype_priors(S&& stream,
                                const HaplotypePriorModel::HaplotypePriorMap& haplotype_priors,
                                const std::size_t n = 5)
    {
        const auto m = std::min(haplotype_priors.size(), n);
        
        stream << "Printing top " << m << " haplotype priors" << '\n';
        
        using HaplotypeReference = std::reference_wrapper<const Haplotype>;
        
        std::vector<std::pair<HaplotypeReference, double>> v {};
        v.reserve(haplotype_priors.size());
        
        std::copy(std::cbegin(haplotype_priors), std::cend(haplotype_priors), std::back_inserter(v));
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&] (const auto& p) {
                          ::debug::print_variant_alleles(stream, p.first);
                          stream << " " << p.second << '\n';
                      });
    }
    void print_haplotype_priors(const HaplotypePriorModel::HaplotypePriorMap& haplotype_priors,
                                std::size_t n = 5);
} // namespace debug
} // namespace Octopus

#endif
