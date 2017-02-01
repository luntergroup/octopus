// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef denovo_model_hpp
#define denovo_model_hpp

#include <cstddef>
#include <unordered_map>

#include <boost/optional.hpp>

#include "core/types/haplotype.hpp"

namespace octopus {

class DeNovoModel
{
public:
    struct Parameters
    {
        double mutation_rate;
    };
    
    enum class CachingStrategy { none, value, address };
    
    DeNovoModel() = delete;
    
    DeNovoModel(Parameters parameters,
                std::size_t num_haplotypes_hint = 1000,
                CachingStrategy caching = CachingStrategy::value);
    
    DeNovoModel(const DeNovoModel&)            = default;
    DeNovoModel& operator=(const DeNovoModel&) = default;
    DeNovoModel(DeNovoModel&&)                 = default;
    DeNovoModel& operator=(DeNovoModel&&)      = default;
    
    ~DeNovoModel() = default;
    
    // ln p(target | given)
    double evaluate(const Haplotype& target, const Haplotype& given) const;
    
private:
    Parameters parameters_;
    std::size_t num_haplotypes_hint_;
    CachingStrategy caching_;
    mutable std::unordered_map<Haplotype, std::unordered_map<Haplotype, double>> value_cache_;
    mutable std::unordered_map<const Haplotype*, std::unordered_map<const Haplotype*, double>> address_cache_;
    
    double evaluate_uncached(const Haplotype& target, const Haplotype& given) const;
    double evaluate_basic_cache(const Haplotype& target, const Haplotype& given) const;
    double evaluate_address_cache(const Haplotype& target, const Haplotype& given) const;
};

} // namespace octopus

 
#endif
