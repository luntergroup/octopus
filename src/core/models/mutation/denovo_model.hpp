// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef denovo_model_hpp
#define denovo_model_hpp

#include <cstddef>
#include <unordered_map>
#include <string>
#include <utility>

#include <boost/optional.hpp>
#include <boost/functional/hash.hpp>

#include "core/types/haplotype.hpp"
#include "../pairhmm/pair_hmm.hpp"

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
    struct AddressPairHash
    {
        std::size_t operator()(const std::pair<const Haplotype*, const Haplotype*>& p) const noexcept
        {
            auto result = boost::hash_value(p.first);
            boost::hash_combine(result, p.second);
            return result;
        }
    };
    
    hmm::BasicMutationModel mutation_model_;
    std::size_t num_haplotypes_hint_;
    CachingStrategy caching_;
    mutable std::unordered_map<Haplotype, std::unordered_map<Haplotype, double>> value_cache_;
    mutable std::unordered_map<std::pair<const Haplotype*, const Haplotype*>, double, AddressPairHash> address_cache_;
    mutable std::string padded_given_;
    
    double evaluate_uncached(const Haplotype& target, const Haplotype& given) const;
    double evaluate_basic_cache(const Haplotype& target, const Haplotype& given) const;
    double evaluate_address_cache(const Haplotype& target, const Haplotype& given) const;
};

} // namespace octopus

 
#endif
