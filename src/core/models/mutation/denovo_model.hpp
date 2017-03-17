// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef denovo_model_hpp
#define denovo_model_hpp

#include <cstddef>
#include <unordered_map>
#include <utility>

#include <boost/optional.hpp>
#include <boost/functional/hash.hpp>

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
    
    template <typename Container> void prime(const Container& haplotypes) const;
    void unprime() const noexcept;
    
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
    
    Parameters parameters_;
    std::size_t num_haplotypes_hint_;
    CachingStrategy caching_;
    std::size_t max_flat_spaces_ = 50'000;
    
    mutable std::unordered_map<Haplotype, std::unordered_map<Haplotype, double>> value_cache_;
    mutable std::unordered_map<std::pair<const Haplotype*, const Haplotype*>, double, AddressPairHash> address_cache_;
    mutable boost::optional<const Haplotype*> first_haplotype_ = boost::none, last_haplotype_ = boost::none;
    mutable std::vector<std::vector<boost::optional<double>>> flat_address_cache_ = {};
    
    double evaluate_uncached(const Haplotype& target, const Haplotype& given) const;
    double evaluate_basic_cache(const Haplotype& target, const Haplotype& given) const;
    double evaluate_address_cache(const Haplotype& target, const Haplotype& given) const;
    double evaluate_flat_address_cache(const Haplotype& target, const Haplotype& given) const noexcept;
    double evaluate_sparse_address_cache(const Haplotype& target, const Haplotype& given) const;
    void do_prime(const std::vector<const Haplotype*>& addresses) const;
};

template <typename Container>
void DeNovoModel::prime(const Container& haplotypes) const
{
    std::vector<const Haplotype*> addresses(haplotypes.size());
    std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::begin(addresses),
                   [] (const Haplotype& haplotype) { return std::addressof(haplotype); });
    do_prime(addresses);
}

} // namespace octopus

#endif
