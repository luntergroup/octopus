// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef coalescent_model_hpp
#define coalescent_model_hpp

#include <vector>
#include <unordered_map>
#include <functional>
#include <utility>
#include <functional>
#include <iterator>
#include <algorithm>
#include <cstddef>
#include <tuple>
#include <cassert>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>

#include "core/types/haplotype.hpp"
#include "core/types/variant.hpp"

namespace octopus {

class CoalescentModel
{
public:
    struct Parameters
    {
        double snp_heterozygosity   = 0.001;
        double indel_heterozygosity = 0.0001;
    };
    
    enum class CachingStrategy { none, value, address };
    
    CoalescentModel() = delete;
    
    CoalescentModel(Haplotype reference,
                    Parameters parameters,
                    std::size_t num_haplotyes_hint = 1024,
                    CachingStrategy caching = CachingStrategy::value);
    
    CoalescentModel(const CoalescentModel&)            = default;
    CoalescentModel& operator=(const CoalescentModel&) = default;
    CoalescentModel(CoalescentModel&&)                 = default;
    CoalescentModel& operator=(CoalescentModel&&)      = default;
    
    ~CoalescentModel() = default;
    
    void set_reference(Haplotype reference);
    
    void prime(std::vector<Haplotype> haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    template <typename Container>
    double evaluate(const Container& haplotypes) const;
    
    double evaluate(const std::vector<unsigned>& haplotype_indices) const;
    
private:
    using VariantReference = std::reference_wrapper<const Variant>;
    using SiteCountTuple = std::tuple<unsigned, unsigned, unsigned>;
    using SiteCountIndelTuple = std::tuple<unsigned, unsigned, unsigned, int>;
    
    struct SiteCountTupleHash
    {
        std::size_t operator()(const SiteCountIndelTuple& t) const noexcept
        {
            return boost::hash_value(t);
        }
    };
    
    Haplotype reference_;
    std::vector<double> reference_base_indel_heterozygosities_;
    Parameters params_;
    std::vector<Haplotype> haplotypes_;
    CachingStrategy caching_;
    
    mutable std::vector<VariantReference> site_buffer1_, site_buffer2_;
    mutable std::unordered_map<Haplotype, std::vector<Variant>> difference_value_cache_;
    mutable std::unordered_map<const Haplotype*, std::vector<Variant>> difference_address_cache_;
    mutable std::vector<boost::optional<std::vector<Variant>>> index_cache_;
    mutable std::vector<bool> index_flag_buffer_;
    mutable std::vector<std::vector<boost::optional<double>>> k_indel_zero_result_cache_;
    mutable std::unordered_map<SiteCountIndelTuple, double, SiteCountTupleHash> k_indel_pos_result_cache_;
    
    double evaluate(const SiteCountTuple& t) const;
    double evaluate(unsigned k_snp, unsigned n) const;
    double evaluate(unsigned k_snp, unsigned k_indel, unsigned n) const;
    
    template <typename Container>
    void fill_site_buffer(const Container& haplotypes) const;
    void fill_site_buffer(const std::vector<unsigned>& haplotype_indices) const;
    void fill_site_buffer_from_value_cache(const Haplotype& haplotype) const;
    void fill_site_buffer_from_address_cache(const Haplotype& haplotype) const;
    
    template <typename Container>
    SiteCountTuple count_segregating_sites(const Container& haplotypes) const;
};

template <typename Container>
double CoalescentModel::evaluate(const Container& haplotypes) const
{
    const auto t = count_segregating_sites(haplotypes);
    return evaluate(t);
}

// private methods

template <typename Container>
void CoalescentModel::fill_site_buffer(const Container& haplotypes) const
{
    assert(site_buffer2_.empty());
    site_buffer1_.clear();
    for (const Haplotype& haplotype : haplotypes) {
        if (caching_ == CachingStrategy::address) {
            fill_site_buffer_from_address_cache(haplotype);
        } else {
            fill_site_buffer_from_value_cache(haplotype);
        }
        site_buffer1_ = std::move(site_buffer2_);
        site_buffer2_.clear();
    }
}

namespace detail {

template <typename Container>
auto size(const Container& haplotypes) noexcept
{
    // Use this because Genotype template does not have a size member method (uses ploidy instead).
    return std::distance(std::cbegin(haplotypes), std::cend(haplotypes));
}

} // namespace detail

template <typename Container>
CoalescentModel::SiteCountTuple
CoalescentModel::count_segregating_sites(const Container& haplotypes) const
{
    fill_site_buffer(haplotypes);
    const auto num_indels = std::count_if(std::cbegin(site_buffer1_), std::cend(site_buffer1_),
                                          [] (const auto& v) noexcept { return is_indel(v); });
    return std::make_tuple(site_buffer1_.size() - num_indels, num_indels,
                           static_cast<unsigned>(detail::size(haplotypes) + 1));
}

} // namespace octopus

#endif
