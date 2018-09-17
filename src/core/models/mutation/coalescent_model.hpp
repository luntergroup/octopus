// Copyright (c) 2015-2018 Daniel Cooke
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
#include "indel_mutation_model.hpp"

namespace octopus {

class CoalescentModel
{
public:
    using LogProbability = double;
    
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
    
    // ln p(haplotype(s))
    LogProbability evaluate(const Haplotype& haplotype) const;
    template <typename Container> double evaluate(const Container& haplotypes) const;
    LogProbability evaluate(const std::vector<unsigned>& haplotype_indices) const;
    
private:
    using VariantReference = std::reference_wrapper<const Variant>;
    using SiteCountTuple = std::tuple<unsigned, unsigned, unsigned>;
    using SiteCountIndelTuple = std::tuple<unsigned, unsigned, unsigned, double>;
    
    struct SiteCountTupleHash
    {
        std::size_t operator()(const SiteCountIndelTuple& t) const noexcept
        {
            return boost::hash_value(t);
        }
    };
    
    Haplotype reference_;
    IndelMutationModel::ContextIndelModel indel_heterozygosity_model_;
    Parameters params_;
    std::vector<Haplotype> haplotypes_;
    CachingStrategy caching_;
    
    mutable std::vector<VariantReference> site_buffer1_, site_buffer2_;
    mutable std::unordered_map<Haplotype, std::vector<Variant>> difference_value_cache_;
    mutable std::unordered_map<const Haplotype*, std::vector<Variant>> difference_address_cache_;
    mutable std::vector<boost::optional<std::vector<Variant>>> index_cache_;
    mutable std::vector<bool> index_flag_buffer_;
    mutable std::vector<std::vector<boost::optional<LogProbability>>> k_indel_zero_result_cache_;
    mutable std::unordered_map<SiteCountIndelTuple, LogProbability, SiteCountTupleHash> k_indel_pos_result_cache_;
    
    LogProbability evaluate(const SiteCountTuple& t) const;
    LogProbability evaluate(unsigned k_snp, unsigned n) const;
    LogProbability evaluate(unsigned k_snp, unsigned k_indel, unsigned n) const;
    
    void fill_site_buffer(const Haplotype& haplotype) const;
    template <typename Container> void fill_site_buffer(const Container& haplotypes) const;
    void fill_site_buffer(const std::vector<unsigned>& haplotype_indices) const;
    void fill_site_buffer_uncached(const Haplotype& haplotype) const;
    void fill_site_buffer_from_value_cache(const Haplotype& haplotype) const;
    void fill_site_buffer_from_address_cache(const Haplotype& haplotype) const;
    
    SiteCountTuple count_segregating_sites(const Haplotype& haplotype) const;
    template <typename Container> SiteCountTuple count_segregating_sites(const Container& haplotypes) const;
    SiteCountTuple count_segregating_sites_in_buffer(unsigned num_haplotypes) const;
    double calculate_buffered_indel_heterozygosity() const;
    double calculate_heterozygosity(const Variant& indel) const;
};

template <typename Container>
CoalescentModel::LogProbability CoalescentModel::evaluate(const Container& haplotypes) const
{
    return evaluate(count_segregating_sites(haplotypes));
}

// private methods

template <typename Container>
void CoalescentModel::fill_site_buffer(const Container& haplotypes) const
{
    assert(site_buffer2_.empty());
    site_buffer1_.clear();
    for (const Haplotype& haplotype : haplotypes) {
        switch (caching_) {
            case CachingStrategy::address: fill_site_buffer_from_address_cache(haplotype); break;
            case CachingStrategy::value: fill_site_buffer_from_value_cache(haplotype); break;
            default: fill_site_buffer_uncached(haplotype);
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
    return count_segregating_sites_in_buffer(detail::size(haplotypes));
}

struct CoalescentProbabilityGreater
{
    CoalescentProbabilityGreater(CoalescentModel model);
    
    bool operator()(const Haplotype& lhs, const Haplotype& rhs) const;

private:
    CoalescentModel model_;
    mutable std::vector<Haplotype> buffer_;
    mutable std::unordered_map<Haplotype, CoalescentModel::LogProbability, std::hash<Haplotype>, HaveSameAlleles> cache_;
};

} // namespace octopus

#endif
