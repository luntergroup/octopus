// Copyright (c) 2015-2021 Daniel Cooke
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
#include <type_traits>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>

#include "concepts/indexed.hpp"
#include "basics/tandem_repeat.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/variant.hpp"
#include "containers/mappable_block.hpp"
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
    
    void prime(MappableBlock<Haplotype> haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    // ln p(haplotype(s))
    LogProbability evaluate(const Haplotype& haplotype) const;
    template <typename Container> double evaluate(const Container& haplotypes) const;
    LogProbability evaluate(const std::vector<unsigned>& haplotype_indices) const;
    
private:
    using VariantReference = std::reference_wrapper<const Variant>;

    struct SegregatingSiteCounts
    {
        unsigned haplotypes;
        unsigned snps;
        unsigned repeat_indels, complex_indels;
    };
    struct SegregatingSiteCountsWithIndelHeterozygosities
    {
        SegregatingSiteCounts counts;
        double repeat_heterozygosity, complex_heterozygosity;
    };
    struct SegregatingSiteCountsHash
    {
        std::size_t operator()(const SegregatingSiteCounts& sites) const noexcept;
    };
    struct SegregatingSiteCountsWithIndelHeterozygositiesHash
    {
        std::size_t operator()(const SegregatingSiteCountsWithIndelHeterozygosities& sites) const noexcept;
    };
    struct SegregatingSiteCountsWithIndelHeterozygositiesEqual
    {
        bool operator()(const SegregatingSiteCountsWithIndelHeterozygosities& lhs,
                        const SegregatingSiteCountsWithIndelHeterozygosities& rhs) const noexcept;
    };

    using SegregatingSiteCountsWithIndelHeterozygositiesCache =  std::unordered_map<
        SegregatingSiteCountsWithIndelHeterozygosities, LogProbability,
        SegregatingSiteCountsWithIndelHeterozygositiesHash,
        SegregatingSiteCountsWithIndelHeterozygositiesEqual
        >;
    
    friend bool operator==(const SegregatingSiteCounts& lhs, const SegregatingSiteCounts& rhs) noexcept;

    Haplotype reference_;
    std::vector<TandemRepeat> reference_repeats_;
    IndelMutationModel::ContextIndelModel indel_heterozygosity_model_;
    Parameters params_;
    MappableBlock<Haplotype> haplotypes_;
    CachingStrategy caching_;
    
    mutable std::vector<VariantReference> site_buffer1_, site_buffer2_;
    mutable std::unordered_map<Haplotype, std::vector<Variant>> difference_value_cache_;
    mutable std::unordered_map<const Haplotype*, std::vector<Variant>> difference_address_cache_;
    mutable std::vector<boost::optional<std::vector<Variant>>> index_cache_;
    mutable std::vector<bool> index_flag_buffer_;
    mutable std::vector<std::vector<boost::optional<LogProbability>>> k_indel_zero_result_cache_;
    mutable SegregatingSiteCountsWithIndelHeterozygositiesCache k_indel_pos_result_cache_;
    
    LogProbability evaluate(const SegregatingSiteCounts& t) const;
    LogProbability evaluate_no_indels(unsigned k_snp, unsigned n) const;

    void fill_site_buffer(const Haplotype& haplotype) const;
    template <typename Range> void fill_site_buffer(const Range& haplotypes) const;
    template <typename Range> void fill_site_buffer(const Range& haplotypes, std::false_type) const;
    template <typename Range> void fill_site_buffer(const Range& haplotypes, std::true_type) const;
    void fill_site_buffer_uncached(const Haplotype& haplotype) const;
    void fill_site_buffer_from_value_cache(const Haplotype& haplotype) const;
    void fill_site_buffer_from_address_cache(const Haplotype& haplotype) const;
    
    SegregatingSiteCounts count_segregating_sites(const Haplotype& haplotype) const;
    template <typename Container> SegregatingSiteCounts count_segregating_sites(const Container& haplotypes) const;
    SegregatingSiteCounts count_segregating_sites_in_buffer(unsigned num_haplotypes) const;
    std::pair<double, double> calculate_buffered_indel_heterozygosities() const;
    double calculate_heterozygosity(const Variant& indel) const;
};

template <typename Container>
CoalescentModel::LogProbability CoalescentModel::evaluate(const Container& haplotypes) const
{
    return evaluate(count_segregating_sites(haplotypes));
}

// private methods

namespace detail {

template <typename T, typename = void> struct is_indexed_or_index : std::false_type {};
template <typename T>
struct is_indexed_or_index<T, std::enable_if_t<is_indexed_v<T> || std::is_integral<T>::value>> : std::true_type {};

} // namespace detail

template <typename Range>
void CoalescentModel::fill_site_buffer(const Range& haplotypes) const
{
    fill_site_buffer(haplotypes, detail::is_indexed_or_index<typename Range::value_type> {});
}

template <typename Range>
void CoalescentModel::fill_site_buffer(const Range& haplotypes, std::false_type) const
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

template <typename T>
inline std::enable_if_t<std::is_integral<T>::value, T> index_of(T i) noexcept { return i; }

template <typename Range>
void CoalescentModel::fill_site_buffer(const Range& haplotypes, std::true_type) const
{
    site_buffer1_.clear();
    std::fill(std::begin(index_flag_buffer_), std::end(index_flag_buffer_), false);
    unsigned num_unique_nonempty_indices {0};
    auto middle = std::begin(site_buffer1_);
    for (auto indexed : haplotypes) {
        const auto index = index_of(indexed);
        if (!index_flag_buffer_[index]) {
            auto& variants = index_cache_[index];
            if (!variants) {
                variants = haplotypes_[index].difference(reference_);
            }
            if (!variants->empty()) {
                middle = site_buffer1_.insert(std::cend(site_buffer1_), std::cbegin(*variants), std::cend(*variants));
                ++num_unique_nonempty_indices;
            }
            index_flag_buffer_[index] = true;
        }
    }
    if (num_unique_nonempty_indices == 2) {
        assert(site_buffer2_.empty());
        std::merge(std::begin(site_buffer1_), middle, middle, std::end(site_buffer1_), std::back_inserter(site_buffer2_));
        site_buffer2_.erase(std::unique(std::begin(site_buffer2_), std::end(site_buffer2_)), std::end(site_buffer2_));
        site_buffer1_ = std::move(site_buffer2_);
        site_buffer2_.clear();
    } else if (num_unique_nonempty_indices > 2) {
        std::sort(std::begin(site_buffer1_), std::end(site_buffer1_));
        site_buffer1_.erase(std::unique(std::begin(site_buffer1_), std::end(site_buffer1_)), std::end(site_buffer1_));
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
CoalescentModel::SegregatingSiteCounts
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
