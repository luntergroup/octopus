// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "coalescent_model.hpp"

#include <memory>
#include <cmath>
#include <complex>
#include <numeric>
#include <stdexcept>

#include <boost/math/special_functions/binomial.hpp>

#include "utils/maths.hpp"
#include "utils/mappable_algorithms.hpp"

namespace octopus {

CoalescentModel::CoalescentModel(Haplotype reference, Parameters params,
                                 std::size_t num_haplotyes_hint, CachingStrategy caching)
: reference_ {std::move(reference)}
, reference_repeats_ {}
, indel_heterozygosity_model_ {make_indel_model(reference_, {params.indel_heterozygosity}, reference_repeats_)}
, params_ {params}
, haplotypes_ {}
, caching_ {caching}
, index_cache_ {}
, index_flag_buffer_ {}
, k_indel_zero_result_cache_ {2 * num_haplotyes_hint, std::vector<boost::optional<LogProbability>> {}}
{
    if (params_.snp_heterozygosity <= 0 || params_.indel_heterozygosity <= 0) {
        throw std::domain_error {"CoalescentModel: snp and indel heterozygosity must be > 0"};
    }
    site_buffer1_.reserve(128);
    site_buffer2_.reserve(128);
    if (caching == CachingStrategy::address) {
        difference_address_cache_.reserve(num_haplotyes_hint);
    } else if (caching_ == CachingStrategy::value) {
        difference_value_cache_.reserve(num_haplotyes_hint);
        difference_value_cache_.emplace(std::piecewise_construct,
                                        std::forward_as_tuple(reference_),
                                        std::forward_as_tuple());
    }
    k_indel_pos_result_cache_.reserve(2 * num_haplotyes_hint);
}

void CoalescentModel::set_reference(Haplotype reference)
{
    reference_ = std::move(reference);
    if (caching_ == CachingStrategy::address) {
        difference_address_cache_.clear();
    } else if (caching_ == CachingStrategy::value) {
        difference_value_cache_.clear();
        difference_value_cache_.emplace(std::piecewise_construct,
                                        std::forward_as_tuple(reference_),
                                        std::forward_as_tuple());
    }
}

void CoalescentModel::prime(MappableBlock<Haplotype> haplotypes)
{
    haplotypes_ = std::move(haplotypes);
    index_cache_.assign(haplotypes_.size(), boost::none);
    index_flag_buffer_.assign(haplotypes_.size(), false);
}

void CoalescentModel::unprime() noexcept
{
    haplotypes_.clear();
    haplotypes_.shrink_to_fit();
    index_cache_.clear();
    index_cache_.shrink_to_fit();
    index_flag_buffer_.clear();
    index_flag_buffer_.shrink_to_fit();
}

bool CoalescentModel::is_primed() const noexcept
{
    return !index_cache_.empty();
}

CoalescentModel::LogProbability CoalescentModel::evaluate(const Haplotype& haplotype) const
{
    return evaluate(count_segregating_sites(haplotype));
}

CoalescentModel::LogProbability CoalescentModel::evaluate(const std::vector<unsigned>& haplotype_indices) const
{
    return evaluate(count_segregating_sites(haplotype_indices));
}

namespace {

auto powm1(const unsigned i) noexcept // std::pow(-1, i)
{
    return (i % 2 == 0) ? 1 : -1;
}

template <typename RealType>
auto coalescent_real_space(const unsigned n, const unsigned k, const RealType theta)
{
    RealType result {0};
    for (unsigned i {2}; i <= n; ++i) {
        result += powm1(i) * maths::binomial_coefficient<RealType>(n - 1, i - 1) * ((i - 1) / (theta + i - 1)) * std::pow(theta / (theta + i - 1), k);
    }
    return std::log(result);
}

template <typename ForwardIt>
auto complex_log_sum_exp(ForwardIt first, ForwardIt last)
{
    using ComplexType = typename std::iterator_traits<ForwardIt>::value_type;
    const auto l = [] (const auto& lhs, const auto& rhs) { return lhs.real() < rhs.real(); };
    const auto max = *std::max_element(first, last, l);
    return max + std::log(std::accumulate(first, last, ComplexType {},
                                          [max] (const auto curr, const auto x) { return curr + std::exp(x - max); }));
}

template <typename Container>
auto complex_log_sum_exp(const Container& logs)
{
    return complex_log_sum_exp(std::cbegin(logs), std::cend(logs));
}

template <typename RealType>
auto coalescent_log_space(const unsigned n, const unsigned k, const RealType theta)
{
    std::vector<std::complex<RealType>> tmp(n - 1, std::log(std::complex<RealType> {-1}));
    for (unsigned i {2}; i <= n; ++i) {
        auto& cur = tmp[i - 2];
        cur *= i;
        cur += maths::log_binomial_coefficient<RealType>(n - 1, i - 1);
        cur += std::log((i - 1) / (theta + i - 1));
        cur += k * std::log(theta / (theta + i - 1));
    }
    return complex_log_sum_exp(tmp).real();
}

template <typename RealType>
auto coalescent(const unsigned n, const unsigned k, const RealType theta)
{
    if (n < 30 && k <= 80) {
        auto result = coalescent_real_space(n, k, theta);
        if (std::isnan(result)) {
            result = coalescent_log_space(n, k, theta);
        }
        return result;
    } else {
        return coalescent_log_space(n, k, theta);
    }
}

template <typename RealType>
auto coalescent(const unsigned n, const unsigned k_snp, const unsigned k_indel,
                const RealType theta_snp, const RealType theta_indel)
{
    const auto theta = theta_snp + theta_indel;
    const auto k_tot = k_snp + k_indel;
    auto result = coalescent(n, k_tot, theta);
    result += k_snp * std::log(theta_snp / theta);
    result += k_indel * std::log(theta_indel / theta);
    result += maths::log_binomial_coefficient<RealType>(k_tot, k_snp);
    return result;
}

template <typename RealType, std::size_t N>
auto coalescent(const unsigned n, 
                const std::array<unsigned, N>& k,
                const std::array<RealType, N>& theta)
{
    const auto theta_tot = std::accumulate(std::cbegin(theta), std::cend(theta), RealType {0});
    const auto k_tot = std::accumulate(std::cbegin(k), std::cend(k), 0u);
    auto result = coalescent(n, k_tot, theta_tot);
    for (std::size_t i {0}; i < N; ++i) {
        result += k[i] * std::log(theta[i] / theta_tot);
    }
    result += maths::log_multinomial_coefficient<RealType>(k);
    return result;
}

} // namespace

CoalescentModel::LogProbability CoalescentModel::evaluate(const SegregatingSiteCounts& counts) const
{
    if (counts.repeat_indels + counts.complex_indels == 0) {
        return evaluate_no_indels(counts.snps, counts.haplotypes);
    } else {
        const auto indel_heterozygosities = calculate_buffered_indel_heterozygosities();
        const auto repeat_indel_heterozygosity = maths::round_sf(indel_heterozygosities.second, 6);
        const auto complex_indel_heterozygosity = maths::round_sf(indel_heterozygosities.first, 6);
        SegregatingSiteCountsWithIndelHeterozygosities sites {counts, repeat_indel_heterozygosity, complex_indel_heterozygosity};
        auto itr = k_indel_pos_result_cache_.find(sites);
        if (itr != std::cend(k_indel_pos_result_cache_)) {
            return itr->second;
        }
        const std::array<unsigned, 3> site_counts {counts.snps, counts.repeat_indels, counts.complex_indels};
        const std::array<double, 3> site_heterozygosities {params_.snp_heterozygosity, repeat_indel_heterozygosity, complex_indel_heterozygosity};
        const auto result = coalescent(counts.haplotypes, site_counts, site_heterozygosities);
        k_indel_pos_result_cache_.emplace(sites, result);
        return result;
    }
}

CoalescentModel::LogProbability CoalescentModel::evaluate_no_indels(const unsigned k_snp, const unsigned n) const
{
    if (k_indel_zero_result_cache_.size() > n) {
        if (k_indel_zero_result_cache_[n].size() > k_snp) {
            auto& result = k_indel_zero_result_cache_[n][k_snp];
            if (!result) {
                result = coalescent(n, k_snp, 0, params_.snp_heterozygosity, params_.indel_heterozygosity);
            }
            return *result;
        } else {
            k_indel_zero_result_cache_[n].resize(k_snp + 1, boost::none);
        }
    } else {
        k_indel_zero_result_cache_.resize(n + 1);
        k_indel_zero_result_cache_[n].assign(k_snp + 1, boost::none);
    }
    const auto result = coalescent(n, k_snp, 0, params_.snp_heterozygosity, params_.indel_heterozygosity);
    k_indel_zero_result_cache_[n][k_snp] = result;
    return result;
}

void CoalescentModel::fill_site_buffer(const Haplotype& haplotype) const
{
    assert(site_buffer2_.empty());
    site_buffer1_.clear();
    if (caching_ == CachingStrategy::address) {
        fill_site_buffer_from_address_cache(haplotype);
    } else {
        fill_site_buffer_from_value_cache(haplotype);
    }
    site_buffer1_ = std::move(site_buffer2_);
    site_buffer2_.clear();
}

void CoalescentModel::fill_site_buffer_uncached(const Haplotype& haplotype) const
{
    // Although we won't retrieve from the cache, we need to make sure all the variants
    // stay in existence as we populate the buffers by reference.
    auto itr = difference_value_cache_.find(reference_);
    if (itr == std::cend(difference_value_cache_)) {
        itr = difference_value_cache_.emplace(reference_, haplotype.difference(reference_)).first;
    } else {
        itr->second = haplotype.difference(reference_);
    }
    std::set_union(std::begin(site_buffer1_), std::end(site_buffer1_),
                   std::cbegin(itr->second), std::cend(itr->second),
                   std::back_inserter(site_buffer2_));
}

void CoalescentModel::fill_site_buffer_from_value_cache(const Haplotype& haplotype) const
{
    auto itr = difference_value_cache_.find(haplotype);
    if (itr == std::cend(difference_value_cache_)) {
        itr = difference_value_cache_.emplace(haplotype, haplotype.difference(reference_)).first;
    }
    std::set_union(std::begin(site_buffer1_), std::end(site_buffer1_),
                   std::cbegin(itr->second), std::cend(itr->second),
                   std::back_inserter(site_buffer2_));
}

void CoalescentModel::fill_site_buffer_from_address_cache(const Haplotype& haplotype) const
{
    auto itr = difference_address_cache_.find(std::addressof(haplotype));
    if (itr == std::cend(difference_address_cache_)) {
        itr = difference_address_cache_.emplace(std::addressof(haplotype), haplotype.difference(reference_)).first;
    }
    std::set_union(std::begin(site_buffer1_), std::end(site_buffer1_),
                   std::cbegin(itr->second), std::cend(itr->second),
                   std::back_inserter(site_buffer2_));
}

CoalescentModel::SegregatingSiteCounts CoalescentModel::count_segregating_sites(const Haplotype& haplotype) const
{
    fill_site_buffer(haplotype);
    return count_segregating_sites_in_buffer(1);
}

CoalescentModel::SegregatingSiteCounts CoalescentModel::count_segregating_sites_in_buffer(const unsigned num_haplotypes) const
{
    unsigned repeat_indels {0}, complex_indels {0};
    for (const auto& site : site_buffer1_) {
        if (is_indel(site)) {
            if (has_overlapped(reference_repeats_, site.get())) {
                ++repeat_indels;
            } else {
                ++complex_indels;
            }
        }
    }
    const auto num_indels = repeat_indels + complex_indels;
    const auto num_snps = static_cast<unsigned>(site_buffer1_.size()) - num_indels;
    return {num_haplotypes + 1, num_snps, repeat_indels, complex_indels};
}

std::pair<double, double> CoalescentModel::calculate_buffered_indel_heterozygosities() const
{
    boost::optional<double> min_heterozygosity {}, max_heterozygosity {};
    for (const auto& site : site_buffer1_) {
        if (is_indel(site)) {
            auto site_heterozygosity = calculate_heterozygosity(site);
            if (min_heterozygosity) {
                min_heterozygosity = std::min(*min_heterozygosity, site_heterozygosity);
            } else {
                min_heterozygosity = site_heterozygosity;
            }
            if (max_heterozygosity) {
                max_heterozygosity = std::max(*max_heterozygosity, site_heterozygosity);
            } else {
                max_heterozygosity = site_heterozygosity;
            }
        }
    }
    if (!min_heterozygosity) min_heterozygosity = params_.indel_heterozygosity;
    if (!max_heterozygosity) max_heterozygosity = params_.indel_heterozygosity;
    return std::make_pair(*min_heterozygosity, *max_heterozygosity);
}

double CoalescentModel::calculate_heterozygosity(const Variant& indel) const
{
    assert(is_indel(indel));
    const auto offset = static_cast<std::size_t>(begin_distance(reference_, indel));
    return calculate_indel_probability(indel_heterozygosity_model_, offset, indel_size(indel));
}

CoalescentProbabilityGreater::CoalescentProbabilityGreater(CoalescentModel model)
: model_ {std::move(model)}
, buffer_ {}
, cache_ {}
{
    buffer_.reserve(1);
    cache_.reserve(100);
}

bool CoalescentProbabilityGreater::operator()(const Haplotype& lhs, const Haplotype& rhs) const
{
    if (have_same_alleles(lhs, rhs)) return true;
    auto cache_itr = cache_.find(lhs);
    if (cache_itr == std::cend(cache_)) {
        buffer_.assign({lhs});
        cache_itr = cache_.emplace(lhs, model_.evaluate(buffer_)).first;
    }
    const auto lhs_probability = cache_itr->second;
    cache_itr = cache_.find(rhs);
    if (cache_itr == std::cend(cache_)) {
        buffer_.assign({rhs});
        cache_itr = cache_.emplace(rhs, model_.evaluate(buffer_)).first;
    }
    const auto rhs_probability = cache_itr->second;
    return lhs_probability > rhs_probability;
}

bool operator==(const CoalescentModel::SegregatingSiteCounts& lhs, const CoalescentModel::SegregatingSiteCounts& rhs) noexcept
{
    return lhs.haplotypes == rhs.haplotypes && lhs.snps == rhs.snps && lhs.repeat_indels == rhs.repeat_indels && lhs.complex_indels == rhs.complex_indels;
}

std::size_t CoalescentModel::SegregatingSiteCountsHash::operator()(const SegregatingSiteCounts& counts) const noexcept
{
    return boost::hash_value(std::make_tuple(counts.haplotypes, counts.snps, counts.repeat_indels, counts.complex_indels));
}

std::size_t CoalescentModel::SegregatingSiteCountsWithIndelHeterozygositiesHash::operator()(const SegregatingSiteCountsWithIndelHeterozygosities& sites) const noexcept
{
    using boost::hash_combine;
    std::size_t result {};
    hash_combine(result, SegregatingSiteCountsHash{}(sites.counts));
    hash_combine(result, std::hash<decltype(sites.repeat_heterozygosity)>()(sites.repeat_heterozygosity));
    hash_combine(result, std::hash<decltype(sites.complex_heterozygosity)>()(sites.complex_heterozygosity));
    return result;
}

bool CoalescentModel::SegregatingSiteCountsWithIndelHeterozygositiesEqual::operator()(const SegregatingSiteCountsWithIndelHeterozygosities& lhs, const SegregatingSiteCountsWithIndelHeterozygosities& rhs) const noexcept
{
    return lhs.counts == rhs.counts && lhs.repeat_heterozygosity == rhs.repeat_heterozygosity && lhs.complex_heterozygosity == rhs.complex_heterozygosity;
}

} // namespace octopus
