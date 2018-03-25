// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "coalescent_model.hpp"

#include <memory>
#include <cmath>
#include <complex>
#include <numeric>
#include <stdexcept>

#include <boost/math/special_functions/binomial.hpp>

#include "utils/maths.hpp"

namespace octopus {

CoalescentModel::CoalescentModel(Haplotype reference, Parameters params,
                                 std::size_t num_haplotyes_hint, CachingStrategy caching)
: reference_ {std::move(reference)}
, indel_heterozygosity_model_ {make_indel_model(reference_, {params.indel_heterozygosity})}
, params_ {params}
, haplotypes_ {}
, caching_ {caching}
, index_cache_ {}
, index_flag_buffer_ {}
, k_indel_zero_result_cache_ {2 * num_haplotyes_hint, std::vector<boost::optional<double>> {}}
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

void CoalescentModel::prime(std::vector<Haplotype> haplotypes)
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

double CoalescentModel::evaluate(const Haplotype& haplotype) const
{
    return evaluate(count_segregating_sites(haplotype));
}

double CoalescentModel::evaluate(const std::vector<unsigned>& haplotype_indices) const
{
    return evaluate(count_segregating_sites(haplotype_indices));
}

namespace {

auto powm1(const unsigned i) noexcept // std::pow(-1, i)
{
    return (i % 2 == 0) ? 1 : -1;
}

auto binom(const unsigned n, const unsigned k)
{
    return boost::math::binomial_coefficient<double>(n, k);
}

auto log_binom(const unsigned n, const unsigned k)
{
    using maths::log_factorial;
    return log_factorial<double>(n) - (log_factorial<double>(k) + log_factorial<double>(n - k));
}

auto coalescent_real_space(const unsigned n, const unsigned k, const double theta)
{
    double result {0};
    for (unsigned i {2}; i <= n; ++i) {
        result += powm1(i) * binom(n - 1, i - 1) * ((i - 1) / (theta + i - 1)) * std::pow(theta / (theta + i - 1), k);
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

auto coalescent_log_space(const unsigned n, const unsigned k, const double theta)
{
    std::vector<std::complex<double>> tmp(n - 1, std::log(std::complex<double> {-1}));
    for (unsigned i {2}; i <= n; ++i) {
        auto& cur = tmp[i - 2];
        cur *= i;
        cur += log_binom(n - 1, i - 1);
        cur += std::log((i - 1) / (theta + i - 1));
        cur += k * std::log(theta / (theta + i - 1));
    }
    return complex_log_sum_exp(tmp).real();
}

auto coalescent(const unsigned n, const unsigned k, const double theta)
{
    if (k <= 80) {
        return coalescent_real_space(n, k, theta);
    } else {
        return coalescent_log_space(n, k, theta);
    }
}

auto coalescent(const unsigned n, const unsigned k_snp, const unsigned k_indel,
                const double theta_snp, const double theta_indel)
{
    const auto theta = theta_snp + theta_indel;
    const auto k_tot = k_snp + k_indel;
    auto result = coalescent(n, k_tot, theta);
    result += k_snp * std::log(theta_snp / theta);
    result += k_indel * std::log(theta_indel / theta);
    result += log_binom(k_tot, k_snp);
    return result;
}

} // namespace

double CoalescentModel::evaluate(const SiteCountTuple& t) const
{
    unsigned k_snp, k_indel, n;
    std::tie(k_snp, k_indel, n) = t;
    if (k_indel == 0) {
        return evaluate(k_snp, n);
    } else {
        return evaluate(k_snp, k_indel, n);
    }
}

double CoalescentModel::evaluate(const unsigned k_snp, const unsigned n) const
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

double CoalescentModel::evaluate(const unsigned k_snp, const unsigned k_indel, const unsigned n) const
{
    const auto indel_heterozygosity = calculate_buffered_indel_heterozygosity();
    const auto t = std::make_tuple(k_snp, k_indel, n, maths::round_sf(indel_heterozygosity, 6));
    auto itr = k_indel_pos_result_cache_.find(t);
    if (itr != std::cend(k_indel_pos_result_cache_)) {
        return itr->second;
    }
    const auto result = coalescent(n, k_snp, k_indel, params_.snp_heterozygosity, indel_heterozygosity);
    k_indel_pos_result_cache_.emplace(t, result);
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

void CoalescentModel::fill_site_buffer(const std::vector<unsigned>& haplotype_indices) const
{
    site_buffer1_.clear();
    std::fill(std::begin(index_flag_buffer_), std::end(index_flag_buffer_), false);
    unsigned num_unique_nonempty_indices {0};
    auto middle = std::begin(site_buffer1_);
    for (auto index : haplotype_indices) {
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

void CoalescentModel::fill_site_buffer_from_value_cache(const Haplotype& haplotype) const
{
    auto itr = difference_value_cache_.find(haplotype);
    if (itr == std::cend(difference_value_cache_)) {
        itr = difference_value_cache_.emplace(std::piecewise_construct,
                                              std::forward_as_tuple(haplotype),
                                              std::forward_as_tuple(haplotype.difference(reference_))).first;
    }
    std::set_union(std::begin(site_buffer1_), std::end(site_buffer1_),
                   std::cbegin(itr->second), std::cend(itr->second),
                   std::back_inserter(site_buffer2_));
}

void CoalescentModel::fill_site_buffer_from_address_cache(const Haplotype& haplotype) const
{
    auto itr = difference_address_cache_.find(std::addressof(haplotype));
    if (itr == std::cend(difference_address_cache_)) {
        itr = difference_address_cache_.emplace(std::piecewise_construct,
                                                std::forward_as_tuple(std::addressof(haplotype)),
                                                std::forward_as_tuple(haplotype.difference(reference_))).first;
    }
    std::set_union(std::begin(site_buffer1_), std::end(site_buffer1_),
                   std::cbegin(itr->second), std::cend(itr->second),
                   std::back_inserter(site_buffer2_));
}

CoalescentModel::SiteCountTuple CoalescentModel::count_segregating_sites(const Haplotype& haplotype) const
{
    fill_site_buffer(haplotype);
    return count_segregating_sites_in_buffer(1);
}

CoalescentModel::SiteCountTuple CoalescentModel::count_segregating_sites_in_buffer(const unsigned num_haplotypes) const
{
    const auto num_indels = std::count_if(std::cbegin(site_buffer1_), std::cend(site_buffer1_),
                                          [] (const auto& v) noexcept { return is_indel(v); });
    return std::make_tuple(site_buffer1_.size() - num_indels, num_indels, num_haplotypes + 1);
}

double CoalescentModel::calculate_buffered_indel_heterozygosity() const
{
    boost::optional<double> result {};
    for (const auto& site : site_buffer1_) {
        if (is_indel(site)) {
            auto site_heterozygosity = calculate_heterozygosity(site);
            if (result) {
                result = std::max(*result, site_heterozygosity);
            } else {
                result = site_heterozygosity;
            }
        }
    }
    return result ? *result : params_.indel_heterozygosity;
}

double CoalescentModel::calculate_heterozygosity(const Variant& indel) const
{
    assert(is_indel(indel));
    const auto offset = static_cast<std::size_t>(begin_distance(reference_, indel));
    const auto indel_length = indel_size(indel);
    assert(offset < indel_heterozygosity_model_.gap_open.size());
    constexpr decltype(indel_length) max_indel_length {50};
    return indel_heterozygosity_model_.gap_open[offset]
           * std::pow(indel_heterozygosity_model_.gap_extend[offset], std::min(indel_length, max_indel_length) - 1);
}

} // namespace octopus
