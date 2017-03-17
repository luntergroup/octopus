// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "coalescent_model.hpp"

#include <memory>
#include <cmath>
#include <complex>
#include <numeric>
#include <stdexcept>

#include <boost/math/special_functions/binomial.hpp>

#include "tandem/tandem.hpp"
#include "utils/maths.hpp"

namespace octopus {

auto percent_of_bases_in_repeat(const Haplotype& haplotype)
{
    const auto repeats = tandem::extract_exact_tandem_repeats(haplotype.sequence(), 1, 6);
    if (repeats.empty()) return 0.0;
    std::vector<unsigned> repeat_counts(sequence_size(haplotype), 0);
    for (const auto& repeat : repeats) {
        const auto itr1 = std::next(std::begin(repeat_counts), repeat.pos);
        const auto itr2 = std::next(itr1, repeat.length);
        std::transform(itr1, itr2, itr1, [] (const auto c) { return c + 1; });
    }
    const auto c = std::count_if(std::cbegin(repeat_counts), std::cend(repeat_counts),
                                 [] (const auto c) { return c > 0; });
    return static_cast<double>(c) / repeat_counts.size();
}

auto calculate_base_indel_heterozygosities(const Haplotype& haplotype,
                                           const double base_indel_heterozygosity)
{
    std::vector<double> result(sequence_size(haplotype), base_indel_heterozygosity);
    const auto repeats = tandem::extract_exact_tandem_repeats(haplotype.sequence(), 1, 3);
    for (const auto& repeat : repeats) {
        const auto itr1 = std::next(std::begin(result), repeat.pos);
        const auto itr2 = std::next(itr1, repeat.length);
        const auto n = repeat.length / repeat.period;
        // TODO: implement a proper model for this
        const auto t = std::min(base_indel_heterozygosity * std::pow(n, 2.6), 1.0);
        std::transform(itr1, itr2, itr1, [t] (const auto h) { return std::max(h, t); });
    }
    return result;
}

CoalescentModel::CoalescentModel(Haplotype reference,
                                 Parameters params,
                                 std::size_t num_haplotyes_hint,
                                 CachingStrategy caching)
: reference_ {std::move(reference)}
, reference_base_indel_heterozygosities_ {}
, params_ {params}
, caching_ {caching}
{
    if (params_.snp_heterozygosity <= 0 || params_.indel_heterozygosity <= 0) {
        throw std::domain_error {"CoalescentModel: snp and indel heterozygosity must be > 0"};
    }
    reference_base_indel_heterozygosities_ = calculate_base_indel_heterozygosities(reference_, params_.indel_heterozygosity);
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
    result_cache_.reserve(num_haplotyes_hint);
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
    const auto l = [](const auto& lhs, const auto& rhs) { return lhs.real() < rhs.real(); };
    const auto max = *std::max_element(first, last, l);
    return max + std::log(std::accumulate(first, last, ComplexType {},
                                          [max](const auto curr, const auto x) {
                                              return curr + std::exp(x - max);
                                          }));
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
        // indel heterozygosity is default in this case
        const auto itr = result_cache_.find(t);
        if (itr != std::cend(result_cache_)) return itr->second;
    }
    auto indel_heterozygosity = params_.indel_heterozygosity;
    if (k_indel > 0) {
        for (const auto& site : site_buffer1_) {
            if (is_indel(site)) {
                const auto offset = begin_distance(reference_, site.get());
                auto itr = std::next(std::cbegin(reference_base_indel_heterozygosities_), offset);
                using S = Variant::MappingDomain::Size;
                itr = std::max_element(itr, std::next(itr, std::max(S {1}, region_size(site.get()))));
                indel_heterozygosity = std::max(*itr, indel_heterozygosity);
            }
        }
    }
    const auto result = coalescent(n, k_snp, k_indel, params_.snp_heterozygosity, indel_heterozygosity);
    if (k_indel > 0) {
        result_cache_.emplace(t, result);
    }
    return result;
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
    if (first_haplotype_) {
        fill_site_buffer_from_flat_address_cache(haplotype);
    } else {
        fill_site_buffer_from_sparse_address_cache(haplotype);
    }
}

void CoalescentModel::fill_site_buffer_from_flat_address_cache(const Haplotype& haplotype) const noexcept
{
    assert(first_haplotype_);
    assert(*first_haplotype_ <= std::addressof(haplotype) && std::addressof(haplotype) <= *last_haplotype_);
    const auto haplotype_index = static_cast<std::size_t>(std::distance(*first_haplotype_, std::addressof(haplotype)));
    assert(haplotype_index < flat_address_cache_.size());
    auto& difference = flat_address_cache_[haplotype_index];
    if (!difference) {
        difference = haplotype.difference(reference_);
    }
    std::set_union(std::begin(site_buffer1_), std::end(site_buffer1_),
                   std::cbegin(*difference), std::cend(*difference),
                   std::back_inserter(site_buffer2_));
}

void CoalescentModel::fill_site_buffer_from_sparse_address_cache(const Haplotype& haplotype) const
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

void CoalescentModel::do_prime(const std::vector<const Haplotype*>& addresses) const
{
    if (!addresses.empty()) {
        const auto p = std::minmax_element(std::cbegin(addresses), std::cend(addresses));
        const auto num_spaces = static_cast<std::size_t>(std::distance(*p.first, *p.second)) + 1;
        if (num_spaces <= max_flat_spaces_) {
            first_haplotype_ = *p.first;
            last_haplotype_  = *p.second;
            flat_address_cache_.assign(num_spaces, boost::none);
        }
    }
}

void CoalescentModel::unprime() const noexcept
{
    first_haplotype_ = boost::none;
    last_haplotype_  = boost::none;
    flat_address_cache_.clear();
    flat_address_cache_.shrink_to_fit();
}

} // namespace octopus
