// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "coalescent_population_prior_model.hpp"

namespace octopus {

CoalescentPopulationPriorModel::CoalescentPopulationPriorModel(CoalescentModel segregation_model)
: segregation_model_ {std::move(segregation_model)}
, genotype_model_ {}
{}

CoalescentPopulationPriorModel::CoalescentPopulationPriorModel(CoalescentModel segregation_model, HardyWeinbergModel genotype_model)
: segregation_model_ {std::move(segregation_model)}
, genotype_model_ {std::move(genotype_model)}
{}

template <typename T>
auto sum_sizes(const std::vector<std::vector<T>>& values) noexcept
{
    return std::accumulate(std::cbegin(values), std::cend(values), std::size_t {0},
                           [] (auto curr, const auto& v) noexcept { return curr + v.size(); });
}

template <typename T>
auto sum_sizes(const std::vector<std::reference_wrapper<const std::vector<T>>>& values) noexcept
{
    return std::accumulate(std::cbegin(values), std::cend(values), std::size_t {0},
                           [] (auto curr, const auto& v) noexcept { return curr + v.get().size(); });
}

double CoalescentPopulationPriorModel::evaluate_segregation_model(const std::vector<std::vector<unsigned>>& indices) const
{
    if (indices.size() == 1) {
        return segregation_model_.evaluate(indices.front());
    }
    const auto num_indices = sum_sizes(indices);
    index_buffer_.resize(num_indices);
    auto itr = std::begin(index_buffer_);
    for (const auto& i : indices) {
        itr = std::copy(std::cbegin(i), std::cend(i), itr);
    }
    return segregation_model_.evaluate(index_buffer_);
}

double CoalescentPopulationPriorModel::evaluate_segregation_model(const std::vector<GenotypeIndiceVectorReference>& indices) const
{
    if (indices.size() == 1) {
        return segregation_model_.evaluate(indices.front().get());
    }
    const auto num_indices = sum_sizes(indices);
    index_buffer_.resize(num_indices);
    auto itr = std::begin(index_buffer_);
    for (const auto& i : indices) {
        itr = std::copy(std::cbegin(i.get()), std::cend(i.get()), itr);
    }
    return segregation_model_.evaluate(index_buffer_);
}

} // namespace
