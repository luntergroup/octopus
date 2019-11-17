// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef k_medoids_hpp
#define k_medoids_hpp

#include <vector>
#include <iterator>
#include <cstddef>
#include <functional>
#include <algorithm>
#include <numeric>
#include <random>
#include <utility>
#include <cassert>

namespace octopus {

struct KMediodsParameters
{
    enum class InitialisationMode { random, max_distance } initialisation = InitialisationMode::max_distance;
    std::size_t max_iterations = 100;
    bool allow_empty_clusters = false;
};

using Cluster = std::vector<std::size_t>;
using ClusterVector = std::vector<Cluster>;

namespace detail {

template <typename T>
using DistanceVector = std::vector<T>;
template <typename T>
using DistanceMatrix = std::vector<DistanceVector<T>>;

template <typename ForwardIterator, typename BinaryFunction>
auto make_distance_matrix(const ForwardIterator first_data_itr, const ForwardIterator last_data_itr,
                          const BinaryFunction& distance)
{
    using DistanceResultType = decltype(distance(*first_data_itr, *first_data_itr));
    const auto N = static_cast<std::size_t>(std::distance(first_data_itr, last_data_itr));
    DistanceMatrix<DistanceResultType> result(N, DistanceVector<DistanceResultType>(N));
    for (std::size_t i {0}; i < N; ++i) {
        for (auto j = i + 1; j < N; ++j) {
            result[i][j] = distance(*std::next(first_data_itr, i), *std::next(first_data_itr, j));
            result[j][i] = result[i][j];
        }
    }
    return result;
}

using MediodVector = std::vector<std::size_t>;

inline auto initialise_mediods_random(const std::size_t k, const std::size_t N)
{
    MediodVector result(N);
    std::iota(std::begin(result), std::end(result), 0);
    if (k < N) {
        static std::mt19937 generator {42};
        std::shuffle(std::begin(result), std::end(result), generator);
        result.resize(k);
    }
    return result;
}

template <typename D>
std::pair<std::size_t, std::size_t>
max_distance(const DistanceMatrix<D>& distances)
{
    std::size_t max_i {0}, max_j {1};
    auto max_distance = distances[0][1];
    for (std::size_t i {0}; i < distances.size(); ++i) {
        for (auto j = i + 1; j < distances.size(); ++j) {
            if (distances[i][j] > max_distance) {
                max_i = i; max_j = j;
                max_distance = distances[i][j];
            }
        }
    }
    return {max_i, max_j};
}

template <typename D>
std::size_t
maxmin_distance_point(const std::vector<std::size_t>& medoids, const DistanceMatrix<D>& distances)
{
    std::size_t result {};
    auto max_distance = distances[0][0];
    for (std::size_t i {0}; i < distances.size(); ++i) {
        if (std::find(std::cbegin(medoids), std::cend(medoids), i) ==  std::cend(medoids)) {
            const auto medoid_less = [&] (auto lhs, auto rhs) { return distances[lhs][i] < distances[rhs][i]; };
            const auto min_medoid = *std::min_element(std::cbegin(medoids), std::cend(medoids), medoid_less);
            if (distances[min_medoid][i] > max_distance) {
                result = i;
                max_distance = distances[min_medoid][i];
            }
        }
    }
    return result;
}

template <typename ForwardIterator, typename D>
auto
initialise_mediods_max_distance(ForwardIterator first, ForwardIterator last,
                                const DistanceMatrix<D>& distances,
                                const std::size_t k)
{
    assert(k > 1);
    MediodVector result(2);
    result.reserve(k);
    std::tie(result[0], result[1]) = max_distance(distances);
    std::generate_n(std::back_inserter(result), k - 2, [&] () { return maxmin_distance_point(result, distances); });
    return result;
}

template <typename ForwardIterator, typename D>
auto
initialise_mediods(ForwardIterator first, ForwardIterator last, const std::size_t k,
                   const DistanceMatrix<D>& distances,
                   const KMediodsParameters& params)
{
    const auto N = static_cast<std::size_t>(std::distance(first, last));
    if (params.initialisation == KMediodsParameters::InitialisationMode::random) {
        return initialise_mediods_random(k, N);
    } else {
        return initialise_mediods_max_distance(first, last, distances, k);
    }
}

template <typename D>
auto choose_medoid(const std::size_t point_idx,
                   const MediodVector& medoids,
                   const DistanceMatrix<D>& distances)
{
    auto min_medoids_itr = std::find(std::cbegin(medoids), std::cend(medoids), point_idx);
    // If the point is already a medoid then use this as must be minimum distance.
    // Also ensures no cluster will be empty.
    if (min_medoids_itr == std::cend(medoids)) {
        const auto distance_less = [&] (const auto& lhs, const auto& rhs) {
            return distances[point_idx][lhs] < distances[point_idx][rhs]; };
        min_medoids_itr = std::min_element(std::cbegin(medoids), std::cend(medoids), distance_less);
    }
    return static_cast<std::size_t>(std::distance(std::cbegin(medoids), min_medoids_itr));
}

template <typename D>
auto find_worst_clustered_point(const ClusterVector& clusters, const MediodVector& medoids,
                                const DistanceMatrix<D>& distances)
{
    int worst_cluster {-1};
    std::size_t worst_point {};
    D max_distance {};
    for (std::size_t cluster_idx {0}; cluster_idx < clusters.size(); ++cluster_idx) {
        for (const auto point : clusters[cluster_idx]) {
            if (worst_cluster == -1 || max_distance < distances[point][medoids[cluster_idx]]) {
                worst_cluster = cluster_idx;
                worst_point = point;
                max_distance = distances[point][medoids[cluster_idx]];
            }
        }
    }
    return std::make_pair(static_cast<std::size_t>(worst_cluster), worst_point);
}

template <typename D>
void assign_empty_clusters(ClusterVector& clusters,
                           const MediodVector& medoids,
                           const DistanceMatrix<D>& distances)
{
    for (auto& cluster : clusters) {
        if (cluster.empty()) {
            auto p = find_worst_clustered_point(clusters, medoids, distances);
            cluster.push_back(p.second);
            clusters[p.first].erase(std::find(std::cbegin(clusters[p.first]), std::cend(clusters[p.first]), p.second));
        }
    }
}

template <typename D>
void initialise_clusters(ClusterVector& clusters,
                         const MediodVector& medoids,
                         const DistanceMatrix<D>& distances,
                         const KMediodsParameters& params)
{
    clusters.resize(medoids.size());
    const auto data_size = distances.size();
    for (auto& cluster : clusters) cluster.reserve(data_size - 1);
    for (std::size_t point_idx {0}; point_idx < data_size; ++point_idx) {
        auto cluster_idx = choose_medoid(point_idx, medoids, distances);
        clusters[cluster_idx].push_back(point_idx);
    }
    if (!params.allow_empty_clusters) {
        assign_empty_clusters(clusters, medoids, distances);
    }
}

template <typename D>
bool update_clusters(ClusterVector& clusters,
                     const MediodVector& medoids,
                     const DistanceMatrix<D>& distances,
                     const KMediodsParameters& params)
{
    assert(medoids.size() == clusters.size());
    const auto data_size = distances.size();
    std::vector<bool> assigned(data_size);
    bool assignments_changed {false};
    std::vector<std::size_t> buffer {};
    for (std::size_t cluster_idx {0}; cluster_idx < clusters.size(); ++cluster_idx) {
        buffer.reserve(clusters[cluster_idx].size());
        for (auto point_idx : clusters[cluster_idx]) {
            if (!assigned[point_idx]) {
                const auto new_cluster_idx = choose_medoid(point_idx, medoids, distances);
                if (new_cluster_idx != cluster_idx) {
                    clusters[new_cluster_idx].push_back(point_idx);
                    assignments_changed = true;
                } else {
                    buffer.push_back(point_idx);
                }
                assigned[point_idx] = true;
            } else {
                buffer.push_back(point_idx);
            }
        }
        clusters[cluster_idx] = std::move(buffer);
        buffer.clear();
    }
    if (!params.allow_empty_clusters && assignments_changed) {
        assign_empty_clusters(clusters, medoids, distances);
    }
    return assignments_changed;
}

template <typename D>
auto sum(const DistanceMatrix<D>& distances, const std::vector<std::size_t>& cluster, const std::size_t medoid)
{
    return std::accumulate(std::cbegin(cluster), std::cend(cluster), D {},
                           [&] (D curr, std::size_t point) { return curr + distances[point][medoid]; });
}

template <typename D>
auto choose_medoid(const std::vector<std::size_t>& cluster,
                   const DistanceMatrix<D>& distances)
{
    std::size_t result {0};
    D min_distance_sum {};
    for (auto idx : cluster) {
        auto distance_sum = sum(distances, cluster, idx);
        if (idx == 0 || distance_sum < min_distance_sum) {
            result = idx;
            min_distance_sum = std::move(distance_sum);
        }
    }
    return result;
}

template <typename D>
void update_mediods(MediodVector& medoids,
                    const ClusterVector& clusters,
                    const DistanceMatrix<D>& distances)
{
    assert(medoids.size() == clusters.size());
    for (std::size_t cluster_idx {0}; cluster_idx < clusters.size(); ++cluster_idx) {
        medoids[cluster_idx] = choose_medoid(clusters[cluster_idx], distances);
    }
}

} // namespace detail

template <typename ForwardIterator, typename BinaryFunction>
auto
k_medoids(const ForwardIterator first_data_itr, const ForwardIterator last_data_itr,
          const std::size_t k,
          ClusterVector& clusters,
          BinaryFunction distance,
          const KMediodsParameters params = KMediodsParameters {})
{
    const auto distances = detail::make_distance_matrix(first_data_itr, last_data_itr, distance);
    auto medoids = detail::initialise_mediods(first_data_itr, last_data_itr, k, distances, params);
    if (static_cast<std::size_t>(std::distance(first_data_itr, last_data_itr)) <= k) {
        clusters.reserve(medoids.size());
        for (auto medoid : medoids) clusters.push_back({medoid});
    }
    detail::initialise_clusters(clusters, medoids, distances, params);
    for (std::size_t n {0}; n < params.max_iterations; ++n) {
        detail::update_mediods(medoids, clusters, distances);
        auto assignments_changed = detail::update_clusters(clusters, medoids, distances, params);
        if (!assignments_changed) break;
    }
    for (auto& cluster : clusters) cluster.shrink_to_fit();
    return medoids;
}

template <typename ForwardIterator>
auto
k_medoids(ForwardIterator first, ForwardIterator last, const std::size_t k,
          ClusterVector& clusters,
          const KMediodsParameters params = KMediodsParameters {})
{
    return k_medoids(first, last, k, clusters, [] (const auto& lhs, const auto& rhs) { return (lhs - rhs) * (lhs - rhs); }, params);
}

template <typename ForwardIterator, typename BinaryFunction>
auto
k_medoids(ForwardIterator first, ForwardIterator last, const std::size_t k,
          BinaryFunction distance,
          const KMediodsParameters params = KMediodsParameters {})
{
    ClusterVector clusters {};
    return k_medoids(first, last, k, clusters, std::move(distance), params);
}

template <typename ForwardIterator>
auto
k_medoids(ForwardIterator first, ForwardIterator last,
          const std::size_t k,
          const KMediodsParameters params = KMediodsParameters {})
{
    return k_medoids(first, last, k, [] (const auto& lhs, const auto& rhs) { return (lhs - rhs) * (lhs - rhs); }, params);
}

template <typename Range>
auto
k_medoids(const Range& values, const std::size_t k, const KMediodsParameters params = KMediodsParameters {})
{
    return k_medoids(std::cbegin(values), std::cend(values), k, params);
}

template <typename Range, typename BinaryFunction>
auto
k_medoids(const Range& values, const std::size_t k,
          BinaryFunction distance,
          const KMediodsParameters params = KMediodsParameters {})
{
    return k_medoids(std::cbegin(values), std::cend(values), k, std::move(distance), params);
}

template <typename Range>
auto
k_medoids(const Range& values, const std::size_t k,
          ClusterVector& clusters,
          const KMediodsParameters params = KMediodsParameters {})
{
    return k_medoids(std::cbegin(values), std::cend(values), k, clusters, params);
}

template <typename Range, typename BinaryFunction>
auto
k_medoids(const Range& values, const std::size_t k,
          ClusterVector& clusters,
          BinaryFunction distance,
          const KMediodsParameters params = KMediodsParameters {})
{
    return k_medoids(std::cbegin(values), std::cend(values), k, clusters, std::move(distance), params);
}

} // namespace octopus

#endif
