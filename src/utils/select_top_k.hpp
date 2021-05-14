// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef select_top_k_hpp
#define select_top_k_hpp

#include <vector>
#include <algorithm>
#include <iterator>
#include <queue>
#include <cstddef>
#include <cmath>
#include <utility>
#include <cassert>

namespace octopus {

template <typename Index = std::size_t>
using IndexTuple = std::vector<Index>;
template <typename Index>
using IndexTupleVector = std::vector<IndexTuple<Index>>;

namespace detail {

template <typename T>
auto index_cref(const std::vector<T>& values)
{
    std::vector<std::pair<std::reference_wrapper<const T>, std::size_t>> result {};
    result.reserve(values.size());
    for (std::size_t idx {0}; idx < values.size(); ++idx) {
        result.emplace_back(values[idx], idx);
    }
    return result;
}

template <typename _, typename T>
auto copy_k_second(const std::vector<std::pair<_, T>>& pairs, std::size_t k)
{
    k = std::min(k, pairs.size());
    std::vector<T> result(k);
    const static auto get_second = [] (const auto& p) { return p.second; };
    std::transform(std::cbegin(pairs), std::next(std::cbegin(pairs), k), std::begin(result), get_second);
    return result;
}

} // namespace detail

template <typename Index, typename T, typename BinaryPredicate = std::greater<T>>
std::vector<Index>
select_top_k_indices(const std::vector<T>& values, const std::size_t k,
                     const bool sorted = true,
                     const BinaryPredicate comp = std::greater<T> {})
{
    if (!sorted && values.size() <= k) {
        std::vector<Index> result(values.size());
        std::iota(std::begin(result), std::end(result), 0);
        return result;
    } else {
        auto indexed_values = detail::index_cref(values);
        const auto ref_comp = [&comp] (const auto& lhs, const auto& rhs) { return comp(lhs.first.get(), rhs.first.get()); };
        if (k < values.size()) {
            const auto kth = std::next(std::begin(indexed_values), k);
            if (sorted) {
                std::partial_sort(std::begin(indexed_values), kth, std::end(indexed_values), ref_comp);
            } else {
                std::nth_element(std::begin(indexed_values), kth, std::end(indexed_values), ref_comp);
            }
        } else {
            std::sort(std::begin(indexed_values), std::end(indexed_values), ref_comp);
        }
        return detail::copy_k_second(indexed_values, k);
    }
}

namespace detail {

template <typename Index, typename T>
auto index(const std::vector<T>& values)
{
    std::vector<std::pair<T, Index>> result(values.size());
    for (std::size_t i {0}; i < values.size(); ++i) {
        result[i] = std::make_pair(values[i], i);
    }
    return result;
}

template <typename Index, typename T>
auto index_and_sort(const std::vector<T>& values, const std::size_t k)
{
    auto result = index<Index>(values);
    const auto kth = std::next(std::begin(result), std::min(k, result.size()));
    std::partial_sort(std::begin(result), kth, std::end(result), std::greater<> {});
    result.erase(kth, std::end(result));
    return result;
}

template <typename Index, typename T>
struct IndexTupleScorePair
{
    IndexTuple<Index> indices;
    T score;
};

template <typename Index, typename T>
bool operator>(const IndexTupleScorePair<Index, T>& lhs, const IndexTupleScorePair<Index, T>& rhs) noexcept
{
    return lhs.score > rhs.score;
}

template <typename Index, typename T>
using IndexTupleScorePairVector = std::vector<IndexTupleScorePair<Index, T>>;

template <typename Index>
using IndexPair = std::pair<Index, Index>;

template <typename Index, typename T>
std::vector<IndexPair<Index>>
find_k_max_pairs(const std::vector<T>& lhs, const std::vector<T>& rhs, const std::size_t k)
{
    // Modified from https://leetcode.com/problems/find-k-pairs-with-smallest-sums/discuss/84607/Clean-16ms-C++-O(N)-Space-O(KlogN)-Time-Solution-using-Priority-queue
    std::vector<IndexPair<Index>> result {};
    if (lhs.empty() || rhs.empty() || k == 0)
        return result;
    auto cmp = [&lhs, &rhs] (const IndexPair<Index>& a, const IndexPair<Index>& b) {
        return lhs[a.first] + rhs[a.second] < lhs[b.first] + rhs[b.second]; };
    std::vector<IndexPair<Index>> heap_storage {};
    heap_storage.reserve(k);
    std::priority_queue<IndexPair<Index>, decltype(heap_storage), decltype(cmp)> max_heap {cmp, std::move(heap_storage)};
    max_heap.emplace(0, 0);
    for (std::size_t i {0}; i < k && !max_heap.empty(); ++i) {
        const auto idx_pair = max_heap.top(); max_heap.pop();
        result.push_back(idx_pair);
        if (idx_pair.first + 1 < lhs.size())
            max_heap.emplace(idx_pair.first + 1, idx_pair.second);
        if (idx_pair.first == 0 && idx_pair.second + 1 < rhs.size())
            max_heap.emplace(idx_pair.first, idx_pair.second + 1);
    }
    return result;
}

template <typename Index, typename T>
void join(const std::vector<std::pair<T, Index>>& values,
          IndexTupleScorePairVector<Index, T>& result,
          const std::size_t k,
          IndexTupleScorePairVector<Index, T>& buffer)
{
    const auto n = std::min(k, values.size());
    if (n == 0) return;
    if (result.empty()) {
        std::transform(std::cbegin(values), std::next(std::cbegin(values), n),
                       std::back_inserter(result), [=] (const auto& p) -> IndexTupleScorePair<Index, T> {
            return {{p.second}, p.first};
        });
    } else {
        if (values.size() > 1) {
            std::vector<T> vals1(values.size()), vals2(result.size());
            std::transform(std::cbegin(values), std::cend(values), std::begin(vals1), [] (const auto& p) { return p.first; });
            std::transform(std::cbegin(result), std::cend(result), std::begin(vals2), [] (const auto& p) { return p.score; });
            const auto max_pairs = find_k_max_pairs<Index>(vals1, vals2, std::min(k, n * result.size()));
            buffer.clear();
            buffer.reserve(max_pairs.size());
            for (const auto& p : max_pairs) {
                buffer.push_back(result[p.second]);
                buffer.back().indices.push_back(values[p.first].second);
                buffer.back().score += values[p.first].first;
            }
            result = std::move(buffer);
            std::sort(std::begin(result), std::end(result), std::greater<> {});
        } else {
            assert(values.size() == 1);
            for (auto& t : result) {
                t.indices.push_back(values[0].second);
                t.score += values[0].first;
            }
        }
    }
}

} // namespace detail

template <typename Index, typename T>
IndexTupleVector<Index>
select_top_k_tuples(const std::vector<std::vector<T>>& values, const std::size_t k)
{
    // Implements a variant of the top-K selection algorithm
    // See Henderson & Eliassi-Rad (http://eliassi.org/papers/henderson-llnltr09.pdf)
    // for comparison.
    // Returns results in descending score order
    detail::IndexTupleScorePairVector<Index, T> joins {}, buffer {};
    joins.reserve(k);
    for (const auto& v : values) {
        detail::join(detail::index_and_sort<Index>(v, k), joins, k, buffer);
    }
    IndexTupleVector<Index> result {};
    result.reserve(k);
    for (auto& p : joins) {
        result.push_back(std::move(p.indices));
        p.indices = {};
    }
    return result;
}

} // namespace octopus

#endif
