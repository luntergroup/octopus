// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "haplotype_filter.hpp"

#include <deque>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <limits>
#include <type_traits>
#include <limits>
#include <cassert>

#include <boost/iterator/transform_iterator.hpp>

#include "config/common.hpp"
#include "core/models/haplotype_likelihood_array.hpp"
#include "utils/maths.hpp"
#include "utils/append.hpp"
#include "logging/logging.hpp"

namespace octopus {

namespace {

struct SampleIntersectTag {};
struct SamplePoolTag {};

template <typename T>
struct FilterGreater
{
    FilterGreater(const std::unordered_map<Haplotype, T>& values) : values_ {values} {}
    
    bool operator()(const Haplotype& lhs, const T rhs) const
    {
        return values_.at(lhs) > rhs;
    }
    bool operator()(const T& lhs, const Haplotype& rhs) const
    {
        return lhs > values_.at(rhs);
    }
    bool operator()(const Haplotype& lhs, const Haplotype& rhs) const
    {
        return values_.at(lhs) > values_.at(rhs);
    }
    
private:
    const std::unordered_map<Haplotype, T>& values_;
};

template <typename T>
bool are_equal(const T& lhs, const T& rhs, std::true_type) noexcept
{
    return maths::almost_equal(lhs, rhs);
}

template <typename T>
bool are_equal(const T& lhs, const T& rhs, std::false_type) noexcept
{
    return lhs == rhs;
}

template <typename T>
bool are_equal(const T& lhs, const T& rhs) noexcept
{
    return are_equal(lhs, rhs, std::is_floating_point<T> {});
}

template <typename F>
std::size_t try_filter(std::vector<Haplotype>& haplotypes,
                       const std::vector<SampleName>& samples,
                       const HaplotypeLikelihoodArray& haplotype_likelihoods,
                       const std::size_t n, std::vector<Haplotype>& result,
                       F filter, SampleIntersectTag)
{
    std::sort(std::begin(haplotypes), std::end(haplotypes));
    using T = std::result_of_t<F(const Haplotype&, decltype(samples), decltype(haplotype_likelihoods))>;
    bool first_sample {true};
    std::vector<Haplotype> new_filtered {};
    for (const auto& sample : samples) {
        std::vector<std::pair<HaplotypeReference, T>> filter_scores {};
        for (const auto& haplotype : haplotypes) {
            filter_scores.emplace_back(haplotype, filter(haplotype, sample, haplotype_likelihoods));
        }
        const auto score_less = [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; };
        const auto first_score_itr = std::begin(filter_scores);
        const auto last_score_itr  = std::end(filter_scores);
        const auto nth_score_itr = std::next(first_score_itr, n);
        std::nth_element(first_score_itr, nth_score_itr, last_score_itr, score_less);
        const auto rlast = std::make_reverse_iterator(first_score_itr);
        const auto itr = std::find_if(std::make_reverse_iterator(std::prev(nth_score_itr)), rlast,
                                     [&] (const auto& p) { return are_equal(p.second, nth_score_itr->second); });
        auto first_filtered_score_itr = nth_score_itr;
        if (itr != rlast) {
            std::sort(nth_score_itr, last_score_itr, score_less);
            first_filtered_score_itr = std::upper_bound(first_score_itr, last_score_itr, nth_score_itr->second,
                                                        [] (const auto& lhs, const auto& rhs) { return lhs > rhs.second; });
        }
        if (first_sample) {
            new_filtered.reserve(std::distance(first_filtered_score_itr, last_score_itr));
            std::transform(first_filtered_score_itr, last_score_itr, std::back_inserter(new_filtered),
                           [] (const auto& p) { return p.first; });
            std::sort(std::begin(new_filtered), std::end(new_filtered)); // was sorted by score
            first_sample = false;
        } else {
            std::vector<Haplotype> filtered_haplotypes {};
            filtered_haplotypes.reserve(std::distance(first_filtered_score_itr, last_score_itr));
            std::transform(first_filtered_score_itr, last_score_itr, std::back_inserter(filtered_haplotypes),
                           [] (const auto& p) { return p.first; });
            std::sort(std::begin(filtered_haplotypes), std::end(filtered_haplotypes));
            std::vector<Haplotype> tmp {};
            tmp.reserve(std::min(filtered_haplotypes.size(), new_filtered.size()));
            assert(std::is_sorted(std::cbegin(new_filtered), std::cend(new_filtered)));
            std::set_intersection(std::cbegin(filtered_haplotypes), std::cend(filtered_haplotypes),
                                  std::cbegin(new_filtered), std::cend(new_filtered),
                                  std::back_inserter(tmp));
            new_filtered = std::move(tmp);
        }
    }
    std::vector<Haplotype> tmp {};
    tmp.reserve(haplotypes.size() - new_filtered.size());
    assert(std::is_sorted(std::cbegin(new_filtered), std::cend(new_filtered)));
    std::set_difference(std::cbegin(haplotypes), std::cend(haplotypes),
                        std::cbegin(new_filtered), std::cend(new_filtered),
                        std::back_inserter(tmp));
    haplotypes = std::move(tmp);
    const auto num_new_filtered = new_filtered.size();
    utils::append(std::move(new_filtered), result);
    return num_new_filtered;
}

template <typename F>
std::size_t try_filter(std::vector<Haplotype>& haplotypes,
                       const std::vector<SampleName>& samples,
                       const HaplotypeLikelihoodArray& haplotype_likelihoods,
                       const std::size_t n, std::vector<Haplotype>& result,
                       F filter, SamplePoolTag)
{
    using T = std::result_of_t<F(const Haplotype&, decltype(samples), decltype(haplotype_likelihoods))>;
    std::unordered_map<Haplotype, T> filter_score {haplotypes.size()};
    for (const auto& haplotype : haplotypes) {
        filter_score.emplace(haplotype, filter(haplotype, samples, haplotype_likelihoods));
    }
    const FilterGreater<T> cmp {filter_score};
    const auto first = std::begin(haplotypes);
    const auto last  = std::end(haplotypes);
    const auto nth   = std::next(first, n);
    std::nth_element(first, nth, last, cmp);
    const auto nth_best = filter_score.at(*nth);
    const auto rlast = std::make_reverse_iterator(first);
    const auto it = std::find_if(std::make_reverse_iterator(std::prev(nth)), rlast,
                                 [nth_best, &filter_score] (const Haplotype& haplotype) {
                                     return are_equal(filter_score.at(haplotype), nth_best);
                                 });
    if (it != rlast) {
        std::sort(first, nth, cmp);
        std::sort(nth, last, cmp);
        const auto er = std::equal_range(first, last, nth_best, cmp);
        result.insert(std::end(result), std::make_move_iterator(er.second), std::make_move_iterator(last));
        const auto num_removed = std::distance(er.second, last);
        haplotypes.erase(er.second, last);
        return num_removed;
    }
    result.insert(std::end(result), std::make_move_iterator(nth), std::make_move_iterator(last));
    const auto num_removed = std::distance(nth, last);
    haplotypes.erase(nth, last);
    return num_removed;
}

template <typename F>
void force_filter(std::vector<Haplotype>& haplotypes,
                  const std::vector<SampleName>& samples,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods,
                  const std::size_t n, std::vector<Haplotype>& result,
                  F filter)
{
    using T = std::result_of_t<F(const Haplotype&, decltype(samples), decltype(haplotype_likelihoods))>;
    std::unordered_map<Haplotype, T> filter_likelihoods {haplotypes.size()};
    for (const auto& haplotype : haplotypes) {
        filter_likelihoods.emplace(haplotype, filter(haplotype, samples, haplotype_likelihoods));
    }
    const auto first = std::begin(haplotypes);
    const auto last  = std::end(haplotypes);
    const auto nth = std::next(first, n);
    const FilterGreater<T> cmp {filter_likelihoods};
    std::nth_element(first, nth, last, cmp);
    result.insert(std::end(result), std::make_move_iterator(nth), std::make_move_iterator(last));
    haplotypes.erase(nth, last);
}

// filters

struct MaxLikelihood
{
    auto operator()(const Haplotype& haplotype, const SampleName& sample,
                    const HaplotypeLikelihoodArray& haplotype_likelihoods) const
    {
        auto result = std::numeric_limits<double>::lowest();
        for (const auto likelihood : haplotype_likelihoods(sample, haplotype)) {
            if (likelihood > result) result = likelihood;
            if (maths::almost_zero(likelihood)) break;
        }
        return result;
    }
    auto operator()(const Haplotype& haplotype, const std::vector<SampleName>& samples,
                    const HaplotypeLikelihoodArray& haplotype_likelihoods) const
    {
        auto result = std::numeric_limits<double>::lowest();
        for (const auto& sample : samples) {
            for (const auto likelihood : haplotype_likelihoods(sample, haplotype)) {
                if (likelihood > result) result = likelihood;
                if (maths::almost_zero(likelihood)) break;
            }
        }
        return result;
    }
};

template <typename Container>
auto count_zeros(const Container& floats) noexcept
{
    return std::count_if(std::cbegin(floats), std::cend(floats),
                         [] (auto value) noexcept { return maths::almost_zero(value); });
}

struct LikelihoodZeroCount
{
    auto operator()(const Haplotype& haplotype, const SampleName& sample,
                    const HaplotypeLikelihoodArray& haplotype_likelihoods) const
    {
        const auto& sample_likelihoods = haplotype_likelihoods(sample, haplotype);
        return count_zeros(sample_likelihoods);
    }
    auto operator()(const Haplotype& haplotype, const std::vector<SampleName>& samples,
                    const HaplotypeLikelihoodArray& haplotype_likelihoods) const
    {
        std::size_t result {0};
        for (const auto& sample : samples) {
            const auto& sample_likelihoods = haplotype_likelihoods(sample, haplotype);
            result += count_zeros(sample_likelihoods);
        }
        return result;
    }
};

class AssignmentCount
{
    std::unordered_map<Haplotype, float> assignments_;
    
public:
    AssignmentCount() = delete;
    
    template <typename ForwardIt>
    AssignmentCount(const ForwardIt first, const ForwardIt last,
                    const std::vector<SampleName>& samples,
                    const HaplotypeLikelihoodArray& haplotype_likelihoods)
    {
        assignments_.reserve(std::distance(first, last));
        std::for_each(first, last, [this] (const auto& haplotype) {
            assignments_.emplace(haplotype, 0.0);
        });
        std::vector<std::reference_wrapper<const Haplotype>> top {};
        for (const auto& sample : samples) {
            const auto n = haplotype_likelihoods.num_likelihoods(sample);
            for (std::size_t i {0}; i < n; ++i) {
                using P = HaplotypeLikelihoodArray::LikelihoodVector::value_type;
                auto cur_max = std::numeric_limits<P>::lowest();
                std::for_each(first, last, [&] (const auto& haplotype) {
                    const auto p = haplotype_likelihoods(sample, haplotype)[i];
                    if (maths::almost_equal(p, cur_max)) {
                        top.emplace_back(haplotype);
                    } else if (p > cur_max) {
                        top.assign({haplotype});
                        cur_max = p;
                    }
                });
                const auto top_score = 1.0f / top.size();
                for (const auto& haplotype : top) {
                    assignments_.at(haplotype) += top_score;
                }
                top.clear();
            }
        }
    }
    
    AssignmentCount(const std::vector<Haplotype>& haplotypes,
                    const std::vector<SampleName>& samples,
                    const HaplotypeLikelihoodArray& haplotype_likelihoods)
    : AssignmentCount {std::cbegin(haplotypes), std::cend(haplotypes), samples, haplotype_likelihoods}
    {}
    
    auto operator()(const Haplotype& haplotype) const
    {
        return assignments_.at(haplotype);;
    }
    
    auto operator()(const Haplotype& haplotype, const std::vector<SampleName>& samples,
                    const HaplotypeLikelihoodArray& haplotype_likelihoods) const
    {
        return assignments_.at(haplotype);
    }
};

struct LikelihoodSum
{
    auto operator()(const Haplotype& haplotype, const std::vector<SampleName>& samples,
                    const HaplotypeLikelihoodArray& haplotype_likelihoods) const
    {
        double result {0};
        for (const auto& sample : samples) {
            const auto& sample_likelihoods = haplotype_likelihoods(sample, haplotype);
            result += std::accumulate(std::cbegin(sample_likelihoods), std::cend(sample_likelihoods), 0.0);
        }
        return result;
    }
};

} // namespace

// main method

std::vector<Haplotype>
filter_to_n(std::vector<Haplotype>& haplotypes, const std::vector<SampleName>& samples,
            const HaplotypeLikelihoodArray& haplotype_likelihoods, const std::size_t n)
{
    std::vector<Haplotype> result {};
    if (haplotypes.size() <= n) {
        return result;
    }
    auto num_to_filter = haplotypes.size() - n;
    result.reserve(haplotypes.size() - n);
    logging::DebugLogger debug_log {};
    logging::TraceLogger trace_log {};
    if (DEBUG_MODE) {
        stream(debug_log) << "Filtering " << num_to_filter << " of "
                          << haplotypes.size() << " haplotypes";
    }
    num_to_filter -= try_filter(haplotypes, samples, haplotype_likelihoods, n, result,
                                MaxLikelihood {}, SamplePoolTag {});
    if (DEBUG_MODE) {
        stream(debug_log) << "There are " << haplotypes.size()
                          << " remaining haplotypes after maximum likelihood filtering";
    }
    if (num_to_filter == 0) {
        return result;
    }
    num_to_filter -= try_filter(haplotypes, samples, haplotype_likelihoods, n, result,
                                AssignmentCount {haplotypes, samples, haplotype_likelihoods},
                                SamplePoolTag {});
    if (DEBUG_MODE) {
        stream(debug_log) << "There are " << haplotypes.size()
                          << " remaining haplotypes after assignment count filtering";
    }
    if (num_to_filter == 0) {
        return result;
    }
    num_to_filter -= try_filter(haplotypes, samples, haplotype_likelihoods, n, result,
                                LikelihoodZeroCount {}, SampleIntersectTag {});
    if (DEBUG_MODE) {
        stream(debug_log) << "There are " << haplotypes.size()
                          << " remaining haplotypes after likelihood zero count filtering";
    }
    if (num_to_filter == 0) {
        return result;
    }
    if (DEBUG_MODE) {
        stream(debug_log) << "Force filtering " << num_to_filter << " of "
                          << haplotypes.size() << " haplotypes with assignment count filtering";
    }
    force_filter(haplotypes, samples, haplotype_likelihoods, n, result,
                 AssignmentCount {haplotypes, samples, haplotype_likelihoods});
    return result;
}

std::vector<HaplotypeReference>
extract_removable(const std::vector<Haplotype>& haplotypes,
                  const HaplotypeReferenceProbabilityMap& haplotype_posteriors,
                  const std::vector<SampleName>& samples,
                  const HaplotypeLikelihoodArray& haplotype_likelihoods,
                  const std::size_t max_to_remove,
                  const double min_posterior)
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::back_inserter;
    std::vector<std::pair<HaplotypeReference, double>> sorted {};
    sorted.reserve(haplotypes.size());
    if (haplotypes.size() == haplotype_posteriors.size()) {
        std::copy(cbegin(haplotype_posteriors), cend(haplotype_posteriors), back_inserter(sorted));
    } else {
        for (const auto& haplotype : haplotypes) {
            sorted.push_back(*haplotype_posteriors.find(haplotype));
        }
    }
    auto first_safe = std::partition(begin(sorted), end(sorted), [=] (const auto& p) { return p.second < min_posterior; });
    const auto num_unsafe = static_cast<std::size_t>(std::distance(begin(sorted), first_safe));
    const auto get_second = [] (const auto& p) { return p.first; };
    if (num_unsafe > max_to_remove) {
        const auto num_unsafe_to_keep = num_unsafe - max_to_remove;
        auto first_unsafe = std::next(begin(sorted), num_unsafe_to_keep);
        std::partial_sort(begin(sorted), first_unsafe, first_safe,
                          [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
        const auto min_safe_posterior = std::prev(first_unsafe)->second;
        const auto equals_min_safe_posterior = [=] (const auto& p) { return maths::almost_equal(p.second, min_safe_posterior); };
        first_unsafe = std::partition(first_unsafe, first_safe, equals_min_safe_posterior);
        if (first_unsafe != first_safe) {
            std::reverse(begin(sorted), first_unsafe);
            const auto assignment_bound_itr = std::prev(first_unsafe);
            const AssignmentCount assignments {boost::make_transform_iterator(assignment_bound_itr, get_second),
                                               boost::make_transform_iterator(first_safe, get_second),
                                               samples, haplotype_likelihoods};
            const auto max_unsafe_assignment = assignments(assignment_bound_itr->first);
            const auto not_greater_than_max_unsafe_assignment = [&assignments, max_unsafe_assignment] (const auto& p) {
                return assignments(p.first) <= max_unsafe_assignment;
            };
            first_safe = std::partition(first_unsafe, first_safe, not_greater_than_max_unsafe_assignment);
            first_safe = std::rotate(begin(sorted), first_unsafe, first_safe);
        } else {
            first_unsafe = std::find_if(begin(sorted), first_unsafe, equals_min_safe_posterior);
            const AssignmentCount assignments {boost::make_transform_iterator(first_unsafe, get_second),
                                               boost::make_transform_iterator(first_safe, get_second),
                                               samples, haplotype_likelihoods};
            const auto assignment_greater = [&] (const auto& lhs, const auto& rhs) {
                const auto lhs_assignments = assignments(lhs.first);
                const auto rhs_assignments = assignments(rhs.first);
                return lhs_assignments != rhs_assignments ? lhs_assignments > rhs_assignments : lhs.first.get() > rhs.first.get();
            };
            std::sort(first_unsafe, first_safe, assignment_greater);
            first_unsafe = std::next(begin(sorted), num_unsafe_to_keep);
            first_safe = std::rotate(begin(sorted), first_unsafe, first_safe);
        }
    }
    std::vector<HaplotypeReference> result {};
    result.reserve(std::distance(begin(sorted), first_safe));
    std::transform(begin(sorted), first_safe, back_inserter(result), get_second);
    return result;
}

} // namespace octopus
