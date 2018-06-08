// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "conditional_threshold_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>

#include "utils/append.hpp"

namespace octopus { namespace csr {

namespace {

auto concat(std::vector<ThresholdVariantCallFilter::ConditionVectorPair> conditions)
{
    ThresholdVariantCallFilter::ConditionVectorPair result {};
    for (auto& p : conditions) {
        utils::append(p.hard, result.hard);
        utils::append(p.soft, result.soft);
    }
    return result;
}

bool are_all_unique(std::vector<std::string> keys)
{
    std::sort(std::begin(keys), std::end(keys));
    return std::adjacent_find(std::cbegin(keys), std::cend(keys)) == std::cend(keys);
}

} // namespace

ConditionalThresholdVariantCallFilter::ConditionalThresholdVariantCallFilter(FacetFactory facet_factory,
                                                                             std::vector<ConditionVectorPair> conditions,
                                                                             std::vector<MeasureWrapper> chooser_measures,
                                                                             std::function<std::size_t(std::vector<Measure::ResultType>)> chooser,
                                                                             OutputOptions output_config,
                                                                             ConcurrencyPolicy threading,
                                                                             boost::optional<ProgressMeter&> progress)
: ThresholdVariantCallFilter {std::move(facet_factory), concat(conditions), output_config, threading, progress}
, hard_ranges_ {}
, soft_ranges_ {}
, chooser_ {std::move(chooser)}
, unique_filter_keys_ {}
, num_chooser_measures_ {chooser_measures.size()}
{
    utils::append(std::move(chooser_measures), measures_);
    measures_.shrink_to_fit();
    hard_ranges_.reserve(conditions.size());
    soft_ranges_.reserve(conditions.size());
    unique_filter_keys_.reserve(conditions.size());
    std::size_t i {0}, j {0}, k {0};
    for (const auto& p : conditions) {
        hard_ranges_.push_back({i, i + p.hard.size(), j});
        i += p.hard.size(); j += p.hard.size();
        soft_ranges_.push_back({i, i + p.soft.size(), k});
        i += p.soft.size();
        std::vector<std::string> filter_keys {std::next(std::cbegin(vcf_filter_keys_), k),
                                              std::next(std::cbegin(vcf_filter_keys_), k + p.soft.size())};
        unique_filter_keys_.push_back(are_all_unique(filter_keys));
        k += p.soft.size();
    }
}

bool ConditionalThresholdVariantCallFilter::passes_all_hard_filters(const MeasureVector& measures) const
{
    return passes_all_hard_filters(measures, hard_ranges_[choose_filter(measures)]);
}

bool ConditionalThresholdVariantCallFilter::passes_all_soft_filters(const MeasureVector& measures) const
{
    return passes_all_soft_filters(measures, soft_ranges_[choose_filter(measures)]);
}

std::vector<std::string> ConditionalThresholdVariantCallFilter::get_failing_vcf_filter_keys(const MeasureVector& measures) const
{
    const auto filter_idx = choose_filter(measures);
    auto result = get_failing_vcf_filter_keys(measures, soft_ranges_[filter_idx]);
    if (!unique_filter_keys_[filter_idx]) {
        std::sort(std::begin(result), std::end(result));
        result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    }
    return result;
}

std::size_t ConditionalThresholdVariantCallFilter::choose_filter(const MeasureVector& measures) const
{
    const MeasureVector chooser_measures(std::prev(std::cend(measures), num_chooser_measures_), std::cend(measures));
    return chooser_(chooser_measures);
}

bool ConditionalThresholdVariantCallFilter::passes_all_hard_filters(const MeasureVector& measures, const MeasureIndexRange range) const
{
    using std::cbegin; using std::next;
    return passes_all_filters(next(cbegin(measures), range.measure_begin), next(cbegin(measures), range.measure_end),
                              next(cbegin(hard_thresholds_), range.threshold_begin));
}

bool ConditionalThresholdVariantCallFilter::passes_all_soft_filters(const MeasureVector& measures, const MeasureIndexRange range) const
{
    using std::cbegin; using std::next;
    return passes_all_filters(next(cbegin(measures), range.measure_begin), next(cbegin(measures), range.measure_end),
                              next(cbegin(soft_thresholds_), range.threshold_begin));
}

std::vector<std::string>
ConditionalThresholdVariantCallFilter::get_failing_vcf_filter_keys(const MeasureVector& measures, const MeasureIndexRange range) const
{
    std::vector<std::string> result {};
    const auto num_conditions = range.measure_end - range.measure_begin;
    result.reserve(num_conditions);
    for (std::size_t i {0}; i < num_conditions; ++i) {
        if (!soft_thresholds_[range.threshold_begin + i](measures[range.measure_begin + i])) {
            result.push_back(vcf_filter_keys_[range.threshold_begin + i]);
        }
    }
    return result;
}

} // namespace csr
} // namespace octopus
