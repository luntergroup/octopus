// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_threshold_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <cassert>

#include <boost/variant.hpp>

#include "utils/concat.hpp"
#include "../measures/is_somatic.hpp"

namespace octopus { namespace csr {

namespace {

bool are_all_unique(std::vector<std::string> keys)
{
    std::sort(std::begin(keys), std::end(keys));
    return std::adjacent_find(std::cbegin(keys), std::cend(keys)) == std::cend(keys);
}

} // namespace

SomaticThresholdVariantCallFilter::SomaticThresholdVariantCallFilter(FacetFactory facet_factory,
                                                                     std::vector<Condition> germline_hard_conditions,
                                                                     std::vector<Condition> germline_soft_conditions,
                                                                     std::vector<Condition> somatic_hard_conditions,
                                                                     std::vector<Condition> somatic_soft_conditions,
                                                                     OutputOptions output_config,
                                                                     ConcurrencyPolicy threading,
                                                                     boost::optional<ProgressMeter&> progress)
: ThresholdVariantCallFilter {std::move(facet_factory),
                              concat(germline_hard_conditions, std::move(somatic_hard_conditions)),
                              concat(germline_soft_conditions, std::move(somatic_soft_conditions)),
                              output_config, threading, progress}
, num_germline_hard_conditions_ {germline_hard_conditions.size()}
, num_germline_soft_conditions_ {germline_soft_conditions.size()}
{
    measures_.emplace_back(make_wrapped_measure<IsSomatic>());
    std::vector<std::string> germline_filter_keys {std::cbegin(vcf_filter_keys_),
                                                   std::next(std::cbegin(vcf_filter_keys_), num_germline_soft_conditions_)};
    all_unique_germline_filter_keys_ = are_all_unique(germline_filter_keys);
    
    std::vector<std::string> somatic_filter_keys {std::next(std::cbegin(vcf_filter_keys_), num_germline_soft_conditions_),
                                                  std::cend(vcf_filter_keys_)};
    all_unique_somatic_filter_keys_ = are_all_unique(somatic_filter_keys);
}

bool SomaticThresholdVariantCallFilter::passes_all_hard_filters(const MeasureVector& measures) const
{
    assert(measures.size() > 1);
    const auto is_somatic = boost::get<bool>(measures.back());
    if (is_somatic) {
        return passes_all_somatic_hard_filters(measures);
    } else {
        return passes_all_germline_hard_filters(measures);
    }
}

bool SomaticThresholdVariantCallFilter::passes_all_soft_filters(const MeasureVector& measures) const
{
    assert(measures.size() > 1);
    const auto is_somatic = boost::get<bool>(measures.back());
    if (is_somatic) {
        return passes_all_somatic_soft_filters(measures);
    } else {
        return passes_all_germline_soft_filters(measures);
    }
}

std::vector<std::string> SomaticThresholdVariantCallFilter::get_failing_vcf_filter_keys(const MeasureVector& measures) const
{
    assert(measures.size() > 1);
    const auto is_somatic = boost::get<bool>(measures.back());
    if (is_somatic) {
        return get_failing_somatic_vcf_filter_keys(measures);
    } else {
        return get_failing_germline_vcf_filter_keys(measures);
    }
}

bool SomaticThresholdVariantCallFilter::passes_all_germline_hard_filters(const MeasureVector& measures) const
{
    return passes_all_filteres(std::cbegin(measures), std::next(std::cbegin(measures), num_germline_hard_conditions_),
                               std::cbegin(hard_thresholds_));
}

bool SomaticThresholdVariantCallFilter::passes_all_somatic_hard_filters(const MeasureVector& measures) const
{
    const auto num_somatic_hard_conditions = hard_thresholds_.size() - num_germline_hard_conditions_;
    auto first_measure_itr = std::next(std::cbegin(measures), num_germline_hard_conditions_);
    return passes_all_filteres(first_measure_itr, std::next(first_measure_itr, num_somatic_hard_conditions),
                               std::next(std::cbegin(hard_thresholds_), num_germline_soft_conditions_));
}

bool SomaticThresholdVariantCallFilter::passes_all_germline_soft_filters(const MeasureVector& measures) const
{
    auto first_measure_itr = std::next(std::cbegin(measures), hard_thresholds_.size());
    return passes_all_filteres(first_measure_itr, std::next(first_measure_itr, num_germline_soft_conditions_),
                               std::cbegin(soft_thresholds_));
}

bool SomaticThresholdVariantCallFilter::passes_all_somatic_soft_filters(const MeasureVector& measures) const
{
    const auto num_somatic_soft_conditions = soft_thresholds_.size() - num_germline_soft_conditions_;
    auto first_measure_itr = std::next(std::cbegin(measures), hard_thresholds_.size() + num_germline_soft_conditions_);
    return passes_all_filteres(first_measure_itr, std::next(first_measure_itr, num_somatic_soft_conditions),
                               std::next(std::cbegin(soft_thresholds_), num_germline_soft_conditions_));
}

std::vector<std::string> SomaticThresholdVariantCallFilter::get_failing_germline_vcf_filter_keys(const MeasureVector& measures) const
{
    std::vector<std::string> result {};
    result.reserve(num_germline_soft_conditions_);
    for (std::size_t i {0}; i < num_germline_soft_conditions_; ++i) {
        if (!soft_thresholds_[i](measures[i + hard_thresholds_.size()])) {
            result.push_back(vcf_filter_keys_[i]);
        }
    }
    return result;
}

std::vector<std::string> SomaticThresholdVariantCallFilter::get_failing_somatic_vcf_filter_keys(const MeasureVector& measures) const
{
    std::vector<std::string> result {};
    const auto num_somatic_soft_conditions = soft_thresholds_.size() - num_germline_soft_conditions_;
    result.reserve(num_somatic_soft_conditions);
    for (std::size_t i {0}; i < num_somatic_soft_conditions; ++i) {
        const auto j = i + num_germline_soft_conditions_;
        if (!soft_thresholds_[j](measures[j + hard_thresholds_.size()])) {
            result.push_back(vcf_filter_keys_[i]);
        }
    }
    return result;
}

} // namespace csr
} // namespace octopus
