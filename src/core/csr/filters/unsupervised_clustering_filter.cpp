// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "unsupervised_clustering_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <numeric>

namespace octopus { namespace csr {

UnsupervisedClusteringFilter::UnsupervisedClusteringFilter(FacetFactory facet_factory,
                                                           std::vector<MeasureWrapper> measures,
                                                           OutputOptions output_config,
                                                           ConcurrencyPolicy threading,
                                                           boost::optional<ProgressMeter&> progress)
: DoublePassVariantCallFilter {std::move(facet_factory), std::move(measures), std::move(output_config), threading, progress}
{}

void UnsupervisedClusteringFilter::annotate(VcfHeader::Builder& header) const
{
    // TODO
}

void UnsupervisedClusteringFilter::record(const std::size_t call_idx, std::size_t sample_idx, MeasureVector measures) const
{
    if (data_.size() == call_idx) {
        data_.push_back(std::move(measures));
    } else if (call_idx > data_.size()) {
        data_.resize(call_idx + 1);
        data_.back() = std::move(measures);
    } else {
        data_[call_idx] = std::move(measures);
    }
}

void UnsupervisedClusteringFilter::prepare_for_classification(boost::optional<Log>& log) const
{
    if (log) {
        stream(*log) << "CSR: clustering " << data_.size() << " records";
    }
    remove_missing_features();
    const auto num_calls = data_.size();
    // TODO
    data_.clear();
    data_.shrink_to_fit();
    classifications_.resize(num_calls);
}

VariantCallFilter::Classification UnsupervisedClusteringFilter::classify(std::size_t call_idx, std::size_t sample_idx) const
{
    return classifications_[call_idx];
}

bool UnsupervisedClusteringFilter::all_missing(const MeasureVector& measures) const noexcept
{
    return std::all_of(std::cbegin(measures), std::cend(measures),
                       [] (const auto& measure) noexcept { return is_missing(measure); });
}

void UnsupervisedClusteringFilter::remove_missing_features() const
{
    if (!data_.empty()) {
        const auto num_dimensions = data_.front().size();
        std::vector<std::size_t> missing_columns(num_dimensions);
        std::iota(std::begin(missing_columns), std::end(missing_columns), 0);
        for (const auto& point : data_) {
            for (auto itr = std::cbegin(missing_columns); itr != std::cend(missing_columns);) {
                assert(*itr < point.size());
                if (is_missing(point[*itr])) {
                    ++itr;
                } else {
                    itr = missing_columns.erase(itr);
                }
            }
            if (missing_columns.empty()) break;
        }
        if (!missing_columns.empty()) {
            if (missing_columns.size() != num_dimensions) {
                assert(std::is_sorted(std::cbegin(missing_columns), std::cend(missing_columns)));
                for (auto& point : data_) {
                    std::for_each(std::crbegin(missing_columns), std::crend(missing_columns),
                                  [&point] (auto idx) { point.erase(std::next(std::cbegin(point), idx)); });
        
                }
            } else {
                data_.clear();
                data_.shrink_to_fit();
            }
        }
    }
}

} // namespace csr
} // namespace octopus
