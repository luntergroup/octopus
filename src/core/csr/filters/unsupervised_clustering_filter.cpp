// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "unsupervised_clustering_filter.hpp"

#include <utility>
#include <iterator>
#include <algorithm>

#include <iostream>
#include <boost/optional/optional_io.hpp>

namespace octopus { namespace csr {

UnsupervisedClusteringFilter::UnsupervisedClusteringFilter(FacetFactory facet_factory,
                                                           std::vector<MeasureWrapper> measures,
                                                           OutputOptions output_config,
                                                           boost::optional<ProgressMeter&> progress)
: DoublePassVariantCallFilter {std::move(facet_factory), std::move(measures), std::move(output_config), progress}
{}

void UnsupervisedClusteringFilter::annotate(VcfHeader::Builder& header) const
{
    // TODO
}

void UnsupervisedClusteringFilter::record(const std::size_t call_idx, MeasureVector measures) const
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

void UnsupervisedClusteringFilter::prepare_for_classification() const
{
    for (const auto& point : data_) {
        std::copy(std::cbegin(point), std::cend(point), std::ostream_iterator<Measure::ResultType> {std::cout, "\t"});
        std::cout << '\n';
    }
    // TODO
    data_.clear();
    data_.shrink_to_fit();
}

VariantCallFilter::Classification UnsupervisedClusteringFilter::classify(std::size_t call_idx) const
{
    // TODO
    return Classification {};
}

} // namespace csr
} // namespace octopus
