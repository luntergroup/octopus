// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef unsupervised_clustering_filter_hpp
#define unsupervised_clustering_filter_hpp

#include <vector>
#include <deque>

#include <boost/optional.hpp>

#include "double_pass_variant_call_filter.hpp"

namespace octopus { namespace csr {

class UnsupervisedClusteringFilter : public DoublePassVariantCallFilter
{
public:
    struct Distance
    {
        double operator()(const Measure::ResultType& lhs, const Measure::ResultType& rhs) const;
    };
    
    
    
    UnsupervisedClusteringFilter() = delete;
    
    UnsupervisedClusteringFilter(FacetFactory facet_factory,
                                 std::vector<MeasureWrapper> measures,
                                 OutputOptions output_config,
                                 boost::optional<ProgressMeter&> progress = boost::none);
    
    UnsupervisedClusteringFilter(const UnsupervisedClusteringFilter&)            = delete;
    UnsupervisedClusteringFilter& operator=(const UnsupervisedClusteringFilter&) = delete;
    UnsupervisedClusteringFilter(UnsupervisedClusteringFilter&&)                 = default;
    UnsupervisedClusteringFilter& operator=(UnsupervisedClusteringFilter&&)      = default;
    
    virtual ~UnsupervisedClusteringFilter() override = default;
    
private:
    mutable std::deque<MeasureVector> data_;
    
    virtual void annotate(VcfHeader::Builder& header) const override;
    void record(std::size_t call_idx, MeasureVector measures) const override;
    void prepare_for_classification() const override;
    Classification classify(std::size_t call_idx) const override;
};

} // namespace csr
} // namespace octopus

#endif
