// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef unsupervised_clustering_filter_factory_hpp
#define unsupervised_clustering_filter_factory_hpp

#include <memory>
#include <vector>
#include <set>
#include <string>

#include <boost/optional.hpp>

#include "logging/progress_meter.hpp"
#include "../measures/measure.hpp"
#include "variant_call_filter_factory.hpp"
#include "variant_call_filter.hpp"
#include "unsupervised_clustering_filter.hpp"

namespace octopus { namespace csr {

class FacetFactory;

class UnsupervisedClusteringFilterFactory : public VariantCallFilterFactory
{
public:
    UnsupervisedClusteringFilterFactory() = default;
    
    UnsupervisedClusteringFilterFactory(const std::set<std::string>& measure_names);
    
    UnsupervisedClusteringFilterFactory(const UnsupervisedClusteringFilterFactory&)            = default;
    UnsupervisedClusteringFilterFactory& operator=(const UnsupervisedClusteringFilterFactory&) = default;
    UnsupervisedClusteringFilterFactory(UnsupervisedClusteringFilterFactory&&)                 = default;
    UnsupervisedClusteringFilterFactory& operator=(UnsupervisedClusteringFilterFactory&&)      = default;
    
    ~UnsupervisedClusteringFilterFactory() = default;

private:
    std::vector<MeasureWrapper> measures_;
    
    std::unique_ptr<VariantCallFilterFactory> do_clone() const override;
    std::unique_ptr<VariantCallFilter> do_make(FacetFactory facet_factory,
                                               VariantCallFilter::OutputOptions output_config,
                                               boost::optional<ProgressMeter&> progress,
                                               VariantCallFilter::ConcurrencyPolicy threading) const override;
};

} // namespace csr

using csr::UnsupervisedClusteringFilterFactory;

} // namespace octopus


#endif
