// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef strand_bias_hpp
#define strand_bias_hpp

#include <string>
#include <vector>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class StrandBias : public Measure
{
    std::unique_ptr<Measure> do_clone() const override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    ResultCardinality do_cardinality() const noexcept override;
    std::string do_name() const override;
    std::vector<std::string> do_requirements() const override;
    
    double min_difference_ = 0.25;
    std::size_t small_sample_size_ = 200, medium_sample_size_ = 1'000,
                big_sample_size_ = 10'000, very_big_sample_size = 100'000;
    double min_medium_trigger_, min_big_trigger_;
    double critical_resample_lb_, critical_resample_ub_;
    bool use_resampling_ = false;
    
public:
    StrandBias() = default;
    StrandBias(double critical_value);
};

} // namespace csr
} // namespace octopus

#endif
