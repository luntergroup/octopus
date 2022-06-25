// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef model_posterior_by_depth_hpp
#define model_posterior_by_depth_hpp

#include <string>
#include <vector>

#include "measure.hpp"
#include "depth.hpp"

namespace octopus {

class VcfRecord;

namespace csr  {

class ModelPosteriorByDepth : public Measure
{
    const static std::string name_;
    Depth depth_;
    std::unique_ptr<Measure> do_clone() const override;
    ValueType get_value_type() const override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    ResultCardinality do_cardinality() const noexcept override;
    const std::string& do_name() const override;
    std::string do_describe() const override;
    std::vector<std::string> do_requirements() const override;
    bool is_equal(const Measure& other) const noexcept override;
public:
    ModelPosteriorByDepth(bool recalculate = false);
};

} // namespace csr
} // namespace octopus

#endif
