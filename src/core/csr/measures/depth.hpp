// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef depth_hpp
#define depth_hpp

#include <string>
#include <vector>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class Depth : public Measure
{
    const static std::string name_;
    bool recalculate_ = false, aggregate_ = false;
    std::unique_ptr<Measure> do_clone() const override;
    ValueType get_value_type() const override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    ResultCardinality do_cardinality() const noexcept override;
    const std::string& do_name() const override;
    std::string do_describe() const override;
    std::vector<std::string> do_requirements() const override;
    bool is_equal(const Measure& other) const noexcept override;
    void do_set_parameters(std::vector<std::string> params) override;
    std::vector<std::string> do_parameters() const override;
public:
    Depth();
    Depth(bool recalculate, bool aggregate_samples);
};

} // namespace csr
} // namespace octopus

#endif
