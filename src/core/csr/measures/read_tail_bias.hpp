// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_tail_bias_hpp
#define read_tail_bias_hpp

#include <string>
#include <vector>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class ReadTailBias : public Measure
{
    const static std::string name_;
    double tail_fraction_ = 0.03;
    std::unique_ptr<Measure> do_clone() const override;
    void do_set_parameters(std::vector<std::string> params) override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    ResultCardinality do_cardinality() const noexcept override;
    const std::string& do_name() const override;
    std::string do_describe() const override;
    std::vector<std::string> do_requirements() const override;
    bool is_equal(const Measure& other) const noexcept override;
};

} // namespace csr
} // namespace octopus

#endif
