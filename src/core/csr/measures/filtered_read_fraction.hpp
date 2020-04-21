// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef filtered_read_fraction_hpp
#define filtered_read_fraction_hpp

#include <string>
#include <vector>

#include "measure.hpp"
#include "depth.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class FilteredReadFraction : public Measure
{
    const static std::string name_;
    Depth calling_depth_, filtering_depth_;
    std::unique_ptr<Measure> do_clone() const override;
    ResultType get_default_result() const override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    ResultCardinality do_cardinality() const noexcept override;
    const std::string& do_name() const override;
    std::string do_describe() const override;
    std::vector<std::string> do_requirements() const override;
    bool is_equal(const Measure& other) const noexcept override;
public:
    FilteredReadFraction(bool aggregate_samples = false);
};

} // namespace csr
} // namespace octopus

#endif
