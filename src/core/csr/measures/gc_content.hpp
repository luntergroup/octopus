// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef gc_content_hpp
#define gc_content_hpp

#include <string>
#include <vector>

#include "measure.hpp"

namespace octopus {

class VcfRecord;

namespace csr {

class GCContent : public Measure
{
    const static std::string name_;
    std::unique_ptr<Measure> do_clone() const override;
    ValueType get_value_type() const override;
    ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const override;
    ResultCardinality do_cardinality() const noexcept override;
    const std::string& do_name() const override;
    std::string do_describe() const override;
    std::vector<std::string> do_requirements() const override;
};

} // namespace csr
} // namespace octopus

#endif
