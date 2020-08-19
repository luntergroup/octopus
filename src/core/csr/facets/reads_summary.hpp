// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef reads_summary_hpp
#define reads_summary_hpp

#include <string>
#include <functional>

#include "facet.hpp"
#include "config/common.hpp"
#include "basics/aligned_read.hpp"

namespace octopus { namespace csr {

class ReadsSummary : public Facet
{
public:
    using ResultType = std::reference_wrapper<const ReadsSummaryMap>;
    
    ReadsSummary() = default;
    
    ReadsSummary(const ReadMap& reads);

private:
    static const std::string name_;
    
    ReadsSummaryMap result_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

} // namespace csr
} // namespace octopus

#endif
