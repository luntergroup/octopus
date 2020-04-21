// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef ploidies_hpp
#define ploidies_hpp

#include <string>
#include <vector>
#include <functional>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "basics/ploidy_map.hpp"
#include "facet.hpp"

namespace octopus { namespace csr {

class Ploidies : public Facet
{
public:
    using ResultType = std::reference_wrapper<const LocalPloidyMap>;
    
    Ploidies() = default;
    Ploidies(const PloidyMap& ploidies, const GenomicRegion& region, const std::vector<SampleName>& samples);

private:
    static const std::string name_;
    
    LocalPloidyMap ploidies_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

} // namespace csr
} // namespace octopus

#endif
