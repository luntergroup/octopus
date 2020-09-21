// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef repeat_context_hpp
#define repeat_context_hpp

#include <string>
#include <functional>
#include <memory>

#include <boost/optional.hpp>

#include "basics/genomic_region.hpp"
#include "io/reference/reference_genome.hpp"
#include "basics/tandem_repeat.hpp"
#include "facet.hpp"

namespace octopus { namespace csr {

class RepeatContext : public Facet
{
public:
    using ResultType = std::reference_wrapper<const std::vector<TandemRepeat>>;
    
    RepeatContext() = default;
    
    RepeatContext(const ReferenceGenome& reference, GenomicRegion region);
    
private:
    static const std::string name_;
    
    std::vector<TandemRepeat> result_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

} // namespace csr
} // namespace octopus

#endif
