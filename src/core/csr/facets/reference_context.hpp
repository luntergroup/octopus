// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef reference_context_hpp
#define reference_context_hpp

#include <string>
#include <functional>
#include <memory>

#include <boost/optional.hpp>

#include "basics/genomic_region.hpp"
#include "core/types/haplotype.hpp"
#include "io/reference/reference_genome.hpp"
#include "facet.hpp"

namespace octopus { namespace csr {

class ReferenceContext : public Facet
{
public:
    using ResultType = std::reference_wrapper<const Haplotype>;
    
    ReferenceContext() = default;
    
    ReferenceContext(const ReferenceGenome& reference, GenomicRegion region);
    
private:
    static const std::string name_;
    
    std::unique_ptr<Haplotype> result_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

} // namespace csr
} // namespace octopus

#endif
