// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef reference_context_hpp
#define reference_context_hpp

#include <functional>
#include <string>

#include <boost/optional.hpp>

#include "basics/genomic_region.hpp"
#include "io/reference/reference_genome.hpp"
#include "facet.hpp"

namespace octopus { namespace csr {

class ReferenceContext : public Facet
{
public:
    using ResultType = ReferenceGenome::GeneticSequence;
    
    ReferenceContext(const ReferenceGenome& reference, GenomicRegion region);
    
private:
    static const std::string name_;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    GenomicRegion region_;
    mutable boost::optional<ResultType> result_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

} // namespace csr
} // namespace octopus

#endif
