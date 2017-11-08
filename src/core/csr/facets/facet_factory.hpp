// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef facet_factory_hpp
#define facet_factory_hpp

#include <string>
#include <vector>
#include <functional>
#include <unordered_map>

#include "facet.hpp"
#include "io/variant/vcf_record.hpp"
#include "readpipe/buffered_read_pipe.hpp"
#include "io/reference/reference_genome.hpp"

namespace octopus { namespace csr {

class FacetFactory
{
public:
    FacetFactory() = delete;
    
    FacetFactory(const ReferenceGenome& reference, BufferedReadPipe read_pipe);
    
    FacetFactory(const FacetFactory&)            = delete;
    FacetFactory& operator=(const FacetFactory&) = delete;
    FacetFactory(FacetFactory&&);
    FacetFactory& operator=(FacetFactory&&);
    
    ~FacetFactory() = default;
    
    FacetWrapper make(const std::string& name, const std::vector<VcfRecord>& records) const;
    
private:
    std::reference_wrapper<const ReferenceGenome> reference_;
    BufferedReadPipe read_pipe_;
    std::unordered_map<std::string, std::function<FacetWrapper(const std::vector<VcfRecord>& records)>> facet_makers_;
    
    void setup_facet_makers();
};

} // namespace csr
} // namespace octopus

#endif
