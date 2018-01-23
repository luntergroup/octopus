// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef facet_factory_hpp
#define facet_factory_hpp

#include <string>
#include <vector>
#include <functional>
#include <unordered_map>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/reference/reference_genome.hpp"
#include "readpipe/buffered_read_pipe.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/thread_pool.hpp"
#include "facet.hpp"

namespace octopus { namespace csr {

class FacetFactory
{
public:
    using CallBlock  = std::vector<VcfRecord>;
    using FacetBlock = std::vector<FacetWrapper>;
    
    FacetFactory() = delete;
    
    FacetFactory(const ReferenceGenome& reference, BufferedReadPipe read_pipe);
    
    FacetFactory(const FacetFactory&)            = delete;
    FacetFactory& operator=(const FacetFactory&) = delete;
    FacetFactory(FacetFactory&&);
    FacetFactory& operator=(FacetFactory&&);
    
    ~FacetFactory() = default;
    
    FacetWrapper make(const std::string& name, const CallBlock& block) const;
    FacetBlock make(const std::vector<std::string>& names, const CallBlock& block) const;
    std::vector<FacetBlock> make(const std::vector<std::string>& names, const std::vector<CallBlock>& blocks,
                                 ThreadPool& workers) const;
    
private:
    struct BlockData
    {
        boost::optional<GenomicRegion> region;
        boost::optional<ReadMap> reads;
        boost::optional<GenotypeMap> genotypes;
    };
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    BufferedReadPipe read_pipe_;
    
    std::unordered_map<std::string, std::function<FacetWrapper(const BlockData& data)>> facet_makers_;
    
    void setup_facet_makers();
    FacetWrapper make(const std::string& name, const BlockData& block) const;
    FacetBlock make(const std::vector<std::string>& names, const BlockData& block) const;
    BlockData make_block_data(const std::vector<std::string>& names, const CallBlock& block) const;
};

} // namespace csr
} // namespace octopus

#endif
