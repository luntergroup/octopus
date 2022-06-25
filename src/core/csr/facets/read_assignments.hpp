// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_assignments_hpp
#define read_assignments_hpp

#include <unordered_map>
#include <string>
#include <functional>

#include <boost/optional.hpp>

#include "facet.hpp"
#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/haplotype_likelihood_model.hpp"
#include "core/tools/read_assigner.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/thread_pool.hpp"

namespace octopus { namespace csr {

class ReadAssignments : public Facet
{
public:
    using ResultType = std::reference_wrapper<const SupportMaps>;

    using OptionalThreadPool = boost::optional<ThreadPool&>;
    
    ReadAssignments() = default;
    
    ReadAssignments(const ReferenceGenome& reference,
                    const GenotypeMap& genotypes,
                    const ReadMap& reads,
                    const std::vector<VcfRecord>& calls,
                    OptionalThreadPool workers = boost::none);
    ReadAssignments(const ReferenceGenome& reference,
                    const GenotypeMap& genotypes,
                    const ReadMap& reads,
                    const std::vector<VcfRecord>& calls,
                    HaplotypeLikelihoodModel model,
                    OptionalThreadPool workers = boost::none);
    
private:
    static const std::string name_;
    
    SupportMaps result_;
    HaplotypeLikelihoodModel likelihood_model_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

} // namespace csr
} // namespace octopus

#endif
