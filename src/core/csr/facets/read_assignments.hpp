// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_assignments_hpp
#define read_assignments_hpp

#include <unordered_map>
#include <string>

#include <boost/optional.hpp>

#include "facet.hpp"
#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/tools/read_assigner.hpp"
#include "io/reference/reference_genome.hpp"

namespace octopus { namespace csr {

class ReadAssignments : public Facet
{
public:
    using GenotypeMap = std::unordered_map<SampleName, MappableFlatSet<Genotype<Haplotype>>>;
    using ResultType = std::unordered_map<SampleName, HaplotypeSupportMap>;
    
    ReadAssignments(const ReferenceGenome& reference, const GenotypeMap& genotypes, const ReadMap& reads);
    
private:
    static const std::string name_;
    
    ResultType result_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

} // namespace csr
} // namespace octopus

#endif
