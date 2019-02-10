// Copyright (c) 2015-2019 Daniel Cooke
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
#include "core/tools/read_assigner.hpp"
#include "io/reference/reference_genome.hpp"

namespace octopus { namespace csr {

class ReadAssignments : public Facet
{
public:
    using ResultType = std::reference_wrapper<const SupportMaps>;
    
    ReadAssignments() = default;
    
    ReadAssignments(const ReferenceGenome& reference, const GenotypeMap& genotypes, const ReadMap& reads);
    
private:
    static const std::string name_;
    
    SupportMaps result_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

AlleleSupportMap
compute_allele_support(const std::vector<Allele>& alleles, const Facet::SupportMaps& assignments, const SampleName& sample);

} // namespace csr
} // namespace octopus

#endif
