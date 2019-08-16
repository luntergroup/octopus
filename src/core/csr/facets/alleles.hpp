// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef alleles_hpp
#define alleles_hpp

#include <vector>
#include <functional>

#include "io/variant/vcf_record.hpp"
#include "facet.hpp"

namespace octopus { namespace csr {

class Alleles : public Facet
{
public:
    using ResultType = std::reference_wrapper<const AlleleMap>;
    
    Alleles() = default;
    Alleles(const std::vector<VcfRecord::SampleName>& samples, const std::vector<VcfRecord>& calls);

private:
    static const std::string name_;
    
    AlleleMap alleles_;
    
    const std::string& do_name() const noexcept override { return name_; }
    Facet::ResultType do_get() const override;
};

std::vector<Allele> copy_overlapped(const MappableFlatSet<Allele>& alleles, const VcfRecord& call);

std::vector<Allele> copy_unique_overlapped(const Facet::AlleleMap& alleles, const VcfRecord& call, const std::vector<SampleName>& samples);

} // namespace csr
} // namespace octopus

#endif
