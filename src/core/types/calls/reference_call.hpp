// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef reference_call_hpp
#define reference_call_hpp

#include <map>
#include <utility>

#include "call.hpp"

namespace octopus {

class ReferenceCall : public Call
{
public:
    ReferenceCall() = delete;
    
    struct GenotypeCall
    {
        unsigned ploidy;
        Phred<double> posterior;
    };
    
    template <typename A>
    ReferenceCall(A&& reference, Phred<double> quality, std::map<SampleName, GenotypeCall> genotypes);
    
    ReferenceCall(const ReferenceCall&)            = default;
    ReferenceCall& operator=(const ReferenceCall&) = default;
    ReferenceCall(ReferenceCall&&)                 = default;
    ReferenceCall& operator=(ReferenceCall&&)      = default;
    
    virtual ~ReferenceCall() = default;
    
    const GenomicRegion& mapped_region() const noexcept override;
    
    const Allele& reference() const noexcept override;
    
    bool is_represented(const Allele& allele) const noexcept override;
    
    void replace(const Allele& old, Allele replacement) override;
    void replace_uncalled_genotype_alleles(const Allele& replacement, char ignore) override;
    
    void decorate(VcfRecord::Builder& record) const override;
    
private:
    Allele reference_;
    
    virtual std::unique_ptr<Call> do_clone() const override;
    void replace_called_alleles(char old_base, char replacement_base) override;
};

template <typename A>
ReferenceCall::ReferenceCall(A&& reference, Phred<double> quality, std::map<SampleName, GenotypeCall> genotypes)
: Call {quality}
, reference_ {std::forward<A>(reference)}
{
    genotype_calls_.reserve(genotypes.size());
    for (const auto& p : genotypes) {
        genotype_calls_.emplace(std::piecewise_construct,
                                std::forward_as_tuple(p.first),
                                std::forward_as_tuple(Genotype<Allele> {p.second.ploidy, reference_}, p.second.posterior));
    }
}

} // namespace octopus

#endif
