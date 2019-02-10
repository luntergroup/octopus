// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_call_hpp
#define variant_call_hpp

#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "call.hpp"

namespace octopus {
    
class ReferenceGenome;

class VariantCall : public Call
{
public:
    using Call::GenotypeCall;
    using Call::PhaseCall;
    
    VariantCall() = delete;
    
    template <typename V, typename T>
    VariantCall(V&& variant, T&& genotype_calls, Phred<double> quality, boost::optional<Phred<double>> posterior = boost::none);
    
    virtual ~VariantCall() = default;
    
    VariantCall(const VariantCall&)            = default;
    VariantCall& operator=(const VariantCall&) = default;
    VariantCall(VariantCall&&)                 = default;
    VariantCall& operator=(VariantCall&&)      = default;
    
    const GenomicRegion& mapped_region() const noexcept override;
    
    const Allele& reference() const noexcept override;
    
    const Allele& alternative() const noexcept;
    
    bool is_represented(const Allele& allele) const noexcept override;
    
    void replace(const Allele& old, Allele replacement) override;
    void replace_uncalled_genotype_alleles(const Allele& replacement, char ignore) override;
    
    virtual bool parsimonise(char dummy_base) override;
    virtual bool parsimonise(const ReferenceGenome& reference) override;
    
protected:
    Variant variant_;
    boost::optional<Phred<double>> posterior_;
    
    bool all_genotypes_are_self_contained() const;
    
private:
    void replace_called_alleles(char old_base, char replacement_base) override;
};

template <typename V, typename T>
VariantCall::VariantCall(V&& variant, T&& genotype_calls, Phred<double> quality, boost::optional<Phred<double>> posterior)
: Call {std::forward<T>(genotype_calls), quality}
, variant_ {std::forward<V>(variant)}
, posterior_ {posterior}
{}

} // namespace octopus

#endif
