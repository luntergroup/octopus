// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef denovo_call_hpp
#define denovo_call_hpp

#include <utility>

#include "variant_call.hpp"

namespace octopus {

class DenovoCall : public VariantCall
{
public:
    using VariantCall::GenotypeCall;
    using VariantCall::PhaseCall;
    
    DenovoCall() = delete;
    
    template <typename V, typename T>
    DenovoCall(V&& variant, T&& genotype_calls, Phred<double> quality, Phred<double> posterior);
    
    DenovoCall(const DenovoCall&)            = default;
    DenovoCall& operator=(const DenovoCall&) = default;
    DenovoCall(DenovoCall&&)                 = default;
    DenovoCall& operator=(DenovoCall&&)      = default;
    
    virtual ~DenovoCall() = default;
    
    virtual void decorate(VcfRecord::Builder& record) const override;
    
    virtual bool requires_model_evaluation() const noexcept override { return true; }
    
private:
    virtual std::unique_ptr<Call> do_clone() const override;
};

template <typename V, typename T>
DenovoCall::DenovoCall(V&& variant, T&& genotype_calls, Phred<double> quality, Phred<double> posterior)
: VariantCall {std::forward<V>(variant), std::forward<T>(genotype_calls), quality, posterior}
{}

} // namespace octopus

#endif
