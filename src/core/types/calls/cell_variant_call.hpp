// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cell_variant_call_hpp
#define cell_variant_call_hpp

#include "variant_call.hpp"

#include <vector>
#include <utility>

#include <boost/optional.hpp>

#include "core/types/allele.hpp"

namespace octopus {

class CellVariantCall : public VariantCall
{
public:
    using VariantCall::GenotypeCall;
    using VariantCall::PhaseCall;
    
    CellVariantCall() = delete;
    
    template <typename V, typename T>
    CellVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality);
    
    CellVariantCall(const CellVariantCall&)            = default;
    CellVariantCall& operator=(const CellVariantCall&) = default;
    CellVariantCall(CellVariantCall&&)                 = default;
    CellVariantCall& operator=(CellVariantCall&&)      = default;
    
    virtual ~CellVariantCall() = default;
    
    void decorate(VcfRecord::Builder& record) const override;

private:
    virtual std::unique_ptr<Call> do_clone() const override;
    
    bool is_somatic() const;
};

template <typename V, typename T>
CellVariantCall::CellVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality)
: VariantCall {std::forward<V>(variant), std::forward<T>(genotype_calls), quality}
{}

} // namespace octopus

#endif
