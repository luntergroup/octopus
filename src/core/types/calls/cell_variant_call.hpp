// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cell_variant_call_hpp
#define cell_variant_call_hpp

#include "variant_call.hpp"

#include <vector>
#include <utility>
#include <unordered_map>

#include <boost/optional.hpp>

#include "basics/phred.hpp"
#include "core/types/allele.hpp"
#include "core/types/phylogeny.hpp"

namespace octopus {

class CellVariantCall : public VariantCall
{
public:
    using VariantCall::GenotypeCall;
    using VariantCall::PhaseCall;
    
    struct PhylogenyInferenceSummary
    {
        Phylogeny<std::size_t> map;
        Phred<double> map_posterior;
        std::vector<Phred<double>> size_posteriors;
        std::unordered_map<SampleName, std::vector<Phred<double>>> sample_node_posteriors;
    };
    
    CellVariantCall() = delete;
    
    template <typename V, typename T>
    CellVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality, PhylogenyInferenceSummary phylogeny_summary);
    
    CellVariantCall(const CellVariantCall&)            = default;
    CellVariantCall& operator=(const CellVariantCall&) = default;
    CellVariantCall(CellVariantCall&&)                 = default;
    CellVariantCall& operator=(CellVariantCall&&)      = default;
    
    virtual ~CellVariantCall() = default;
    
    void decorate(VcfRecord::Builder& record) const override;

private:
    PhylogenyInferenceSummary phylogeny_summary_;
    
    virtual std::unique_ptr<Call> do_clone() const override;
    
    bool is_somatic() const;
};

template <typename V, typename T>
CellVariantCall::CellVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality, PhylogenyInferenceSummary phylogeny_summary)
: VariantCall {std::forward<V>(variant), std::forward<T>(genotype_calls), quality}
, phylogeny_summary_ {std::move(phylogeny_summary)}
{}

} // namespace octopus

#endif
