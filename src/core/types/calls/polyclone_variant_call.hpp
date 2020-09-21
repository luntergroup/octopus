// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef polyclone_variant_call_hpp
#define polyclone_variant_call_hpp

#include "variant_call.hpp"

#include <vector>
#include <utility>

#include <boost/optional.hpp>

#include "core/types/allele.hpp"

namespace octopus {

class PolycloneVariantCall : public VariantCall
{
public:
    using VariantCall::GenotypeCall;
    using VariantCall::PhaseCall;

    struct HaplotypeFrequencyStats
    {
        double pseudo_count, map;
    };
    using HaplotypeFrequencyStatsVector = std::vector<HaplotypeFrequencyStats>;
    
    PolycloneVariantCall() = delete;
    
    template <typename V, typename T>
    PolycloneVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality);
    template <typename V, typename T>
    PolycloneVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality,
                        HaplotypeFrequencyStatsVector haplotype_frequency_stats);
    
    PolycloneVariantCall(const PolycloneVariantCall&)            = default;
    PolycloneVariantCall& operator=(const PolycloneVariantCall&) = default;
    PolycloneVariantCall(PolycloneVariantCall&&)                 = default;
    PolycloneVariantCall& operator=(PolycloneVariantCall&&)      = default;
    
    virtual ~PolycloneVariantCall() = default;
    
    void decorate(VcfRecord::Builder& record) const override;
    
private:
    boost::optional<HaplotypeFrequencyStatsVector> haplotype_frequency_stats_;

    virtual std::unique_ptr<Call> do_clone() const override;
};

template <typename V, typename T>
PolycloneVariantCall::PolycloneVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality)
: VariantCall {std::forward<V>(variant), std::forward<T>(genotype_calls), quality}
, haplotype_frequency_stats_ {}
{}

template <typename V, typename T>
PolycloneVariantCall::PolycloneVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality,
                                           HaplotypeFrequencyStatsVector haplotype_frequency_stats)
: VariantCall {std::forward<V>(variant), std::forward<T>(genotype_calls), quality}
, haplotype_frequency_stats_ {std::move(haplotype_frequency_stats)}
{}

} // namespace octopus

#endif
