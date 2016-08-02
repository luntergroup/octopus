//
//  somatic_call.hpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef somatic_call_hpp
#define somatic_call_hpp

#include <utility>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include <config/common.hpp>
#include "variant_call.hpp"
#include <core/types/allele.hpp>
#include <core/types/cancer_genotype.hpp>

namespace octopus {

class SomaticCall : public VariantCall
{
public:
    using VariantCall::GenotypeCall;
    using VariantCall::PhaseCall;
    
    using CredibleRegion = std::pair<double, double>;
    
    struct GenotypeCredibleRegions
    {
        std::vector<CredibleRegion> germline;
        boost::optional<CredibleRegion> somatic;
    };
    
    SomaticCall() = delete;
    
    template <typename V, typename C>
    SomaticCall(V&& variant, const CancerGenotype<Allele>& genotype_call,
                Phred<double> genotype_posterior, C&& credible_regions, Phred<double> quality);
    
    virtual ~SomaticCall() = default;
    
    SomaticCall(const SomaticCall&)            = default;
    SomaticCall& operator=(const SomaticCall&) = default;
    SomaticCall(SomaticCall&&)                 = default;
    SomaticCall& operator=(SomaticCall&&)      = default;
    
    virtual void decorate(VcfRecord::Builder& record) const override;
    
protected:
    std::unordered_map<SampleName, GenotypeCredibleRegions> credible_regions_;
};

template <typename V, typename C>
SomaticCall::SomaticCall(V&& variant, const CancerGenotype<Allele>& genotype_call,
                         Phred<double> genotype_posterior, C&& credible_regions, Phred<double> quality)
:
VariantCall {std::forward<V>(variant), decltype(genotype_calls_) {}, quality},
credible_regions_ {std::forward<C>(credible_regions)}
{
    if (variant_.ref_allele() == variant_.alt_allele()) {
        Allele::NucleotideSequence missing_sequence(ref_sequence_size(variant_), 'N');
        using octopus::mapped_region;
        variant_ = Variant {
            Allele {mapped_region(variant_), std::move(missing_sequence)},
            variant_.alt_allele()
        };
    }
    
    genotype_calls_.reserve(credible_regions_.size()); // num samples
    
    for (const auto& p : credible_regions_) {
        if (p.second.somatic) {
            genotype_calls_.emplace(p.first, GenotypeCall {convert(genotype_call), genotype_posterior});
        } else {
            genotype_calls_.emplace(p.first, GenotypeCall {genotype_call.germline_genotype(), genotype_posterior});
        }
    }
}

} // namespace octopus

#endif /* somatic_call_hpp */
