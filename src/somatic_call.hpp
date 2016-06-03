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

#include "common.hpp"
#include "variant_call.hpp"
#include "allele.hpp"
#include "cancer_genotype.hpp"

namespace Octopus
{
    class SomaticCall : public VariantCall
    {
    public:
        using VariantCall::GenotypeCall;
        using VariantCall::PhaseCall;
        
        using CredibleRegion = std::pair<double, double>;
        
        struct GenotypeCredibleRegions
        {
            std::vector<CredibleRegion> germline_credible_regions;
            CredibleRegion somatic_credible_region;
        };
        
        SomaticCall() = delete;
        
        template <typename V, typename C>
        SomaticCall(V&& variant, const CancerGenotype<Allele>& genotype_call,
                    double genotype_posteriors, C&& credible_regions, double quality);
        
        virtual ~SomaticCall() = default;
        
        SomaticCall(const SomaticCall&)            = default;
        SomaticCall& operator=(const SomaticCall&) = default;
        SomaticCall(SomaticCall&&)                 = default;
        SomaticCall& operator=(SomaticCall&&)      = default;
        
        virtual void decorate(VcfRecord::Builder& record) const override;
        
    protected:
        std::unordered_map<SampleIdType, GenotypeCredibleRegions> credible_regions_;
    };
    
    template <typename V, typename C>
    SomaticCall::SomaticCall(V&& variant, const CancerGenotype<Allele>& genotype_call,
                             double genotype_posteriors, C&& credible_regions, double quality)
    :
    VariantCall {std::forward<V>(variant), decltype(genotype_calls_) {}, quality},
    credible_regions_ {std::forward<C>(credible_regions)}
    {
        if (variant_.ref_allele() == variant_.alt_allele()) {
            Allele::SequenceType missing_sequence(ref_sequence_size(variant_), 'N');
            variant_ = Variant {
                Allele {::mapped_region(variant_), std::move(missing_sequence)},
                variant_.alt_allele()
            };
        }
        
        genotype_calls_.reserve(credible_regions_.size()); // num samples
        
        for (const auto& p : credible_regions_) {
            genotype_calls_.emplace(p.first, GenotypeCall {convert(genotype_call), genotype_posteriors});
        }
    }
} // namespace Octopus

#endif /* somatic_call_hpp */
