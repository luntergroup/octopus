// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef somatic_call_hpp
#define somatic_call_hpp

#include <utility>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>
#include <boost/container/flat_map.hpp>

#include "config/common.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/cancer_genotype.hpp"
#include "variant_call.hpp"

namespace octopus {

class SomaticCall : public VariantCall
{
public:
    using VariantCall::GenotypeCall;
    using VariantCall::PhaseCall;
    
    struct AlleleStats
    {
        using CredibleRegion = std::pair<double, double>;
        CredibleRegion vaf_credible_region;
        double map_vaf, pseudo_count;
    };
    
    struct GenotypeAlleleStats
    {
        std::vector<AlleleStats> germline, somatic;
    };
    
    using GenotypeStatsMap = std::unordered_map<SampleName, GenotypeAlleleStats>;
    
    SomaticCall() = delete;
    
    SomaticCall(Variant variant,
                const CancerGenotype<Allele>& genotype_call,
                Phred<double> quality, Phred<double> genotype_posterior,
                GenotypeStatsMap stats,
                boost::optional<Phred<double>> classification_posterior = boost::none);
    
    SomaticCall(const SomaticCall&)            = default;
    SomaticCall& operator=(const SomaticCall&) = default;
    SomaticCall(SomaticCall&&)                 = default;
    SomaticCall& operator=(SomaticCall&&)      = default;
    
    virtual ~SomaticCall() = default;
    virtual void decorate(VcfRecord::Builder& record) const override;
    virtual bool requires_model_evaluation() const noexcept override { return true; }
    
protected:
    struct ExtendedAlleleStats
    {
        AlleleStats stats;
        bool is_somatic_haplotype;
    };

    using SquashedGenotypeStatsMap = boost::container::flat_map<SampleName, std::vector<ExtendedAlleleStats>>;
    
    SquashedGenotypeStatsMap genotype_stats_;
    
private:
    virtual std::unique_ptr<Call> do_clone() const override;
    virtual void reorder_genotype_fields(const SampleName& sample, const std::vector<unsigned>& order) override;
};

} // namespace octopus

#endif
