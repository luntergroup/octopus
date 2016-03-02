//
//  cancer_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 16/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__cancer_caller__
#define __Octopus__cancer_caller__

#include <string>
#include <memory>

#include "common.hpp"
#include "variant_caller.hpp"
#include "cancer_genotype_model.hpp"

class GenomicRegion;
class ReadPipe;
class Variant;
class VcfRecord;

namespace Octopus
{
    class CancerVariantCaller : public VariantCaller
    {
    public:
        CancerVariantCaller() = delete;
        
        explicit CancerVariantCaller(const ReferenceGenome& reference,
                                     ReadPipe& read_pipe,
                                     CandidateVariantGenerator&& candidate_generator,
                                     RefCallType refcalls,
                                     double min_variant_posterior,
                                     double min_somatic_posterior,
                                     double min_refcall_posterior,
                                     const SampleIdType& normal_sample,
                                     bool call_somatics_only);
        
        ~CancerVariantCaller() = default;
        
        CancerVariantCaller(const CancerVariantCaller&)            = delete;
        CancerVariantCaller& operator=(const CancerVariantCaller&) = delete;
        CancerVariantCaller(CancerVariantCaller&&)                 = delete;
        CancerVariantCaller& operator=(CancerVariantCaller&&)      = delete;
        
    private:
        using VariantCaller::HaplotypePriorMap;
        
        struct Latents : public CallerLatents
        {
            using CallerLatents::HaplotypePosteiorMap;
            using CallerLatents::GenotypePosteriorMap;
            
            template <typename T1, typename T2, typename T3>
            Latents(T1 g, T2 c, T3 m) : germline_genotype_posteriors {g}, cancer_haplotype_posteriors {c}, genotype_mixtures {m} {}
            
            std::shared_ptr<HaplotypePosteiorMap> get_haplotype_posteriors() const override
            {
                return std::make_shared<HaplotypePosteiorMap>();
            }
            
            std::shared_ptr<GenotypePosteriorMap> get_genotype_posteriors() const override
            {
                ProbabilityMatrix<Genotype<Haplotype>> result {};
                assign_keys(extract_keys(germline_genotype_posteriors), result);
                insert_sample("test", extract_values(germline_genotype_posteriors), result);
                return std::make_shared<GenotypePosteriorMap>(std::move(result));
            }
            
            std::unordered_map<Genotype<Haplotype>, double> germline_genotype_posteriors;
            std::unordered_map<Haplotype, double> cancer_haplotype_posteriors;
            GenotypeModel::Cancer::GenotypeMixtures genotype_mixtures;
        };
        
        mutable GenotypeModel::Cancer genotype_model_;
        
        const SampleIdType normal_sample_;
        const double min_variant_posterior_          = 0.95;
        const double min_somatic_mutation_posterior_ = 0.9;
        const double min_refcall_posterior_          = 0.5;
        
        bool call_somatics_only_;
        
//        std::vector<VcfRecord> call_variants(const GenomicRegion& region,
//                                             const std::vector<Variant>& candidates,
//                                             const ReadMap& reads) const override;
        
        std::unique_ptr<CallerLatents>
        infer_latents(const std::vector<SampleIdType>& samples,
                      const std::vector<Haplotype>& haplotypes,
                      const HaplotypePriorMap& haplotype_priors,
                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const override;
        
        std::vector<VcfRecord::Builder>
        call_variants(const std::vector<Variant>& candidates,
                      const std::vector<Allele>& callable_alleles,
                      CallerLatents* latents,
                      const Phaser::PhaseSet& phase_set,
                      const ReadMap& reads) const override;
    };
    
} // namespace Octopus

#endif /* defined(__Octopus__cancer_caller__) */
