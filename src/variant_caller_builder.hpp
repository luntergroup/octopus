//
//  variant_caller_builder.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef variant_caller_builder_hpp
#define variant_caller_builder_hpp

#include <unordered_map>
#include <string>
#include <memory>
#include <functional>

#include <boost/optional.hpp>

#include "common.hpp"
#include "variant_caller.hpp"
#include "read_pipe.hpp"
#include "candidate_generator_builder.hpp"
#include "haplotype_generator.hpp"

#include "pedigree.hpp"

namespace Octopus {
    
    class VariantCallerBuilder
    {
    public:
        using RefCallType = VariantCaller::RefCallType;
        
        VariantCallerBuilder()  = delete;
        
        explicit VariantCallerBuilder(const ReferenceGenome& reference,
                                      ReadPipe& read_pipe,
                                      const CandidateGeneratorBuilder& candidate_generator_builder,
                                      HaplotypeGenerator::Builder haplotype_generator_builder);
        
        ~VariantCallerBuilder() = default;
        
        VariantCallerBuilder(const VariantCallerBuilder&);
        VariantCallerBuilder& operator=(const VariantCallerBuilder&);
        VariantCallerBuilder(VariantCallerBuilder&&);
        VariantCallerBuilder& operator=(VariantCallerBuilder&&);
        
        // common
        VariantCallerBuilder& set_reference(const ReferenceGenome& reference) noexcept;
        VariantCallerBuilder& set_read_pipe(ReadPipe& read_pipe) noexcept;
        VariantCallerBuilder& set_candidate_generator_builder(const CandidateGeneratorBuilder& generator) noexcept;        
        VariantCallerBuilder& set_ploidy(unsigned ploidy) noexcept;
        VariantCallerBuilder& set_caller(std::string caller);
        VariantCallerBuilder& set_refcall_type(VariantCaller::RefCallType refcall_type) noexcept;
        VariantCallerBuilder& set_sites_only() noexcept;
        VariantCallerBuilder& set_min_variant_posterior(double min_posterior) noexcept;
        VariantCallerBuilder& set_min_refcall_posterior(double min_posterior) noexcept;
        VariantCallerBuilder& set_max_haplotypes(unsigned max_haplotypes) noexcept;
        VariantCallerBuilder& set_min_haplotype_posterior(double p) noexcept;
        VariantCallerBuilder& set_flank_scoring(bool allow_flank_scoring) noexcept;
        VariantCallerBuilder& set_min_phase_score(double min_phase_score) noexcept;
        
        // cancer
        VariantCallerBuilder& set_normal_sample(SampleIdType normal_sample);
        VariantCallerBuilder& set_somatic_mutation_rate(double somatic_mutation_rate);
        VariantCallerBuilder& set_min_somatic_posterior(double min_posterior) noexcept;
        VariantCallerBuilder& set_somatic_only_calls() noexcept;
        VariantCallerBuilder& set_somatic_and_variant_calls() noexcept;
        VariantCallerBuilder& set_somatic_and_variant_and_refcalls_calls() noexcept;
        
        // trio
        
        VariantCallerBuilder& set_maternal_sample(SampleIdType mother);
        VariantCallerBuilder& set_paternal_sample(SampleIdType father);
        
        // pedigree
        
        VariantCallerBuilder& set_pedigree(Pedigree pedigree);
        
        // build
        
        std::unique_ptr<VariantCaller> build() const;
        
    private:
        struct Parameters
        {
            Parameters()  = delete;
            explicit Parameters(const ReferenceGenome& reference,
                                ReadPipe& read_pipe,
                                const CandidateGeneratorBuilder& candidate_generator_builder,
                                HaplotypeGenerator::Builder haplotype_generator_builder);
            ~Parameters() = default;
            
            Parameters(const Parameters&)            = default;
            Parameters& operator=(const Parameters&) = default;
            Parameters(Parameters&&)                 = default;
            Parameters& operator=(Parameters&&)      = default;
            
            // common
            std::reference_wrapper<const ReferenceGenome> reference;
            std::reference_wrapper<ReadPipe> read_pipe;
            
            unsigned ploidy;
            std::string caller;
            std::reference_wrapper<const CandidateGeneratorBuilder> candidate_generator_builder;
            HaplotypeGenerator::Builder haplotype_generator_builder;
            VariantCaller::RefCallType refcall_type = VariantCaller::RefCallType::None;
            bool call_sites_only = false;
            double min_variant_posterior;
            double min_refcall_posterior;
            unsigned max_haplotypes;
            double min_haplotype_posterior;
            bool allow_flank_scoring;
            double min_phase_score;
            
            // cancer
            
            boost::optional<SampleIdType> normal_sample;
            double somatic_mutation_rate;
            double min_somatic_posterior;
            bool call_somatics_only;
            
            // trio
            
            boost::optional<SampleIdType> maternal_sample, paternal_sample;
            
            // pedigree
            
            boost::optional<Pedigree> pedigree;
        };
        
        using CallerFactoryMap = std::unordered_map<std::string, std::function<std::unique_ptr<VariantCaller>()>>;
        
        Parameters parameters_;
        CallerFactoryMap factory_;
        
        CallerFactoryMap generate_factory() const;
    };
} // namespace Octopus

#endif /* variant_caller_builder_hpp */
