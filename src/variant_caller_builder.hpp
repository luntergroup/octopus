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
#include "read_pipe.hpp"
#include "composer.hpp"
#include "haplotype_generator.hpp"
#include "variant_caller.hpp"
#include "trio.hpp"
#include "pedigree.hpp"

namespace octopus
{
    using core::generators::Composer;
    
    class VariantCallerBuilder
    {
    public:
        using RefCallType = VariantCaller::RefCallType;
        
        VariantCallerBuilder()  = delete;
        
        VariantCallerBuilder(const ReferenceGenome& reference, const ReadPipe& read_pipe,
                             const Composer::Builder& candidate_variant_generator_builder,
                             HaplotypeGenerator::Builder haplotype_generator_builder);
        
        ~VariantCallerBuilder() = default;
        
        VariantCallerBuilder(const VariantCallerBuilder&);
        VariantCallerBuilder& operator=(const VariantCallerBuilder&);
        VariantCallerBuilder(VariantCallerBuilder&&);
        VariantCallerBuilder& operator=(VariantCallerBuilder&&);
        
        // common
        VariantCallerBuilder& set_reference(const ReferenceGenome& reference) noexcept;
        VariantCallerBuilder& set_read_pipe(const ReadPipe& read_pipe) noexcept;
        VariantCallerBuilder& set_candidate_variant_generator_builder(const Composer::Builder& generator) noexcept;
        VariantCallerBuilder& set_ploidy(unsigned ploidy) noexcept;
        VariantCallerBuilder& set_caller(std::string caller);
        VariantCallerBuilder& set_refcall_type(VariantCaller::RefCallType type) noexcept;
        VariantCallerBuilder& set_sites_only() noexcept;
        VariantCallerBuilder& set_min_variant_posterior(Phred<double> posterior) noexcept;
        VariantCallerBuilder& set_min_refcall_posterior(Phred<double> posterior) noexcept;
        VariantCallerBuilder& set_max_haplotypes(unsigned n) noexcept;
        VariantCallerBuilder& set_min_haplotype_posterior(double p) noexcept;
        VariantCallerBuilder& set_flank_scoring(bool b) noexcept;
        VariantCallerBuilder& set_model_filtering(bool b) noexcept;
        VariantCallerBuilder& set_min_phase_score(Phred<double> score) noexcept;
        
        VariantCallerBuilder& set_snp_heterozygosity(double heterozygosity) noexcept;
        VariantCallerBuilder& set_indel_heterozygosity(double heterozygosity) noexcept;
        
        // cancer
        VariantCallerBuilder& set_normal_sample(SampleName normal_sample);
        VariantCallerBuilder& set_somatic_mutation_rate(double rate) noexcept;
        VariantCallerBuilder& set_credible_mass(double mass) noexcept;
        VariantCallerBuilder& set_min_somatic_frequency(double frequency) noexcept;
        VariantCallerBuilder& set_min_somatic_posterior(Phred<double> posterior) noexcept;
        
        // trio
        
        VariantCallerBuilder& set_maternal_sample(SampleName mother);
        VariantCallerBuilder& set_paternal_sample(SampleName father);
        
        // pedigree
        
        VariantCallerBuilder& set_pedigree(Pedigree pedigree);
        
        // build
        
        std::unique_ptr<VariantCaller> build() const;
        
    private:
        struct Parameters
        {
            Parameters()  = delete;
            
            Parameters(const ReferenceGenome& reference, const ReadPipe& read_pipe,
                       const Composer::Builder& candidate_variant_generator_builder,
                       HaplotypeGenerator::Builder haplotype_generator_builder);
            
            ~Parameters() = default;
            
            Parameters(const Parameters&)            = default;
            Parameters& operator=(const Parameters&) = default;
            Parameters(Parameters&&)                 = default;
            Parameters& operator=(Parameters&&)      = default;
            
            // common
            std::reference_wrapper<const ReferenceGenome> reference;
            std::reference_wrapper<const ReadPipe> read_pipe;
            
            unsigned ploidy;
            std::string caller;
            std::reference_wrapper<const Composer::Builder> candidate_variant_generator_builder;
            HaplotypeGenerator::Builder haplotype_generator_builder;
            VariantCaller::RefCallType refcall_type = VariantCaller::RefCallType::None;
            bool call_sites_only = false;
            Phred<double> min_variant_posterior;
            Phred<double> min_refcall_posterior;
            unsigned max_haplotypes;
            double min_haplotype_posterior;
            bool allow_flank_scoring;
            bool allow_model_filtering;
            
            double snp_heterozygosity;
            double indel_heterozygosity;
            
            Phred<double> min_phase_score;
            
            // cancer
            
            boost::optional<SampleName> normal_sample;
            double somatic_mutation_rate;
            double min_somatic_frequency;
            double credible_mass;
            Phred<double> min_somatic_posterior;
            bool call_somatics_only;
            
            // trio
            
            boost::optional<SampleName> maternal_sample, paternal_sample;
            
            // pedigree
            
            boost::optional<Pedigree> pedigree;
        };
        
        using CallerFactoryMap = std::unordered_map<std::string, std::function<std::unique_ptr<VariantCaller>()>>;
        
        Parameters parameters_;
        CallerFactoryMap factory_;
        
        CallerFactoryMap generate_factory() const;
    };
} // namespace octopus

#endif /* variant_caller_builder_hpp */
