// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef caller_builder_hpp
#define caller_builder_hpp

#include <unordered_map>
#include <string>
#include <memory>
#include <functional>

#include <boost/optional.hpp>

#include <config/common.hpp>
#include <readpipe/read_pipe.hpp>
#include <core/tools/coretools.hpp>
#include <core/types/trio.hpp>
#include <core/types/pedigree.hpp>

#include "caller.hpp"

namespace octopus {

class CallerBuilder
{
public:
    using RefCallType = Caller::RefCallType;
    
    CallerBuilder() = delete;
    
    CallerBuilder(const ReferenceGenome& reference,
                  const ReadPipe& read_pipe,
                  VariantGenerator::Builder vgb,
                  HaplotypeGenerator::Builder hgb);
    
    CallerBuilder(const CallerBuilder&);
    CallerBuilder& operator=(const CallerBuilder&);
    CallerBuilder(CallerBuilder&&);
    CallerBuilder& operator=(CallerBuilder&&);
    
    ~CallerBuilder() = default;
    
    // common
    CallerBuilder& set_reference(const ReferenceGenome& reference) noexcept;
    CallerBuilder& set_read_pipe(const ReadPipe& read_pipe) noexcept;
    CallerBuilder& set_variant_generator(const VariantGenerator::Builder& vb) noexcept;
    CallerBuilder& set_ploidy(unsigned ploidy) noexcept;
    CallerBuilder& set_caller(std::string caller);
    CallerBuilder& set_refcall_type(Caller::RefCallType type) noexcept;
    CallerBuilder& set_sites_only() noexcept;
    CallerBuilder& set_min_variant_posterior(Phred<double> posterior) noexcept;
    CallerBuilder& set_min_refcall_posterior(Phred<double> posterior) noexcept;
    CallerBuilder& set_max_haplotypes(unsigned n) noexcept;
    CallerBuilder& set_min_haplotype_posterior(double p) noexcept;
    CallerBuilder& set_flank_scoring(bool b) noexcept;
    CallerBuilder& set_model_filtering(bool b) noexcept;
    CallerBuilder& set_min_phase_score(Phred<double> score) noexcept;
    
    CallerBuilder& set_snp_heterozygosity(double heterozygosity) noexcept;
    CallerBuilder& set_indel_heterozygosity(double heterozygosity) noexcept;
    
    // cancer
    CallerBuilder& set_normal_sample(SampleName normal_sample);
    CallerBuilder& set_somatic_mutation_rate(double rate) noexcept;
    CallerBuilder& set_credible_mass(double mass) noexcept;
    CallerBuilder& set_min_somatic_frequency(double frequency) noexcept;
    CallerBuilder& set_min_somatic_posterior(Phred<double> posterior) noexcept;
    
    // trio
    
    CallerBuilder& set_maternal_sample(SampleName mother);
    CallerBuilder& set_paternal_sample(SampleName father);
    
    // pedigree
    
    CallerBuilder& set_pedigree(Pedigree pedigree);
    
    std::unique_ptr<Caller> build() const;
    
private:
    struct Components
    {
        std::reference_wrapper<const ReferenceGenome> reference;
        std::reference_wrapper<const ReadPipe> read_pipe;
        VariantGenerator::Builder variant_generator_builder;
        HaplotypeGenerator::Builder haplotype_generator_builder;
        Phaser phaser;
    };
    
    struct Parameters
    {
        // common
        unsigned ploidy;
        
        Caller::RefCallType refcall_type = Caller::RefCallType::None;
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
    
    using CallerFactoryMap = std::unordered_map<std::string, std::function<std::unique_ptr<Caller>()>>;
    
    std::string caller_;
    
    Components components_;
    
    Parameters params_;
    
    CallerFactoryMap factory_;
    
    Caller::Components make_components() const;
    
    CallerFactoryMap generate_factory() const;
};

} // namespace octopus

#endif
