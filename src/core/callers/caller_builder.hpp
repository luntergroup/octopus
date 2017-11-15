// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef caller_builder_hpp
#define caller_builder_hpp

#include <unordered_map>
#include <string>
#include <memory>
#include <functional>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/ploidy_map.hpp"
#include "readpipe/read_pipe.hpp"
#include "core/tools/coretools.hpp"
#include "basics/trio.hpp"
#include "basics/pedigree.hpp"
#include "caller.hpp"

namespace octopus {

class CallerBuilder
{
public:
    using RefCallType = Caller::RefCallType;
    
    CallerBuilder() = delete;
    
    CallerBuilder(const ReferenceGenome& reference, const ReadPipe& read_pipe,
                  VariantGeneratorBuilder vgb, HaplotypeGenerator::Builder hgb);
    
    CallerBuilder(const CallerBuilder&);
    CallerBuilder& operator=(const CallerBuilder&);
    CallerBuilder(CallerBuilder&&);
    CallerBuilder& operator=(CallerBuilder&&);
    
    ~CallerBuilder() = default;
    
    // common
    CallerBuilder& set_reference(const ReferenceGenome& reference) noexcept;
    CallerBuilder& set_read_pipe(const ReadPipe& read_pipe) noexcept;
    CallerBuilder& set_variant_generator(const VariantGeneratorBuilder& vb) noexcept;
    CallerBuilder& set_ploidies(PloidyMap ploidies) noexcept;
    CallerBuilder& set_caller(std::string caller);
    CallerBuilder& set_refcall_type(Caller::RefCallType type) noexcept;
    CallerBuilder& set_sites_only() noexcept;
    CallerBuilder& set_min_variant_posterior(Phred<double> posterior) noexcept;
    CallerBuilder& set_min_refcall_posterior(Phred<double> posterior) noexcept;
    CallerBuilder& set_max_haplotypes(unsigned n) noexcept;
    CallerBuilder& set_haplotype_extension_threshold(Phred<double> p) noexcept;
    CallerBuilder& set_flank_scoring(bool b) noexcept;
    CallerBuilder& set_model_filtering(bool b) noexcept;
    CallerBuilder& set_min_phase_score(Phred<double> score) noexcept;
    CallerBuilder& set_snp_heterozygosity(double heterozygosity) noexcept;
    CallerBuilder& set_indel_heterozygosity(double heterozygosity) noexcept;
    CallerBuilder& set_max_joint_genotypes(unsigned max) noexcept;
    CallerBuilder& set_sequencer(std::string sequencer) noexcept;
    CallerBuilder& set_model_mapping_quality(bool b) noexcept;
    
    // cancer
    CallerBuilder& set_normal_sample(SampleName normal_sample);
    CallerBuilder& set_somatic_mutation_rate(double rate) noexcept;
    CallerBuilder& set_min_expected_somatic_frequency(double frequency) noexcept;
    CallerBuilder& set_credible_mass(double mass) noexcept;
    CallerBuilder& set_min_credible_somatic_frequency(double frequency) noexcept;
    CallerBuilder& set_min_somatic_posterior(Phred<double> posterior) noexcept;
    
    // trio
    CallerBuilder& set_trio(Trio trio);
    CallerBuilder& set_min_denovo_posterior(Phred<double> posterior) noexcept;
    CallerBuilder& set_snv_denovo_mutation_rate(double rate) noexcept;
    CallerBuilder& set_indel_denovo_mutation_rate(double rate) noexcept;
    
    // pedigree
    CallerBuilder& set_pedigree(Pedigree pedigree);
    
    std::unique_ptr<Caller> build(const ContigName& contig) const;
    
private:
    struct Components
    {
        std::reference_wrapper<const ReferenceGenome> reference;
        std::reference_wrapper<const ReadPipe> read_pipe;
        VariantGeneratorBuilder variant_generator_builder;
        HaplotypeGenerator::Builder haplotype_generator_builder;
        Phaser phaser;
    };
    
    struct Parameters
    {
        // common
        Caller::Parameters general;
        PloidyMap ploidies;
        Phred<double> min_variant_posterior, min_refcall_posterior;
        boost::optional<double> snp_heterozygosity, indel_heterozygosity;
        Phred<double> min_phase_score;
        unsigned max_joint_genotypes;
        
        // cancer
        boost::optional<SampleName> normal_sample;
        double somatic_mutation_rate;
        double min_expected_somatic_frequency;
        double credible_mass;
        double min_credible_somatic_frequency;
        Phred<double> min_somatic_posterior;
        bool call_somatics_only;
        
        // trio
        boost::optional<Trio> trio;
        Phred<double> min_denovo_posterior;
        boost::optional<double> snv_denovo_mutation_rate, indel_denovo_mutation_rate;
        
        // pedigree
        boost::optional<Pedigree> pedigree;
    };
    
    using CallerFactoryMap = std::unordered_map<std::string, std::function<std::unique_ptr<Caller>()>>;
    
    std::string caller_;
    Components components_;
    Parameters params_;
    CallerFactoryMap factory_;
    
    mutable boost::optional<ContigName> requested_contig_;
    
    Caller::Components make_components() const;
    CallerFactoryMap generate_factory() const;
};

} // namespace octopus

#endif
