//
//  variant_caller_builder.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "variant_caller_builder.hpp"

#include "individual_caller.hpp"
#include "population_caller.hpp"
#include "cancer_caller.hpp"
#include "phaser.hpp"

namespace octopus {

// public methods

VariantCallerBuilder::VariantCallerBuilder(const ReferenceGenome& reference,
                                           const ReadPipe& read_pipe,
                                           VariantGenerator::Builder vgb,
                                           HaplotypeGenerator::Builder hgb)
:
components_ {reference, read_pipe, std::move(vgb), std::move(hgb)},
factory_ {generate_factory()}
{}

VariantCallerBuilder::VariantCallerBuilder(const VariantCallerBuilder& other)
:
caller_ {other.caller_},
components_ {other.components_},
params_ {other.params_},
factory_    {generate_factory()}
{}

VariantCallerBuilder& VariantCallerBuilder::operator=(const VariantCallerBuilder& other)
{
    caller_     = other.caller_;
    components_ = other.components_;
    params_     = other.params_;
    factory_    = generate_factory();
    return *this;
}

VariantCallerBuilder::VariantCallerBuilder(VariantCallerBuilder&& other)
:
caller_ {std::move(other.caller_)},
components_ {std::move(other.components_)},
params_ {std::move(other.params_)},
factory_    {generate_factory()}
{}

VariantCallerBuilder& VariantCallerBuilder::operator=(VariantCallerBuilder&& other)
{
    using std::swap;
    swap(caller_, other.caller_);
    swap(components_, other.components_);
    swap(params_, other.params_);
    factory_ = generate_factory();
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_reference(const ReferenceGenome& reference) noexcept
{
    components_.reference = reference;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_read_pipe(const ReadPipe& read_pipe) noexcept
{
    components_.read_pipe = read_pipe;
    return *this;
}

VariantCallerBuilder&
VariantCallerBuilder::set_variant_generator(const VariantGenerator::Builder& vgb) noexcept
{
    components_.variant_generator_builder = vgb;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_ploidy(unsigned ploidy) noexcept
{
    params_.ploidy = ploidy;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_caller(std::string caller)
{
    caller_ = std::move(caller);
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_refcall_type(VariantCaller::RefCallType type) noexcept
{
    params_.refcall_type = type;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_sites_only() noexcept
{
    params_.call_sites_only = true;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_min_variant_posterior(Phred<double> posterior) noexcept
{
    params_.min_variant_posterior = posterior;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_min_refcall_posterior(Phred<double> posterior) noexcept
{
    params_.min_refcall_posterior = posterior;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_max_haplotypes(unsigned n) noexcept
{
    params_.max_haplotypes = n;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_min_haplotype_posterior(double p) noexcept
{
    params_.min_haplotype_posterior = p;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_flank_scoring(bool b) noexcept
{
    params_.allow_flank_scoring = b;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_model_filtering(bool b) noexcept
{
    params_.allow_model_filtering = b;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_min_phase_score(Phred<double> score) noexcept
{
    params_.min_phase_score = score;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_snp_heterozygosity(double heterozygosity) noexcept
{
    params_.snp_heterozygosity = heterozygosity;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_indel_heterozygosity(double heterozygosity) noexcept
{
    params_.indel_heterozygosity = heterozygosity;
    return *this;
}

// cancer
VariantCallerBuilder& VariantCallerBuilder::set_normal_sample(SampleName normal_sample)
{
    params_.normal_sample = std::move(normal_sample);
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_somatic_mutation_rate(double rate) noexcept
{
    params_.somatic_mutation_rate = rate;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_min_somatic_frequency(double frequency) noexcept
{
    params_.min_somatic_frequency = frequency;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_credible_mass(double mass) noexcept
{
    params_.credible_mass = mass;
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_min_somatic_posterior(Phred<double> posterior) noexcept
{
    params_.min_somatic_posterior = posterior;
    return *this;
}

// trio

VariantCallerBuilder& VariantCallerBuilder::set_maternal_sample(SampleName mother)
{
    params_.maternal_sample = std::move(mother);
    return *this;
}

VariantCallerBuilder& VariantCallerBuilder::set_paternal_sample(SampleName father)
{
    params_.paternal_sample = std::move(father);
    return *this;
}

// pedigree

VariantCallerBuilder& VariantCallerBuilder::set_pedigree(Pedigree pedigree)
{
    params_.pedigree = std::move(pedigree);
    return *this;
}

// build

std::unique_ptr<VariantCaller> VariantCallerBuilder::build() const
{
    if (factory_.count(caller_) == 0) {
        throw std::runtime_error {"VariantCallerBuilder: unknown caller " + caller_};
    }
    return factory_.at(caller_)();
}

// private methods

VariantCaller::Components VariantCallerBuilder::make_components() const
{
    return {
        components_.reference,
        components_.read_pipe,
        components_.variant_generator_builder.build(components_.reference),
        components_.haplotype_generator_builder,
        Phaser {params_.min_phase_score}

    };
}

VariantCallerBuilder::CallerFactoryMap VariantCallerBuilder::generate_factory() const
{
    VariantCaller::Parameters general_parameters {
        params_.refcall_type,
        params_.call_sites_only,
        params_.max_haplotypes,
        params_.min_haplotype_posterior,
        params_.allow_flank_scoring,
        params_.allow_model_filtering
    };
    
    return CallerFactoryMap {
        {"individual", [this, general_parameters = std::move(general_parameters)] () {
            return std::make_unique<IndividualVariantCaller>(make_components(),
                                                             std::move(general_parameters),
                                                             IndividualVariantCaller::Parameters {
                                                                 params_.min_variant_posterior,
                                                                 params_.min_refcall_posterior,
                                                                 params_.ploidy,
                                                                 params_.snp_heterozygosity,
                                                                 params_.indel_heterozygosity
                                                             });
        }},
        {"population", [this, general_parameters = std::move(general_parameters)] () {
            return std::make_unique<PopulationVariantCaller>(make_components(),
                                                             std::move(general_parameters),
                                                             PopulationVariantCaller::Parameters {
                                                                 params_.min_variant_posterior,
                                                                 params_.min_refcall_posterior,
                                                                 params_.ploidy
                                                             });
        }},
        {"cancer", [this, general_parameters = std::move(general_parameters)] () {
            return std::make_unique<CancerVariantCaller>(make_components(),
                                                         std::move(general_parameters),
                                                         CancerVariantCaller::Parameters {
                                                             params_.min_variant_posterior,
                                                             params_.min_somatic_posterior,
                                                             params_.min_refcall_posterior,
                                                             params_.ploidy,
                                                             params_.normal_sample,
                                                             params_.somatic_mutation_rate,
                                                             params_.min_somatic_frequency,
                                                             params_.credible_mass,
                                                             50'000
                                                         });
        }}
    };
}

} // namespace octopus
