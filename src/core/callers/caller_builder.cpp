// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "caller_builder.hpp"

#include "core/tools/phaser/phaser.hpp"
#include "individual_caller.hpp"
#include "population_caller.hpp"
#include "cancer_caller.hpp"
#include "trio_caller.hpp"

namespace octopus {

// public methods

CallerBuilder::CallerBuilder(const ReferenceGenome& reference, const ReadPipe& read_pipe,
                             VariantGeneratorBuilder vgb, HaplotypeGenerator::Builder hgb)
: components_ {reference, read_pipe, std::move(vgb), std::move(hgb), Phaser {}}
, factory_ {generate_factory()}
{}

CallerBuilder::CallerBuilder(const CallerBuilder& other)
: caller_ {other.caller_}
, components_ {other.components_}
, params_ {other.params_}
, factory_    {generate_factory()}
{}

CallerBuilder& CallerBuilder::operator=(const CallerBuilder& other)
{
    caller_     = other.caller_;
    components_ = other.components_;
    params_     = other.params_;
    factory_    = generate_factory();
    return *this;
}

CallerBuilder::CallerBuilder(CallerBuilder&& other)
: caller_ {std::move(other.caller_)}
, components_ {std::move(other.components_)}
, params_ {std::move(other.params_)}
, factory_    {generate_factory()}
{}

CallerBuilder& CallerBuilder::operator=(CallerBuilder&& other)
{
    using std::swap;
    swap(caller_, other.caller_);
    swap(components_, other.components_);
    swap(params_, other.params_);
    factory_ = generate_factory();
    return *this;
}

CallerBuilder& CallerBuilder::set_reference(const ReferenceGenome& reference) noexcept
{
    components_.reference = reference;
    return *this;
}

CallerBuilder& CallerBuilder::set_read_pipe(const ReadPipe& read_pipe) noexcept
{
    components_.read_pipe = read_pipe;
    return *this;
}

CallerBuilder& CallerBuilder::set_variant_generator(const VariantGeneratorBuilder& vgb) noexcept
{
    components_.variant_generator_builder = vgb;
    return *this;
}

CallerBuilder& CallerBuilder::set_ploidies(PloidyMap ploidies) noexcept
{
    params_.ploidies = std::move(ploidies);
    return *this;
}

CallerBuilder& CallerBuilder::set_caller(std::string caller)
{
    caller_ = std::move(caller);
    return *this;
}

CallerBuilder& CallerBuilder::set_refcall_type(Caller::RefCallType type) noexcept
{
    params_.refcall_type = type;
    return *this;
}

CallerBuilder& CallerBuilder::set_sites_only() noexcept
{
    params_.call_sites_only = true;
    return *this;
}

CallerBuilder& CallerBuilder::set_min_variant_posterior(Phred<double> posterior) noexcept
{
    params_.min_variant_posterior = posterior;
    return *this;
}

CallerBuilder& CallerBuilder::set_min_refcall_posterior(Phred<double> posterior) noexcept
{
    params_.min_refcall_posterior = posterior;
    return *this;
}

CallerBuilder& CallerBuilder::set_max_haplotypes(unsigned n) noexcept
{
    params_.max_haplotypes = n;
    return *this;
}

CallerBuilder& CallerBuilder::set_haplotype_extension_threshold(Phred<double> p) noexcept
{
    params_.haplotype_extension_threshold = p;
    return *this;
}

CallerBuilder& CallerBuilder::set_flank_scoring(bool b) noexcept
{
    params_.allow_flank_scoring = b;
    return *this;
}

CallerBuilder& CallerBuilder::set_model_filtering(bool b) noexcept
{
    params_.allow_model_filtering = b;
    return *this;
}

CallerBuilder& CallerBuilder::set_min_phase_score(Phred<double> score) noexcept
{
    params_.min_phase_score = score;
    return *this;
}

CallerBuilder& CallerBuilder::set_snp_heterozygosity(double heterozygosity) noexcept
{
    params_.snp_heterozygosity = heterozygosity;
    return *this;
}

CallerBuilder& CallerBuilder::set_indel_heterozygosity(double heterozygosity) noexcept
{
    params_.indel_heterozygosity = heterozygosity;
    return *this;
}

// cancer
CallerBuilder& CallerBuilder::set_normal_sample(SampleName normal_sample)
{
    params_.normal_sample = std::move(normal_sample);
    return *this;
}

CallerBuilder& CallerBuilder::set_somatic_mutation_rate(double rate) noexcept
{
    params_.somatic_mutation_rate = rate;
    return *this;
}

CallerBuilder& CallerBuilder::set_min_expected_somatic_frequency(double frequency) noexcept
{
    params_.min_expected_somatic_frequency = frequency;
    return *this;
}

CallerBuilder& CallerBuilder::set_credible_mass(double mass) noexcept
{
    params_.credible_mass = mass;
    return *this;
}

CallerBuilder& CallerBuilder::set_min_credible_somatic_frequency(double frequency) noexcept
{
    params_.min_credible_somatic_frequency = frequency;
    return *this;
}

CallerBuilder& CallerBuilder::set_min_somatic_posterior(Phred<double> posterior) noexcept
{
    params_.min_somatic_posterior = posterior;
    return *this;
}

CallerBuilder& CallerBuilder::set_trio(Trio trio)
{
    params_.trio = std::move(trio);
    return *this;
}

CallerBuilder& CallerBuilder::set_denovo_mutation_rate(double rate) noexcept
{
    params_.denovo_mutation_rate = rate;
    return *this;
}

std::unique_ptr<Caller> CallerBuilder::build(const ContigName& contig) const
{
    if (factory_.count(caller_) == 0) {
        throw std::runtime_error {"CallerBuilder: unknown caller " + caller_};
    }
    requested_contig_ = contig;
    return factory_.at(caller_)();
}

// private methods

Caller::Components CallerBuilder::make_components() const
{
    return {
        components_.reference,
        components_.read_pipe,
        components_.variant_generator_builder.build(components_.reference),
        components_.haplotype_generator_builder,
        Phaser {params_.min_phase_score}

    };
}

auto make_individual_prior_model(boost::optional<double> snp_heterozygosity,
                                 boost::optional<double> indel_heterozygosity)
{
    using T = decltype(IndividualCaller::Parameters::prior_model_params);
    if (snp_heterozygosity && indel_heterozygosity) {
        return T {T::value_type {*snp_heterozygosity, *indel_heterozygosity}};
    } else {
        return T {};
    }
}

auto make_population_prior_model(boost::optional<double> snp_heterozygosity,
                                 boost::optional<double> indel_heterozygosity)
{
    using T = decltype(PopulationCaller::Parameters::prior_model_params);
    if (snp_heterozygosity && indel_heterozygosity) {
        return T {T::value_type {*snp_heterozygosity, *indel_heterozygosity}};
    } else {
        return T {};
    }
}

auto make_cancer_prior_model(boost::optional<double> snp_heterozygosity,
                             boost::optional<double> indel_heterozygosity)
{
    using T = decltype(CancerCaller::Parameters::germline_prior_model_params);
    if (snp_heterozygosity && indel_heterozygosity) {
        return T {T::value_type {*snp_heterozygosity, *indel_heterozygosity}};
    } else {
        return T {};
    }
}

auto make_trio_prior_model(boost::optional<double> snp_heterozygosity,
                           boost::optional<double> indel_heterozygosity)
{
    using T = decltype(TrioCaller::Parameters::germline_prior_model_params);
    if (snp_heterozygosity && indel_heterozygosity) {
        return T {T::value_type {*snp_heterozygosity, *indel_heterozygosity}};
    } else {
        return T {};
    }
}

CallerBuilder::CallerFactoryMap CallerBuilder::generate_factory() const
{
    Caller::Parameters general_parameters {
        params_.refcall_type,
        params_.call_sites_only,
        params_.max_haplotypes,
        params_.haplotype_extension_threshold,
        params_.allow_flank_scoring,
        params_.allow_model_filtering
    };
    const auto& samples = components_.read_pipe.get().samples();
    return CallerFactoryMap {
        {"individual", [this, general_parameters = std::move(general_parameters), &samples] () {
            return std::make_unique<IndividualCaller>(make_components(),
                                                      std::move(general_parameters),
                                                      IndividualCaller::Parameters {
                                                          params_.ploidies.of(samples.front(), *requested_contig_),
                                                          make_individual_prior_model(params_.snp_heterozygosity, params_.indel_heterozygosity),
                                                          params_.min_variant_posterior,
                                                          params_.min_refcall_posterior
                                                      });
        }},
        {"population", [this, general_parameters = std::move(general_parameters), &samples] () {
            return std::make_unique<PopulationCaller>(make_components(),
                                                      std::move(general_parameters),
                                                      PopulationCaller::Parameters {
                                                          params_.min_variant_posterior,
                                                          params_.min_refcall_posterior,
                                                          make_population_prior_model(params_.snp_heterozygosity, params_.indel_heterozygosity),
                                                          get_ploidies(samples, *requested_contig_, params_.ploidies),
                                                      });
        }},
        {"cancer", [this, general_parameters = std::move(general_parameters), &samples] () {
            return std::make_unique<CancerCaller>(make_components(),
                                                  std::move(general_parameters),
                                                  CancerCaller::Parameters {
                                                      params_.min_variant_posterior,
                                                      params_.min_somatic_posterior,
                                                      params_.min_refcall_posterior,
                                                      params_.ploidies.of(samples.front(), *requested_contig_),
                                                      params_.normal_sample,
                                                      make_cancer_prior_model(params_.snp_heterozygosity, params_.indel_heterozygosity),
                                                      {params_.somatic_mutation_rate},
                                                      params_.min_expected_somatic_frequency,
                                                      params_.credible_mass,
                                                      params_.min_credible_somatic_frequency
                                                  });
        }},
        {"trio", [this, general_parameters = std::move(general_parameters)] () {
            return std::make_unique<TrioCaller>(make_components(),
                                                std::move(general_parameters),
                                                TrioCaller::Parameters {
                                                    *params_.trio,
                                                    params_.ploidies.of(params_.trio->mother(), *requested_contig_),
                                                    params_.ploidies.of(params_.trio->father(), *requested_contig_),
                                                    params_.ploidies.of(params_.trio->child(), *requested_contig_),
                                                    make_trio_prior_model(params_.snp_heterozygosity, params_.indel_heterozygosity),
                                                    {*params_.denovo_mutation_rate},
                                                    params_.min_variant_posterior,
                                                    params_.min_refcall_posterior,
                                                });
        }}
    };
}

} // namespace octopus
