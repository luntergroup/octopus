// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "caller_builder.hpp"

#include "core/tools/phaser/phaser.hpp"
#include "individual_caller.hpp"
#include "population_caller.hpp"
#include "trio_caller.hpp"
#include "polyclone_caller.hpp"
#include "cell_caller.hpp"

namespace octopus {

// public methods

CallerBuilder::CallerBuilder(const ReferenceGenome& reference, const ReadPipe& read_pipe,
                             VariantGeneratorBuilder vgb, HaplotypeGenerator::Builder hgb)
: components_ {reference, read_pipe, std::move(vgb), std::move(hgb), HaplotypeLikelihoodModel {}, Phaser {}}
, params_ {}
, factory_ {}
{
    params_.general.refcall_type = Caller::RefCallType::none;
    params_.general.refcall_block_merge_threshold = boost::none;
    params_.general.call_sites_only = false;
    params_.general.model_posterior_policy = Caller::ModelPosteriorPolicy::special;
    params_.general.haplotype_extension_threshold = 1e-10;
    params_.general.saturation_limit = 0.9;
    params_.general.max_haplotypes = 200;
    factory_ = generate_factory();
}

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
    params_.general.refcall_type = type;
    return *this;
}

CallerBuilder& CallerBuilder::set_refcall_merge_block_threshold(Phred<double> threshold) noexcept
{
    params_.general.refcall_block_merge_threshold = threshold;
    return *this;
}

CallerBuilder& CallerBuilder::set_sites_only() noexcept
{
    params_.general.call_sites_only = true;
    return *this;
}

CallerBuilder& CallerBuilder::set_reference_haplotype_protection(bool b) noexcept
{
    params_.general.protect_reference_haplotype = b;
    return *this;
}

CallerBuilder& CallerBuilder::set_target_memory_footprint(MemoryFootprint memory) noexcept
{
    params_.general.target_max_memory = memory;
    return *this;
}

CallerBuilder& CallerBuilder::set_execution_policy(ExecutionPolicy policy) noexcept
{
    params_.general.execution_policy = policy;
    return *this;
}

CallerBuilder& CallerBuilder::set_read_linkage(ReadLinkageType linkage) noexcept
{
    params_.general.read_linkage = linkage;
    return *this;
}

CallerBuilder& CallerBuilder::set_bad_region_detector(BadRegionDetector detector) noexcept
{
    components_.bad_region_detector = std::move(detector);
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
    params_.general.max_haplotypes = n;
    return *this;
}

CallerBuilder& CallerBuilder::set_haplotype_extension_threshold(double p) noexcept
{
    params_.general.haplotype_extension_threshold = p;
    return *this;
}

CallerBuilder& CallerBuilder::set_model_posterior_policy(Caller::ModelPosteriorPolicy policy) noexcept
{
    params_.general.model_posterior_policy = policy;
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

CallerBuilder& CallerBuilder::set_max_genotypes(boost::optional<std::size_t> max) noexcept
{
    params_.max_genotypes = max;
    return *this;
}

CallerBuilder& CallerBuilder::set_max_joint_genotypes(boost::optional<std::size_t> max) noexcept
{
    params_.max_joint_genotypes = max;
    return *this;
}

CallerBuilder& CallerBuilder::set_likelihood_model(HaplotypeLikelihoodModel model) noexcept
{
    components_.likelihood_model = std::move(model);
    return *this;
}

CallerBuilder& CallerBuilder::set_model_based_haplotype_dedup(bool use) noexcept
{
    params_.deduplicate_haplotypes_with_caller_model = use;
    return *this;
}

CallerBuilder& CallerBuilder::set_independent_genotype_prior_flag(bool use_independent) noexcept
{
    params_.use_independent_genotype_priors = use_independent;
    return *this;
}

CallerBuilder& CallerBuilder::set_max_vb_seeds(unsigned n) noexcept
{
    params_.max_vb_seeds = n;
    return *this;
}

// cancer

CallerBuilder& CallerBuilder::add_normal_sample(SampleName normal_sample)
{
    params_.normal_samples.push_back(std::move(normal_sample));
    return *this;
}

CallerBuilder& CallerBuilder::set_normal_sample(SampleName normal_sample)
{
    params_.normal_samples = {std::move(normal_sample)};
    return *this;
}

CallerBuilder& CallerBuilder::set_max_somatic_haplotypes(unsigned n) noexcept
{
    params_.max_somatic_haplotypes = n;
    return *this;
}

CallerBuilder& CallerBuilder::set_somatic_snv_mutation_rate(double rate) noexcept
{
    params_.somatic_snv_mutation_rate = rate;
    return *this;
}

CallerBuilder& CallerBuilder::set_somatic_indel_mutation_rate(double rate) noexcept
{
    params_.somatic_indel_mutation_rate = rate;
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

CallerBuilder& CallerBuilder::set_tumour_germline_concentration(double concentration) noexcept
{
    params_.tumour_germline_concentration = concentration;
    return *this;
}

CallerBuilder& CallerBuilder::set_min_somatic_posterior(Phred<double> posterior) noexcept
{
    params_.min_somatic_posterior = posterior;
    return *this;
}

CallerBuilder& CallerBuilder::set_normal_contamination_risk(NormalContaminationRisk risk) noexcept
{
    params_.normal_contamination_risk = risk;
    return *this;
}

CallerBuilder& CallerBuilder::set_trio(Trio trio)
{
    params_.trio = std::move(trio);
    return *this;
}

CallerBuilder& CallerBuilder::set_min_denovo_posterior(Phred<double> posterior) noexcept
{
    params_.min_denovo_posterior = posterior;
    return *this;
}

CallerBuilder& CallerBuilder::set_snv_denovo_mutation_rate(double rate) noexcept
{
    params_.snv_denovo_mutation_rate = rate;
    return *this;
}

CallerBuilder& CallerBuilder::set_indel_denovo_mutation_rate(double rate) noexcept
{
    params_.indel_denovo_mutation_rate = rate;
    return *this;
}

CallerBuilder& CallerBuilder::set_max_clones(unsigned n) noexcept
{
    params_.max_clones = n;
    return *this;
}

CallerBuilder& CallerBuilder::set_dropout_concentration(double concentration) noexcept
{
    params_.dropout_concentration = concentration;
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
        components_.likelihood_model,
        Phaser {Phaser::Config {Phaser::GenotypeMatchType::exact, params_.min_phase_score}},
        components_.bad_region_detector
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
    const auto& samples = components_.read_pipe.get().samples();
    return CallerFactoryMap {
        {"individual", [this, &samples] () {
            return std::make_unique<IndividualCaller>(make_components(),
                                                      params_.general,
                                                      IndividualCaller::Parameters {
                                                          params_.ploidies.of(samples.front(), *requested_contig_),
                                                          make_individual_prior_model(params_.snp_heterozygosity, params_.indel_heterozygosity),
                                                          params_.min_variant_posterior,
                                                          params_.min_refcall_posterior,
                                                          params_.deduplicate_haplotypes_with_caller_model,
                                                          params_.max_genotypes
                                                      });
        }},
        {"population", [this, &samples] () {
            return std::make_unique<PopulationCaller>(make_components(),
                                                      params_.general,
                                                      PopulationCaller::Parameters {
                                                          params_.min_variant_posterior,
                                                          params_.min_refcall_posterior,
                                                          get_ploidies(samples, *requested_contig_, params_.ploidies),
                                                          make_population_prior_model(params_.snp_heterozygosity, params_.indel_heterozygosity),
                                                          params_.max_joint_genotypes,
                                                          params_.use_independent_genotype_priors,
                                                          params_.deduplicate_haplotypes_with_caller_model
                                                      });
        }},
        {"cancer", [this, &samples] () {
            boost::optional<SampleName> normal {};
            if (!params_.normal_samples.empty()) {
                normal = params_.normal_samples.front();
            }
            CancerCaller::Parameters cancer_params {
                params_.min_variant_posterior,
                params_.min_somatic_posterior,
                params_.min_refcall_posterior,
                params_.ploidies.of(samples.front(), *requested_contig_),
                std::move(normal),
                make_cancer_prior_model(params_.snp_heterozygosity, params_.indel_heterozygosity),
                {params_.somatic_snv_mutation_rate, params_.somatic_indel_mutation_rate},
                params_.min_expected_somatic_frequency,
                params_.credible_mass,
                params_.min_credible_somatic_frequency,
                params_.max_genotypes,
                params_.max_somatic_haplotypes,
                params_.normal_contamination_risk,
                params_.deduplicate_haplotypes_with_caller_model,
                params_.max_vb_seeds
            };
            cancer_params.concentrations.somatic.tumour_germline = params_.tumour_germline_concentration;
            return std::make_unique<CancerCaller>(make_components(), params_.general, std::move(cancer_params));
        }},
        {"trio", [this] () {
            return std::make_unique<TrioCaller>(make_components(),
                                                params_.general,
                                                TrioCaller::Parameters {
                                                    *params_.trio,
                                                    params_.ploidies.of(params_.trio->mother(), *requested_contig_),
                                                    params_.ploidies.of(params_.trio->father(), *requested_contig_),
                                                    params_.ploidies.of(params_.trio->child(), *requested_contig_),
                                                    make_trio_prior_model(params_.snp_heterozygosity, params_.indel_heterozygosity),
                                                    {*params_.snv_denovo_mutation_rate, *params_.indel_denovo_mutation_rate},
                                                    params_.min_variant_posterior,
                                                    params_.min_denovo_posterior,
                                                    params_.min_refcall_posterior,
                                                    params_.max_joint_genotypes,
                                                    params_.deduplicate_haplotypes_with_caller_model
                                                });
        }},
        {"polyclone", [this] () {
            return std::make_unique<PolycloneCaller>(make_components(),
                                                     params_.general,
                                                     PolycloneCaller::Parameters {
                                                         make_individual_prior_model(params_.snp_heterozygosity, params_.indel_heterozygosity),
                                                         params_.min_variant_posterior,
                                                         params_.min_refcall_posterior,
                                                         params_.deduplicate_haplotypes_with_caller_model,
                                                         params_.max_clones,
                                                         params_.max_genotypes
                                                     });
        }},
        {"cell", [this, &samples] () {
            return std::make_unique<CellCaller>(make_components(),
                                                params_.general,
                                                CellCaller::Parameters {
                                                    params_.ploidies.of(samples.front(), *requested_contig_),
                                                    make_individual_prior_model(params_.snp_heterozygosity, params_.indel_heterozygosity),
                                                    params_.min_variant_posterior,
                                                    params_.min_refcall_posterior,
                                                    params_.deduplicate_haplotypes_with_caller_model,
                                                    params_.max_clones,
                                                    params_.max_genotypes,
                                                    params_.max_joint_genotypes,
                                                    params_.dropout_concentration,
                                                    {params_.somatic_snv_mutation_rate, params_.somatic_indel_mutation_rate},
                                                    params_.max_vb_seeds,
                                                    params_.normal_samples
                                                });
        }}
    };
}

} // namespace octopus
