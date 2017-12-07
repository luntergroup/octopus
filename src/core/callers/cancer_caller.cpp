// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cancer_caller.hpp"

#include <typeinfo>
#include <string>
#include <utility>
#include <algorithm>
#include <numeric>
#include <deque>
#include <unordered_set>
#include <stdexcept>
#include <iostream>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "basics/genomic_region.hpp"
#include "containers/probability_matrix.hpp"
#include "readpipe/read_pipe.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/genotype.hpp"
#include "core/models/genotype/uniform_genotype_prior_model.hpp"
#include "core/models/genotype/coalescent_genotype_prior_model.hpp"
#include "utils/read_stats.hpp"
#include "utils/sequence_utils.hpp"
#include "utils/merge_transform.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"
#include "logging/logging.hpp"
#include "core/types/calls/germline_variant_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "core/types/calls/somatic_call.hpp"

namespace octopus {

// public methods

CancerCaller::CancerCaller(Caller::Components&& components,
                           Caller::Parameters general_parameters,
                           Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {std::move(specific_parameters)}
{
    if (parameters_.ploidy == 0) {
        throw std::logic_error {"CancerCaller: ploidy must be > 0"};
    }
    if (parameters_.max_genotypes == 0) {
        throw std::logic_error {"CancerCaller: max genotypes must be > 0"};
    }
    if (has_normal_sample()) {
        if (std::find(std::cbegin(samples_), std::cend(samples_), normal_sample()) == std::cend(samples_)) {
            throw std::invalid_argument {"CancerCaller: normal sample is not a valid sample"};
        }
    }
    if (parameters_.min_variant_posterior == Phred<double> {0}) {
        logging::WarningLogger wlog {};
        wlog << "Having no germline variant posterior threshold means no somatic variants will be called";
    }
    if (debug_log_) {
        if (has_normal_sample()) {
            stream(*debug_log_) << "Normal sample is " << *parameters_.normal_sample;
        } else {
            *debug_log_ << "There is no normal sample";
        }
    }
}

std::string CancerCaller::do_name() const
{
    return "cancer";
}

CancerCaller::CallTypeSet CancerCaller::do_call_types() const
{
    return {
        std::type_index(typeid(GermlineVariantCall)),
        std::type_index(typeid(SomaticCall))
    };
}

// private methods

bool CancerCaller::has_normal_sample() const noexcept
{
    return static_cast<bool>(parameters_.normal_sample);
}

const SampleName& CancerCaller::normal_sample() const
{
    return *parameters_.normal_sample;
}

std::unique_ptr<CancerCaller::Caller::Latents>
CancerCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                            const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    // Store any intermediate results in Latents for reuse, so the order of model evaluation matters!
    auto result = std::make_unique<Latents>(haplotypes, samples_, get_model_priors());
    generate_germline_genotypes(*result, haplotypes);
    if (debug_log_) stream(*debug_log_) << "There are " << result->germline_genotypes_.size() << " candidate germline genotypes";
    evaluate_germline_model(*result, haplotype_likelihoods);
    evaluate_cnv_model(*result, haplotype_likelihoods);
    generate_cancer_genotypes(*result, haplotype_likelihoods);
    if (debug_log_) stream(*debug_log_) << "There are " << result->cancer_genotypes_.size() << " candidate cancer genotypes";
    if (has_normal_sample()) result->normal_sample_ = std::cref(normal_sample());
    evaluate_tumour_model(*result, haplotype_likelihoods);
    evaluate_noise_model(*result, haplotype_likelihoods);
    return result;
}

boost::optional<double>
CancerCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                        const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                        const Caller::Latents& latents) const
{
    return calculate_model_posterior(haplotypes, haplotype_likelihoods,
                                     dynamic_cast<const Latents&>(latents));
}

static double calculate_model_posterior(const double normal_germline_model_log_evidence,
                                        const double normal_dummy_model_log_evidence)
{
    constexpr double normalModelPrior {0.999};
    constexpr double dummyModelPrior {1.0 - normalModelPrior};
    const auto normal_model_ljp = std::log(normalModelPrior) + normal_germline_model_log_evidence;
    const auto dummy_model_ljp  = std::log(dummyModelPrior) + normal_dummy_model_log_evidence;
    const auto norm = maths::log_sum_exp(normal_model_ljp, dummy_model_ljp);
    return std::exp(normal_model_ljp - norm);
}

static double calculate_model_posterior(const double germline_model_log_evidence,
                                        const double dummy_model_log_evidence,
                                        const double noise_model_log_evidence)
{
    constexpr double normalModelPrior {0.999};
    constexpr double dummyModelPrior {1.0 - normalModelPrior};
    const auto normal_model_ljp = std::log(normalModelPrior) + germline_model_log_evidence;
    const auto dummy_model_ljp  = std::log(dummyModelPrior) + dummy_model_log_evidence;
    const auto noise_model_ljp  = std::log(dummyModelPrior) + noise_model_log_evidence;
    const auto norm = maths::log_sum_exp(normal_model_ljp, std::max(dummy_model_ljp, noise_model_ljp));
    return std::exp(normal_model_ljp - norm);
}

boost::optional<double>
CancerCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                        const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                        const Latents& latents) const
{
    if (has_normal_sample()) {
        assert(latents.germline_model_);
        const auto& germline_model = *latents.germline_model_;
        haplotype_likelihoods.prime(normal_sample());
        GermlineModel::InferredLatents normal_inferences;
        if (latents.normal_germline_inferences_) {
            normal_inferences = *latents.normal_germline_inferences_;
        } else {
            normal_inferences = germline_model.evaluate(latents.germline_genotypes_, haplotype_likelihoods);
        }
        const auto dummy_genotypes = generate_all_genotypes(haplotypes, parameters_.ploidy + 1);
        const auto dummy_inferences = germline_model.evaluate(dummy_genotypes, haplotype_likelihoods);
        auto noise_model_priors = get_normal_noise_model_priors(germline_model.prior_model());
        const CNVModel noise_model {{normal_sample()}, std::move(noise_model_priors)};
        auto noise_inferences = noise_model.evaluate(latents.germline_genotypes_, haplotype_likelihoods);
        return octopus::calculate_model_posterior(normal_inferences.log_evidence,
                                                  dummy_inferences.log_evidence,
                                                  noise_inferences.approx_log_evidence);
    } else {
        // TODO
        return boost::none;
    }
}

auto pool_likelihood(const std::vector<SampleName>& samples,
                     const std::vector<Haplotype>& haplotypes,
                     const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    static const SampleName pooled_sample {"pool"};
    auto result = merge_samples(samples, pooled_sample, haplotypes, haplotype_likelihoods);
    result.prime(pooled_sample);
    return result;
}

void CancerCaller::generate_germline_genotypes(Latents& latents, const std::vector<Haplotype>& haplotypes) const
{
    if (haplotypes.size() < 4) {
        latents.germline_genotypes_ = generate_all_genotypes(haplotypes, parameters_.ploidy);
    } else {
        std::vector<std::vector<unsigned>> germline_genotype_indices {};
        latents.germline_genotypes_ = generate_all_genotypes(haplotypes, parameters_.ploidy, germline_genotype_indices);
        latents.germline_genotype_indices_ = std::move(germline_genotype_indices);
    }
}

namespace {

template <typename Genotype_>
auto zip_cref(const std::vector<Genotype_>& genotypes, const std::vector<double>& probabilities)
{
    using GenotypeReference = std::reference_wrapper<const Genotype_>;
    std::vector<std::pair<GenotypeReference, double>> result {};
    result.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(probabilities), std::back_inserter(result),
                   [] (const auto& g, const auto& p) noexcept { return std::make_pair(std::cref(g), p); });
    return result;
}

template <typename Genotype_>
auto extract_greatest_probability_genotypes(const std::vector<Genotype_>& genotypes,
                                            const std::vector<double>& probabilities,
                                            const std::size_t n,
                                            const boost::optional<double> min_include_probability = boost::none,
                                            const boost::optional<double> max_exclude_probability = boost::none)
{
    assert(genotypes.size() == probabilities.size());
    if (genotypes.size() <= n) return genotypes;
    auto genotype_probabilities = zip_cref(genotypes, probabilities);
    auto last_include_itr = std::next(std::begin(genotype_probabilities), n);
    const auto probability_greater = [] (const auto& lhs, const auto& rhs) noexcept { return lhs.second > rhs.second; };
    std::partial_sort(std::begin(genotype_probabilities), last_include_itr, std::end(genotype_probabilities), probability_greater);
    if (min_include_probability) {
        last_include_itr = std::upper_bound(std::begin(genotype_probabilities), last_include_itr, *min_include_probability,
                                            [] (auto lhs, const auto& rhs) noexcept { return lhs > rhs.second; });
        if (last_include_itr == std::begin(genotype_probabilities)) ++last_include_itr;
    }
    if (max_exclude_probability) {
        last_include_itr = std::partition(last_include_itr, std::end(genotype_probabilities),
                                          [&] (const auto& p) noexcept { return p.second > *max_exclude_probability; });
    }
    std::vector<Genotype_> result {};
    result.reserve(std::distance(std::begin(genotype_probabilities), last_include_itr));
    std::transform(std::begin(genotype_probabilities), last_include_itr, std::back_inserter(result),
                   [] (const auto& p) { return p.first.get(); });
    return result;
}

auto extract_greatest_probability_genotypes(const std::vector<Genotype<Haplotype>>& genotypes,
                                            const std::vector<std::vector<unsigned>>& genotype_indices,
                                            const std::vector<double>& probabilities,
                                            const std::size_t n,
                                            const boost::optional<double> min_include_probability = boost::none,
                                            const boost::optional<double> max_exclude_probability = boost::none)
{
    assert(genotypes.size() == genotype_indices.size());
    using GenotypeReference      = std::reference_wrapper<const Genotype<Haplotype>>;
    using GenotypeIndexReference = std::reference_wrapper<const std::vector<unsigned>>;
    std::vector<std::pair<GenotypeReference, GenotypeIndexReference>> zipped {};
    zipped.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(genotype_indices), std::back_inserter(zipped),
                   [] (const auto& g, const auto& g_idx) { return std::make_pair(std::cref(g), std::cref(g_idx)); });
    auto tmp = extract_greatest_probability_genotypes(zipped, probabilities, n, min_include_probability, max_exclude_probability);
    std::vector<Genotype<Haplotype>> result_genotypes {};
    result_genotypes.reserve(tmp.size());
    std::vector<std::vector<unsigned>> result_indices {};
    result_indices.reserve(tmp.size());
    for (const auto& p : tmp) {
        result_genotypes.push_back(p.first.get());
        result_indices.push_back(p.second.get());
    }
    return std::make_pair(std::move(result_genotypes), std::move(result_indices));
}

} // namespace

void CancerCaller::generate_cancer_genotypes(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    const auto& germline_genotypes = latents.germline_genotypes_;
    const auto num_haplotypes = latents.haplotypes_.get().size();
    const auto num_germline_genotypes = germline_genotypes.size();
    const auto max_possible_cancer_genotypes = num_haplotypes * num_germline_genotypes;
    if (max_possible_cancer_genotypes <= parameters_.max_genotypes) {
        generate_cancer_genotypes(latents, latents.germline_genotypes_);
    } else if (has_normal_sample()) {
        if (has_high_normal_contamination_risk(latents)) {
            generate_cancer_genotypes_with_contaminated_normal(latents, haplotype_likelihoods);
        } else {
            generate_cancer_genotypes_with_clean_normal(latents, haplotype_likelihoods);
        }
    } else {
        generate_cancer_genotypes_with_no_normal(latents, haplotype_likelihoods);
    }
}

void CancerCaller::generate_cancer_genotypes_with_clean_normal(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    const auto& germline_genotypes = latents.germline_genotypes_;
    assert(latents.germline_model_);
    haplotype_likelihoods.prime(normal_sample());
    if (latents.germline_genotype_indices_) {
        latents.normal_germline_inferences_ = latents.germline_model_->evaluate(germline_genotypes,
                                                                                *latents.germline_genotype_indices_,
                                                                                haplotype_likelihoods);
    } else {
        latents.normal_germline_inferences_ = latents.germline_model_->evaluate(germline_genotypes, haplotype_likelihoods);
    }
    const auto& germline_normal_posteriors = latents.normal_germline_inferences_->posteriors.genotype_probabilities;
    const auto max_germline_genotype_bases = parameters_.max_genotypes / latents.haplotypes_.get().size();
    if (latents.germline_genotype_indices_) {
        std::vector<Genotype<Haplotype>> germline_bases;
        std::vector<std::vector<unsigned>> germline_bases_indices;
        std::tie(germline_bases, germline_bases_indices) = extract_greatest_probability_genotypes(germline_genotypes,
                                                                                                  *latents.germline_genotype_indices_,
                                                                                                  germline_normal_posteriors,
                                                                                                  max_germline_genotype_bases,
                                                                                                  1e-100, 1e-2);
        std::vector<std::pair<std::vector<unsigned>, unsigned>> cancer_genotype_indices {};
        latents.cancer_genotypes_ = generate_all_cancer_genotypes(germline_bases, germline_bases_indices,
                                                                  latents.haplotypes_, cancer_genotype_indices);
        latents.cancer_genotype_indices_ = std::move(cancer_genotype_indices);
    } else {
        auto germline_bases = extract_greatest_probability_genotypes(germline_genotypes, germline_normal_posteriors,
                                                                     max_germline_genotype_bases, 1e-100, 1e-2);
        latents.cancer_genotypes_ = generate_all_cancer_genotypes(germline_bases, latents.haplotypes_);
    }
}

void CancerCaller::generate_cancer_genotypes_with_contaminated_normal(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    // TODO
    generate_cancer_genotypes(latents, latents.germline_genotypes_);
}

void CancerCaller::generate_cancer_genotypes_with_no_normal(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    // TODO
    generate_cancer_genotypes(latents, latents.germline_genotypes_);
}

void CancerCaller::generate_cancer_genotypes(Latents& latents, const std::vector<Genotype<Haplotype>>& germline_genotypes) const
{
    if (latents.germline_genotype_indices_) {
        std::vector<std::pair<std::vector<unsigned>, unsigned>> cancer_genotype_indices {};
        latents.cancer_genotypes_ = generate_all_cancer_genotypes(germline_genotypes, *latents.germline_genotype_indices_,
                                                                  latents.haplotypes_, cancer_genotype_indices);
        latents.cancer_genotype_indices_ = std::move(cancer_genotype_indices);
    } else {
        latents.cancer_genotypes_ = generate_all_cancer_genotypes(germline_genotypes, latents.haplotypes_);
    }
}

bool CancerCaller::has_high_normal_contamination_risk(const Latents& latents) const
{
    return parameters_.normal_contamination_risk == Parameters::NormalContaminationRisk::high;
}

void CancerCaller::evaluate_germline_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!(latents.haplotypes_.get().empty() || latents.germline_genotypes_.empty()));
    latents.germline_prior_model_ = make_germline_prior_model(latents.haplotypes_);
    latents.germline_model_ = std::make_unique<GermlineModel>(*latents.germline_prior_model_);
    const auto pooled_likelihoods = pool_likelihood(samples_,  latents.haplotypes_, haplotype_likelihoods);
    if (latents.germline_genotype_indices_) {
        latents.germline_prior_model_->prime(latents.haplotypes_);
        latents.germline_model_inferences_ = latents.germline_model_->evaluate(latents.germline_genotypes_,
                                                                               *latents.germline_genotype_indices_,
                                                                               pooled_likelihoods);
    } else {
        latents.germline_model_inferences_ = latents.germline_model_->evaluate(latents.germline_genotypes_, pooled_likelihoods);
    }
}

void CancerCaller::evaluate_cnv_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!latents.germline_genotypes_.empty() && latents.germline_prior_model_);
    auto cnv_model_priors = get_cnv_model_priors(*latents.germline_prior_model_);
    const CNVModel cnv_model {samples_, cnv_model_priors};
    if (latents.germline_genotype_indices_) {
        latents.cnv_model_inferences_ = cnv_model.evaluate(latents.germline_genotypes_, *latents.germline_genotype_indices_,
                                                           haplotype_likelihoods);
    } else {
        latents.cnv_model_inferences_ = cnv_model.evaluate(latents.germline_genotypes_, haplotype_likelihoods);
    }
}

void CancerCaller::evaluate_tumour_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(latents.germline_prior_model_ && !latents.cancer_genotypes_.empty());
    SomaticMutationModel mutation_model {parameters_.somatic_mutation_model_params};
    latents.cancer_genotype_prior_model_ = CancerGenotypePriorModel {*latents.germline_prior_model_, std::move(mutation_model)};
    auto somatic_model_priors = get_somatic_model_priors(*latents.cancer_genotype_prior_model_);
    const TumourModel somatic_model {samples_, somatic_model_priors};
    if (latents.cancer_genotype_indices_) {
        latents.cancer_genotype_prior_model_->mutation_model().prime(latents.haplotypes_);
        latents.tumour_model_inferences_ = somatic_model.evaluate(latents.cancer_genotypes_, *latents.cancer_genotype_indices_,
                                                                  haplotype_likelihoods);
    } else {
        latents.tumour_model_inferences_ = somatic_model.evaluate(latents.cancer_genotypes_, haplotype_likelihoods);
    }
}

auto get_high_posterior_genotypes(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                  const model::TumourModel::InferredLatents& latents)
{
    return extract_greatest_probability_genotypes(genotypes, latents.posteriors.genotype_probabilities, 10, 1e-3);
}

void CancerCaller::evaluate_noise_model(Latents& latents, const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    if (has_normal_sample()) {
        assert(latents.cancer_genotype_prior_model_);
        auto noise_model_priors = get_noise_model_priors(*latents.cancer_genotype_prior_model_);
        const TumourModel noise_model {samples_, noise_model_priors};
        auto noise_genotypes = get_high_posterior_genotypes(latents.cancer_genotypes_, latents.tumour_model_inferences_);
        latents.noise_model_inferences_ = noise_model.evaluate(noise_genotypes, haplotype_likelihoods);
    }
}

CancerCaller::CNVModel::Priors
CancerCaller::get_cnv_model_priors(const GenotypePriorModel& prior_model) const
{
    using Priors = CNVModel::Priors;
    Priors::GenotypeMixturesDirichletAlphaMap cnv_alphas {};
    cnv_alphas.reserve(samples_.size());
    for (const auto& sample : samples_) {
        if (has_normal_sample() && sample == normal_sample()) {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas(parameters_.ploidy, parameters_.cnv_normal_alpha);
            cnv_alphas.emplace(sample, std::move(sample_alphas));
        } else {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas(parameters_.ploidy, parameters_.cnv_tumour_alpha);
            cnv_alphas.emplace(sample, std::move(sample_alphas));
        }
    }
    return Priors {prior_model, std::move(cnv_alphas)};
}

CancerCaller::TumourModel::Priors
CancerCaller::get_somatic_model_priors(const CancerGenotypePriorModel& prior_model) const
{
    using Priors = TumourModel::Priors;
    Priors::GenotypeMixturesDirichletAlphaMap alphas {};
    alphas.reserve(samples_.size());
    for (const auto& sample : samples_) {
        if (has_normal_sample() && sample == normal_sample()) {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas(parameters_.ploidy + 1, parameters_.somatic_normal_germline_alpha);
            sample_alphas.back() = parameters_.somatic_normal_somatic_alpha;
            alphas.emplace(sample, std::move(sample_alphas));
        } else {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas(parameters_.ploidy + 1, parameters_.somatic_tumour_germline_alpha);
            sample_alphas.back() = parameters_.somatic_tumour_somatic_alpha;
            alphas.emplace(sample, std::move(sample_alphas));
        }
    }
    return Priors {prior_model, std::move(alphas)};
}

CancerCaller::TumourModel::Priors
CancerCaller::get_noise_model_priors(const CancerGenotypePriorModel& prior_model) const
{
    // The noise model is intended to capture noise that may also be present in the normal sample,
    // hence all samples have the same prior alphas.
    using Priors = TumourModel::Priors;
    Priors::GenotypeMixturesDirichletAlphas noise_alphas(parameters_.ploidy + 1, parameters_.somatic_tumour_germline_alpha);
    noise_alphas.back() = parameters_.somatic_tumour_somatic_alpha;
    Priors::GenotypeMixturesDirichletAlphaMap alphas {};
    alphas.reserve(samples_.size());
    for (const auto& sample : samples_) {
        alphas.emplace(sample, noise_alphas);
    }
    return Priors {prior_model, std::move(alphas)};
}

CancerCaller::CNVModel::Priors
CancerCaller::get_normal_noise_model_priors(const GenotypePriorModel& prior_model) const
{
    using Priors = CNVModel::Priors;
    Priors::GenotypeMixturesDirichletAlphaMap cnv_alphas {};
    if (has_normal_sample()) {
        Priors::GenotypeMixturesDirichletAlphas sample_alphas(parameters_.ploidy, 0.5);
        cnv_alphas.emplace(normal_sample(), std::move(sample_alphas));
    }
    return Priors {prior_model, std::move(cnv_alphas)};
}

std::vector<std::unique_ptr<VariantCall>>
CancerCaller::call_variants(const std::vector<Variant>& candidates,
                                   const Caller::Latents& latents) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

namespace {

using VariantReference  = std::reference_wrapper<const Variant>;
using VariantPosteriorVector = std::vector<std::pair<VariantReference, Phred<double>>>;

struct GermlineVariantCall : Mappable<GermlineVariantCall>
{
    GermlineVariantCall() = delete;
    GermlineVariantCall(const std::pair<VariantReference, Phred<double>>& p)
    : variant {p.first}
    , posterior {p.second}
    {}
    GermlineVariantCall(const Variant& variant, Phred<double> posterior)
    : variant {variant}
    , posterior {posterior}
    {}
    
    const GenomicRegion& mapped_region() const noexcept { return octopus::mapped_region(variant.get()); }
    
    VariantReference variant;
    Phred<double> posterior;
};

using GermlineVariantCalls = std::vector<GermlineVariantCall>;

struct SomaticVariantCall : Mappable<SomaticVariantCall>
{
    SomaticVariantCall() = delete;
    SomaticVariantCall(const std::pair<VariantReference, Phred<double>>& p)
    : variant {p.first}, posterior {p.second} {}
    SomaticVariantCall(const Variant& variant, Phred<double> posterior)
    : variant {variant}, posterior {posterior} {}
    
    const GenomicRegion& mapped_region() const noexcept { return octopus::mapped_region(variant.get()); }
    
    VariantReference variant;
    Phred<double> posterior;
};

using SomaticVariantCalls = std::vector<SomaticVariantCall>;

struct GermlineGenotypeCall
{
    template <typename T>
    GermlineGenotypeCall(T&& genotype, Phred<double> posterior)
    : genotype {std::forward<T>(genotype)}
    , somatic {}
    , posterior {posterior}
    {}
    template <typename T, typename A>
    GermlineGenotypeCall(T&& genotype, A&& somatic, Phred<double> posterior)
    : genotype {std::forward<T>(genotype)}
    , somatic {std::forward<A>(somatic)}
    , posterior {posterior}
    {}
    
    Genotype<Allele> genotype;
    boost::optional<Allele> somatic;
    Phred<double> posterior;
};

using GermlineGenotypeCalls = std::vector<GermlineGenotypeCall>;

struct CancerGenotypeCall
{
    template <typename T>
    CancerGenotypeCall(T&& genotype, Phred<double> posterior)
    : genotype {std::forward<T>(genotype)}, posterior {posterior} {}
    
    CancerGenotype<Allele> genotype;
    Phred<double> posterior;
    std::unordered_map<SampleName, std::vector<std::pair<double, double>>> credible_regions;
};

using CancerGenotypeCalls = std::vector<CancerGenotypeCall>;

template <typename L>
auto find_map_genotype(const L& posteriors)
{
    return std::max_element(std::cbegin(posteriors), std::cend(posteriors),
                            [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
}

// germline variant posterior calculations

template <typename M>
Phred<double> marginalise(const Allele& allele, const M& genotype_posteriors)
{
    auto p = std::accumulate(std::cbegin(genotype_posteriors),
                             std::cend(genotype_posteriors),
                             0.0, [&allele] (const auto curr, const auto& p) {
                                 return curr + (contains(p.first, allele) ? 0.0 : p.second);
                             });
    return probability_to_phred(p);
}

template <typename M>
VariantPosteriorVector compute_candidate_posteriors(const std::vector<Variant>& candidates,
                                                    const M& genotype_posteriors)
{
    VariantPosteriorVector result {};
    result.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        result.emplace_back(candidate, marginalise(candidate.alt_allele(), genotype_posteriors));
    }
    return result;
}

// germline variant calling

bool contains_alt(const Genotype<Haplotype>& genotype_call, const VariantReference& candidate)
{
    return includes(genotype_call, candidate.get().alt_allele());
}

auto call_candidates(const VariantPosteriorVector& candidate_posteriors,
                     const Genotype<Haplotype>& genotype_call,
                     const Phred<double> min_posterior)
{
    GermlineVariantCalls calls {};
    calls.reserve(candidate_posteriors.size());
    std::vector<VariantReference> uncalled {};
    for (const auto& p : candidate_posteriors) {
        if (p.second >= min_posterior) {
            if (contains_alt(genotype_call, p.first)) {
                calls.emplace_back(p.first, p.second);
            }
        } else {
            uncalled.emplace_back(p.first);
        }
    }
    return std::make_pair(std::move(calls), std::move(uncalled));
}

// somatic variant posterior

auto compute_somatic_variant_posteriors(const std::vector<VariantReference>& candidates,
                                        const std::vector<CancerGenotype<Haplotype>>& cancer_genotypes,
                                        const std::vector<double>& cancer_genotype_posteriors,
                                        const double somatic_posterior,
                                        const double somatic_model_posterior)
{
    VariantPosteriorVector result {};
    result.reserve(candidates.size());
    
    for (const auto& candidate : candidates) {
        const auto& allele = candidate.get().alt_allele();
        const auto p = std::inner_product(std::cbegin(cancer_genotypes), std::cend(cancer_genotypes),
                                          std::cbegin(cancer_genotype_posteriors), 0.0, std::plus<> {},
                                          [&allele] (const auto& genotype, auto posterior) {
                                              if (genotype.somatic_element().contains(allele)
                                                  && !contains(genotype.germline_genotype(), allele)) {
                                                  return posterior;
                                              } else {
                                                  return 0.0;
                                              }
                                          });
        const auto complement = std::min(somatic_model_posterior * p * somatic_posterior, 1.0);
        result.emplace_back(candidate, probability_to_phred(1.0 - complement));
    }
    
    return result;
}

auto call_somatic_variants(const VariantPosteriorVector& somatic_variant_posteriors,
                           const CancerGenotype<Haplotype>& called_genotype,
                           const Phred<double> min_posterior)
{
    SomaticVariantCalls result {};
    result.reserve(somatic_variant_posteriors.size());
    std::copy_if(std::begin(somatic_variant_posteriors), std::end(somatic_variant_posteriors), std::back_inserter(result),
                 [min_posterior, &called_genotype] (const auto& p) {
                     return p.second >= min_posterior && includes(called_genotype, p.first.get().alt_allele());
                 });
    return result;
}

template <typename T>
auto compute_marginal_credible_interval(const T& alphas, const double mass)
{
    const auto a0 = std::accumulate(std::cbegin(alphas), std::cend(alphas), 0.0);
    std::vector<std::pair<double, double>> result {};
    result.reserve(alphas.size());
    for (const auto& alpha : alphas) {
        result.push_back(maths::beta_hdi(alpha, a0 - alpha, mass));
    }
    return result;
}

using CredibleRegionMap = std::unordered_map<SampleName, std::vector<std::pair<double, double>>>;

template <typename M>
auto compute_marginal_credible_intervals(const M& alphas, const double mass)
{
    CredibleRegionMap result {};
    result.reserve(alphas.size());
    for (const auto& p : alphas) {
        result.emplace(p.first, compute_marginal_credible_interval(p.second, mass));
    }
    return result;
}

template <typename T>
auto compute_somatic_mass(const T& alphas, const double c = 0.05)
{
    const auto a0 = std::accumulate(std::cbegin(alphas), std::cend(alphas), 0.0);
    return maths::beta_cdf_complement(alphas.back(), a0 - alphas.back(), c);
}

template <typename T>
auto call_somatic_genotypes(const CancerGenotype<Haplotype>& called_genotype,
                            const std::vector<GenomicRegion>& called_somatic_regions,
                            const std::vector<CancerGenotype<Haplotype>>& genotypes,
                            const std::vector<double>& genotype_posteriors,
                            const T& credible_regions)
{
    CancerGenotypeCalls result {};
    result.reserve(called_somatic_regions.size());
    for (const auto& region : called_somatic_regions) {
        auto genotype_chunk = copy<Allele>(called_genotype, region);
        const auto inv_posterior = std::inner_product(std::cbegin(genotypes), std::cend(genotypes),
                                                      std::cbegin(genotype_posteriors), 0.0, std::plus<> {},
                                                      [&genotype_chunk] (const auto& g, auto p) {
                                                          return contains(g, genotype_chunk) ? 0.0 : p;
                                                      });
        result.emplace_back(std::move(genotype_chunk), probability_to_phred(inv_posterior));
        result.back().credible_regions = credible_regions;
    }
    return result;
}

// output

octopus::VariantCall::GenotypeCall convert(GermlineGenotypeCall call)
{
    return octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
}

std::unique_ptr<octopus::VariantCall>
transform_germline_call(GermlineVariantCall&& variant_call, GermlineGenotypeCall&& genotype_call,
                        const std::vector<SampleName>& samples,
                        const std::vector<SampleName>& somatic_samples)
{
    std::vector<std::pair<SampleName, Call::GenotypeCall>> genotypes {};
    for (const auto& sample : samples) {
        if (std::find(std::cbegin(somatic_samples), std::cend(somatic_samples), sample) == std::cend(somatic_samples)) {
            genotypes.emplace_back(sample, convert(genotype_call));
        } else {
            auto copy = genotype_call;
            copy.genotype.emplace(*copy.somatic);
            genotypes.emplace_back(sample, convert(std::move(copy)));
        }
    }
    return std::make_unique<octopus::GermlineVariantCall>(variant_call.variant.get(),
                                                          std::move(genotypes),
                                                          variant_call.posterior);
}

auto transform_somatic_calls(SomaticVariantCalls&& somatic_calls, CancerGenotypeCalls&& genotype_calls,
                             const std::vector<SampleName>& somatic_samples)
{
    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
    result.reserve(somatic_calls.size());
    std::transform(std::make_move_iterator(std::begin(somatic_calls)),
                   std::make_move_iterator(std::end(somatic_calls)),
                   std::make_move_iterator(std::begin(genotype_calls)),
                   std::back_inserter(result),
                   [&somatic_samples] (auto&& variant_call, auto&& genotype_call) -> std::unique_ptr<octopus::VariantCall> {
                       std::unordered_map<SampleName, SomaticCall::GenotypeCredibleRegions> credible_regions {};
                       for (const auto& p : genotype_call.credible_regions) {
                           SomaticCall::GenotypeCredibleRegions sample_credible_regions {};
                           sample_credible_regions.germline.reserve(p.second.size() - 1);
                           std::copy(std::cbegin(p.second), std::prev(std::cend(p.second)),
                                     std::back_inserter(sample_credible_regions.germline));
                           if (std::find(std::cbegin(somatic_samples), std::cend(somatic_samples), p.first) != std::cend(somatic_samples)) {
                               sample_credible_regions.somatic = p.second.back();
                           }
                           credible_regions.emplace(p.first, std::move(sample_credible_regions));
                       }
                       return std::make_unique<SomaticCall>(variant_call.variant.get(),
                                                            std::move(genotype_call.genotype),
                                                            genotype_call.posterior,
                                                            std::move(credible_regions),
                                                            variant_call.posterior);
                   });
    return result;
}

} // namespace

std::vector<std::unique_ptr<VariantCall>>
CancerCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    // TODO: refactor this into smaller methods!
    const auto model_posteriors = calculate_model_posteriors(latents);
    log(model_posteriors);
    const auto germline_genotype_posteriors = calculate_germline_genotype_posteriors(latents, model_posteriors);
    const auto& cancer_genotype_posteriors = latents.tumour_model_inferences_.posteriors.genotype_probabilities;
    log(latents.germline_genotypes_, germline_genotype_posteriors, latents.cnv_model_inferences_,
        latents.cancer_genotypes_, latents.tumour_model_inferences_);
    const auto sample_somatic_inv_posteriors = calculate_probability_samples_not_somatic(latents);
    const auto somatic_posterior = calculate_somatic_probability(sample_somatic_inv_posteriors, model_posteriors);
    const auto germline_candidate_posteriors = compute_candidate_posteriors(candidates, germline_genotype_posteriors);
    Genotype<Haplotype> called_germline_genotype {};
    boost::optional<CancerGenotype<Haplotype>> called_cancer_genotype {};
    if (model_posteriors.somatic > model_posteriors.germline && somatic_posterior >= parameters_.min_somatic_posterior) {
        if (debug_log_) *debug_log_ << "Using cancer genotype for germline genotype call";
        if (!called_cancer_genotype) {
            auto cancer_posteriors = zip_cref(latents.cancer_genotypes_, cancer_genotype_posteriors);
            called_cancer_genotype = find_map_genotype(cancer_posteriors)->first;
        }
        called_germline_genotype = called_cancer_genotype->germline_genotype();
    } else {
        called_germline_genotype = find_map_genotype(germline_genotype_posteriors)->first;
    }
    GermlineVariantCalls germline_variant_calls;
    std::vector<VariantReference> uncalled_germline_candidates;
    std::tie(germline_variant_calls, uncalled_germline_candidates) = call_candidates(germline_candidate_posteriors,
                                                                                     called_germline_genotype,
                                                                                     parameters_.min_variant_posterior);
    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
    boost::optional<Haplotype> called_somatic_haplotype {};
    std::vector<SampleName> somatic_samples {};
    if (somatic_posterior >= parameters_.min_somatic_posterior) {
        auto somatic_allele_posteriors = compute_somatic_variant_posteriors(uncalled_germline_candidates,
                                                                            latents.cancer_genotypes_,
                                                                            cancer_genotype_posteriors,
                                                                            somatic_posterior.probability_true(),
                                                                            model_posteriors.somatic);
        if (!called_cancer_genotype) {
            auto cancer_posteriors = zip_cref(latents.cancer_genotypes_, cancer_genotype_posteriors);
            called_cancer_genotype = find_map_genotype(cancer_posteriors)->first.get();
        }
        if (called_cancer_genotype->germline_genotype() == called_germline_genotype) {
            auto somatic_variant_calls = call_somatic_variants(somatic_allele_posteriors, *called_cancer_genotype,
                                                               parameters_.min_somatic_posterior);
            const auto& somatic_alphas = latents.tumour_model_inferences_.posteriors.alphas;
            const auto credible_regions = compute_marginal_credible_intervals(somatic_alphas, parameters_.credible_mass);
            if (!somatic_variant_calls.empty()) {
                for (const auto& p : credible_regions) {
                    if (p.second.back().first >= parameters_.min_credible_somatic_frequency) {
                        if (has_normal_sample() && p.first == normal_sample()) {
                            somatic_samples.clear();
                            break;
                        }
                        somatic_samples.push_back(p.first);
                    }
                }
                if (has_normal_sample() && latents.noise_model_inferences_) {
                    // Does the normal sample contain the called somatic variant?
                    const auto& noisy_alphas = latents.noise_model_inferences_->posteriors.alphas.at(normal_sample());
                    const auto noise_credible_region = compute_marginal_credible_interval(noisy_alphas, parameters_.credible_mass).back();
                    const auto somatic_mass = compute_somatic_mass(noisy_alphas, parameters_.min_expected_somatic_frequency);
                    if (noise_credible_region.first >= parameters_.min_credible_somatic_frequency || somatic_mass > 0.5) {
                        somatic_samples.clear();
                    }
                }
                if (somatic_samples.empty()) {
                    somatic_variant_calls.clear();
                    somatic_variant_calls.shrink_to_fit();
                } else {
                    called_somatic_haplotype = called_cancer_genotype->somatic_element();
                }
            }
            const auto called_somatic_regions = extract_regions(somatic_variant_calls);
            auto cancer_genotype_calls = call_somatic_genotypes(*called_cancer_genotype, called_somatic_regions,
                                                                latents.cancer_genotypes_, cancer_genotype_posteriors,
                                                                credible_regions);
            result = transform_somatic_calls(std::move(somatic_variant_calls), std::move(cancer_genotype_calls),
                                             somatic_samples);
        } else if (debug_log_) {
            stream(*debug_log_) << "Conflict between called germline genotype and called cancer genotype. Not calling somatics";
        }
    }
    
    const auto called_germline_regions = extract_regions(germline_variant_calls);
    GermlineGenotypeCalls germline_genotype_calls {};
    germline_genotype_calls.reserve(called_germline_regions.size());
    for (const auto& region : called_germline_regions) {
        auto genotype_chunk = copy<Allele>(called_germline_genotype, region);
        const auto inv_posterior = std::accumulate(std::cbegin(germline_genotype_posteriors),
                                                   std::cend(germline_genotype_posteriors), 0.0,
                                                   [&called_germline_genotype] (const double curr, const auto& p) {
                                                       return curr + (contains(p.first, called_germline_genotype) ? 0.0 : p.second);
                                                   });
        if (called_somatic_haplotype) {
            germline_genotype_calls.emplace_back(std::move(genotype_chunk),
                                                 copy<Allele>(*called_somatic_haplotype, region),
                                                 probability_to_phred(inv_posterior));
        } else {
            germline_genotype_calls.emplace_back(std::move(genotype_chunk), probability_to_phred(inv_posterior));
        }
    }
    result.reserve(result.size() + germline_variant_calls.size());
    const auto itr = std::end(result);
    std::transform(std::make_move_iterator(std::begin(germline_variant_calls)),
                   std::make_move_iterator(std::end(germline_variant_calls)),
                   std::make_move_iterator(std::begin(germline_genotype_calls)),
                   std::back_inserter(result),
                   [this, &somatic_samples] (auto&& variant_call, auto&& genotype_call) {
                       return transform_germline_call(std::move(variant_call), std::move(genotype_call),
                                                      samples_, somatic_samples);
                   });
    std::inplace_merge(std::begin(result), itr, std::end(result),
                       [] (const auto& lhs, const auto& rhs) { return *lhs < *rhs; });
    return result;
}

CancerCaller::ModelPriors CancerCaller::get_model_priors() const
{
    const auto s = parameters_.germline_weight + parameters_.cnv_weight + parameters_.somatic_weight;
    return {parameters_.germline_weight / s, parameters_.cnv_weight / s, parameters_.somatic_weight / s};
}

CancerCaller::ModelPosteriors
CancerCaller::calculate_model_posteriors(const Latents& latents) const
{
    const auto& germline_inferences = latents.germline_model_inferences_;
    const auto& cnv_inferences      = latents.cnv_model_inferences_;
    const auto& somatic_inferences  = latents.tumour_model_inferences_;
    const auto& model_priors        = latents.model_priors_;
    const auto germline_model_jlp = std::log(model_priors.germline) + germline_inferences.log_evidence;
    const auto cnv_model_jlp      = std::log(model_priors.cnv) + cnv_inferences.approx_log_evidence;
    const auto somatic_model_jlp  = std::log(model_priors.somatic) + somatic_inferences.approx_log_evidence;
    const auto norm = maths::log_sum_exp(germline_model_jlp, cnv_model_jlp, somatic_model_jlp);
    auto germline_model_posterior = std::exp(germline_model_jlp - norm);
    auto cnv_model_posterior      = std::exp(cnv_model_jlp - norm);
    auto somatic_model_posterior  = std::exp(somatic_model_jlp - norm);
    return {germline_model_posterior, cnv_model_posterior, somatic_model_posterior};
}

CancerCaller::GermlineGenotypeProbabilityMap
CancerCaller::calculate_germline_genotype_posteriors(const Latents& latents,
                                                     const ModelPosteriors& model_posteriors) const
{
    const auto& germline_genotypes = latents.germline_genotypes_;
    GermlineGenotypeProbabilityMap result {germline_genotypes.size()};
    
    std::transform(std::cbegin(germline_genotypes), std::cend(germline_genotypes),
                   std::cbegin(latents.germline_model_inferences_.posteriors.genotype_probabilities),
                   std::inserter(result, std::begin(result)),
                   [&model_posteriors] (const auto& genotype, const auto& posterior) {
                       return std::make_pair(genotype, model_posteriors.germline * posterior);
                   });
    const auto& cnv_posteriors = latents.cnv_model_inferences_.posteriors.genotype_probabilities;
    for (std::size_t i {0}; i < latents.germline_genotypes_.size(); ++i) {
        result[germline_genotypes[i]] += model_posteriors.cnv * cnv_posteriors[i];
    }
    const auto& cancer_genotypes = latents.cancer_genotypes_;
    const auto& tumour_posteriors = latents.tumour_model_inferences_.posteriors.genotype_probabilities;
    for (std::size_t i {0}; i < cancer_genotypes.size(); ++i) {
        result[cancer_genotypes[i].germline_genotype()] += model_posteriors.somatic * tumour_posteriors[i];
    }
    
    return result;
}

CancerCaller::ProbabilityVector
CancerCaller::calculate_probability_samples_not_somatic(const Latents& inferences) const
{
    std::vector<double> result(samples_.size());
    const auto& posterior_alphas = inferences.tumour_model_inferences_.posteriors.alphas;
    std::transform(std::cbegin(posterior_alphas), std::cend(posterior_alphas),
                   std::begin(result), [this] (const auto& p) {
                       const auto a0 = std::accumulate(std::cbegin(p.second), std::prev(std::cend(p.second)), 0.0);
                       return maths::beta_cdf(p.second.back(), a0, parameters_.min_expected_somatic_frequency);
                   });
    return result;
}

Phred<double> CancerCaller::calculate_somatic_probability(const ProbabilityVector& sample_somatic_posteriors,
                                                          const ModelPosteriors& model_posteriors) const
{
    auto result = 1.0 - std::accumulate(std::cbegin(sample_somatic_posteriors),
                                        std::cend(sample_somatic_posteriors),
                                        1.0, std::multiplies<> {});
    result *= model_posteriors.somatic;
    return probability_to_phred(1 - result);
}

std::vector<std::unique_ptr<ReferenceCall>>
CancerCaller::call_reference(const std::vector<Allele>& alleles,
                             const Caller::Latents& latents,
                             const ReadMap& reads) const
{
    return {};
}

std::unique_ptr<GenotypePriorModel> CancerCaller::make_germline_prior_model(const std::vector<Haplotype>& haplotypes) const
{
    if (parameters_.germline_prior_model_params) {
        return std::make_unique<CoalescentGenotypePriorModel>(CoalescentModel {
        Haplotype {octopus::mapped_region(haplotypes.front()), reference_},
        *parameters_.germline_prior_model_params
        });
    } else {
        return std::make_unique<UniformGenotypePriorModel>();
    }
}

// CancerCaller::Latents

CancerCaller::Latents::Latents(const std::vector<Haplotype>& haplotypes,
                               const std::vector<SampleName>& samples,
                               CancerCaller::ModelPriors model_priors)
: haplotypes_ {haplotypes}
, samples_ {samples}
, model_priors_ {model_priors}
{}

std::shared_ptr<CancerCaller::Latents::HaplotypeProbabilityMap>
CancerCaller::Latents::haplotype_posteriors() const
{
    if (haplotype_posteriors_ == nullptr) {
        compute_haplotype_posteriors();
    }
    return haplotype_posteriors_;
}

std::shared_ptr<CancerCaller::Latents::GenotypeProbabilityMap>
CancerCaller::Latents::genotype_posteriors() const
{
    if (genotype_posteriors_ == nullptr) {
        compute_genotype_posteriors();
    }
    return genotype_posteriors_;
}

void CancerCaller::Latents::compute_genotype_posteriors() const
{
    // TODO: properly
    GenotypeProbabilityMap genotype_posteriors {std::begin(germline_genotypes_), std::end(germline_genotypes_)};
    for (const auto& sample : samples_.get()) {
        insert_sample(sample, germline_model_inferences_.posteriors.genotype_probabilities, genotype_posteriors);
    }
    genotype_posteriors_ = std::make_shared<Latents::GenotypeProbabilityMap>(std::move(genotype_posteriors));
}

namespace {

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
    auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end   = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}

}

void CancerCaller::Latents::compute_haplotype_posteriors() const
{
    Latents::HaplotypeProbabilityMap result {haplotypes_.get().size()};
    // Contribution from germline model
    for (const auto& haplotype : haplotypes_.get()) {
        result.emplace(haplotype, 0.0);
    }
    for (const auto& p : zip(germline_genotypes_, germline_model_inferences_.posteriors.genotype_probabilities)) {
        for (const auto& haplotype : p.get<0>().copy_unique_ref()) {
            result.at(haplotype) += p.get<1>();
        }
    }
    // Contribution from CNV model
    Latents::HaplotypeProbabilityMap cnv_result {haplotypes_.get().size()};
    for (const auto& haplotype : haplotypes_.get()) {
        cnv_result.emplace(haplotype, 0.0);
    }
    for (const auto& p : zip(germline_genotypes_, cnv_model_inferences_.posteriors.genotype_probabilities)) {
        for (const auto& haplotype : p.get<0>().copy_unique_ref()) {
            cnv_result.at(haplotype) += p.get<1>();
        }
    }
    // Contribution from tumour model
    Latents::HaplotypeProbabilityMap somatic_result {haplotypes_.get().size()};
    for (const auto& haplotype : haplotypes_.get()) {
        somatic_result.emplace(haplotype, 0.0);
    }
    for (const auto& p : zip(cancer_genotypes_, tumour_model_inferences_.posteriors.genotype_probabilities)) {
        for (const auto& haplotype : p.get<0>().germline_genotype().copy_unique_ref()) {
            somatic_result.at(haplotype) += p.get<1>();
        }
        somatic_result.at(p.get<0>().somatic_element()) += p.get<1>();
    }
    for (auto& p : result) {
        p.second *= model_priors_.germline;
        p.second += model_priors_.cnv * cnv_result.at(p.first);
        p.second += model_priors_.somatic * somatic_result.at(p.first);
    }
    haplotype_posteriors_ = std::make_shared<Latents::HaplotypeProbabilityMap>(std::move(result));
}

// logging

void CancerCaller::log(const ModelPosteriors& model_posteriors) const
{
    if (debug_log_) {
        stream(*debug_log_) << "Germline model posterior: " << model_posteriors.germline;
        stream(*debug_log_) << "CNV model posterior:      " << model_posteriors.cnv;
        stream(*debug_log_) << "Somatic model posterior:  " << model_posteriors.somatic;
    }
}

void CancerCaller::log(const GenotypeVector& germline_genotypes,
                       const GermlineGenotypeProbabilityMap& germline_genotype_posteriors,
                       const CNVModel::InferredLatents& cnv_inferences,
                       const CancerGenotypeVector& cancer_genotypes,
                       const TumourModel::InferredLatents& tumour_inferences) const
{
    if (debug_log_) {
        auto map_germline = find_map_genotype(germline_genotype_posteriors);
        auto germline_log = stream(*debug_log_);
        germline_log << "MAP germline genotype: ";
        debug::print_variant_alleles(germline_log, map_germline->first);
        auto cnv_posteriors = zip_cref(germline_genotypes, cnv_inferences.posteriors.genotype_probabilities);
        auto map_cnv = find_map_genotype(cnv_posteriors);
        auto cnv_log = stream(*debug_log_);
        cnv_log << "MAP CNV genotype: ";
        debug::print_variant_alleles(cnv_log, map_cnv->first);
        auto somatic_log = stream(*debug_log_);
        auto cancer_posteriors = zip_cref(cancer_genotypes, tumour_inferences.posteriors.genotype_probabilities);
        auto map_somatic = find_map_genotype(cancer_posteriors);
        auto map_cancer_genotype = map_somatic->first.get();
        somatic_log << "MAP cancer genotype: ";
        debug::print_variant_alleles(somatic_log, map_cancer_genotype);
        somatic_log << ' ' << map_somatic->second;
    }
}

} // namespace octopus
