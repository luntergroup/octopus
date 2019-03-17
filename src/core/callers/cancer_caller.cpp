// Copyright (c) 2015-2019 Daniel Cooke
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
#include <boost/multiprecision/gmp.hpp>

#include "basics/genomic_region.hpp"
#include "containers/probability_matrix.hpp"
#include "readpipe/read_pipe.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/genotype.hpp"
#include "core/models/genotype/uniform_genotype_prior_model.hpp"
#include "core/models/genotype/coalescent_genotype_prior_model.hpp"
#include "core/models/genotype/constant_mixture_genotype_likelihood_model.hpp"
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
    if (parameters_.concentrations.cnv.normal <= 0.0
        || parameters_.concentrations.cnv.tumour <= 0.0
        || parameters_.concentrations.somatic.normal_germline <= 0.0
        || parameters_.concentrations.somatic.normal_somatic <= 0.0
        || parameters_.concentrations.somatic.tumour_germline <= 0.0
        || parameters_.concentrations.somatic.tumour_somatic <= 0.0) {
        throw std::invalid_argument {"CancerCaller: concentration parameters must be positive"};
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
    if (!has_normal_sample()) {
        parameters_.concentrations.cnv.tumour = parameters_.concentrations.somatic.tumour_germline;
    }
}

// private methods

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

unsigned CancerCaller::do_min_callable_ploidy() const
{
    return parameters_.ploidy;
}

unsigned CancerCaller::do_max_callable_ploidy() const
{
    return parameters_.ploidy + parameters_.max_somatic_haplotypes;
}

bool CancerCaller::has_normal_sample() const noexcept
{
    return static_cast<bool>(parameters_.normal_sample);
}

const SampleName& CancerCaller::normal_sample() const
{
    return *parameters_.normal_sample;
}

std::size_t CancerCaller::do_remove_duplicates(std::vector<Haplotype>& haplotypes) const
{
    if (parameters_.deduplicate_haplotypes_with_germline_model) {
        if (haplotypes.size() < 2) return 0;
        CoalescentModel::Parameters model_params {};
        if (parameters_.germline_prior_model_params) model_params = *parameters_.germline_prior_model_params;
        Haplotype reference {mapped_region(haplotypes.front()), reference_.get()};
        CoalescentModel model {std::move(reference), model_params, haplotypes.size(), CoalescentModel::CachingStrategy::none};
        const CoalescentProbabilityGreater cmp {std::move(model)};
        return octopus::remove_duplicates(haplotypes, cmp);
    } else {
        return Caller::do_remove_duplicates(haplotypes);
    }
}

std::unique_ptr<CancerCaller::Caller::Latents>
CancerCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                            const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    // Store any intermediate results in Latents for reuse, so the order of model evaluation matters!
    auto result = std::make_unique<Latents>(haplotypes, samples_, parameters_);
    set_model_priors(*result);
    generate_germline_genotypes(*result, haplotypes);
    if (debug_log_) stream(*debug_log_) << "There are " << result->germline_genotypes_.size() << " candidate germline genotypes";
    evaluate_germline_model(*result, haplotype_likelihoods);
    evaluate_cnv_model(*result, haplotype_likelihoods);
    if (haplotypes.size() > 1) {
        fit_somatic_model(*result, haplotype_likelihoods);
        evaluate_noise_model(*result, haplotype_likelihoods);
        set_model_posteriors(*result);
    }
    return result;
}

boost::optional<double>
CancerCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                        const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                        const Caller::Latents& latents) const
{
    return calculate_model_posterior(haplotypes, haplotype_likelihoods,
                                     dynamic_cast<const Latents&>(latents));
}

void CancerCaller::set_cancer_genotype_prior_model(Latents& latents) const
{
    SomaticMutationModel mutation_model {parameters_.somatic_mutation_model_params};
    latents.cancer_genotype_prior_model_ = CancerGenotypePriorModel {*latents.germline_prior_model_, std::move(mutation_model)};
}

void CancerCaller::fit_somatic_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    set_cancer_genotype_prior_model(latents);
    SomaticModel::InferredLatents prev_latents;
    std::vector<CancerGenotype<Haplotype>> prev_cancer_genotypes;
    boost::optional<std::vector<CancerGenotypeIndex>> prev_cancer_genotype_indices;
    for (unsigned somatic_ploidy {1}; somatic_ploidy <= parameters_.max_somatic_haplotypes; ++somatic_ploidy) {
        if (debug_log_) stream(*debug_log_) << "Fitting somatic model with somatic ploidy " << somatic_ploidy;
        latents.somatic_ploidy_ = somatic_ploidy;
        generate_cancer_genotypes(latents, haplotype_likelihoods);
        if (debug_log_) stream(*debug_log_) << "There are " << latents.cancer_genotypes_.size() << " candidate cancer genotypes";
        evaluate_somatic_model(latents, haplotype_likelihoods);
        if (somatic_ploidy > 1) {
            if (latents.somatic_model_inferences_.approx_log_evidence <= prev_latents.approx_log_evidence) {
                break;
            }
        } else {
            set_model_posteriors(latents);
            if (latents.model_posteriors_.somatic < std::max(latents.model_posteriors_.germline, latents.model_posteriors_.cnv)) {
                break;
            }
        }
        if (latents.haplotypes_.get().size() <= somatic_ploidy + 1) break;
        if (somatic_ploidy < parameters_.max_somatic_haplotypes) {
            // save previous state, but don't move as next cancer genotype generation may use this information
            prev_latents = latents.somatic_model_inferences_;
            prev_cancer_genotypes = latents.cancer_genotypes_;
            prev_cancer_genotype_indices = latents.cancer_genotype_indices_;
        }
    }
    if (latents.somatic_ploidy_ > 1) {
        if (latents.cancer_genotypes_.empty()
            || latents.somatic_model_inferences_.approx_log_evidence <= prev_latents.approx_log_evidence) {
            // load previous state
            --latents.somatic_ploidy_;
            latents.somatic_model_inferences_ = std::move(prev_latents);
            latents.cancer_genotypes_ = std::move(prev_cancer_genotypes);
            latents.cancer_genotype_indices_ = std::move(prev_cancer_genotype_indices);
        }
    }
    if (debug_log_) stream(*debug_log_) << "Using somatic model with somatic ploidy " << latents.somatic_ploidy_;
}

static double calculate_model_posterior(const double normal_germline_model_log_evidence,
                                        const double normal_dummy_model_log_evidence)
{
    constexpr double normalModelPrior {0.99};
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
    constexpr double normalModelPrior {0.99};
    constexpr double dummyModelPrior {1.0 - normalModelPrior};
    const auto normal_model_ljp = std::log(normalModelPrior) + germline_model_log_evidence;
    const auto dummy_model_ljp  = std::log(dummyModelPrior) + dummy_model_log_evidence;
    const auto noise_model_ljp  = std::log(dummyModelPrior) + noise_model_log_evidence;
    const auto norm = maths::log_sum_exp(normal_model_ljp, std::max(dummy_model_ljp, noise_model_ljp));
    return std::exp(normal_model_ljp - norm);
}

namespace {

auto demote_each(const std::vector<CancerGenotype<Haplotype>>& genotypes)
{
    std::vector<Genotype<Haplotype>> result {};
    result.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::back_inserter(result),
                   [] (const auto& genotype) { return demote(genotype); });
    return result;
}

} // namespace

boost::optional<double>
CancerCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                        const HaplotypeLikelihoodArray& haplotype_likelihoods,
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
        const auto dummy_genotypes = demote_each(latents.cancer_genotypes_);
        const auto dummy_inferences = germline_model.evaluate(dummy_genotypes, haplotype_likelihoods);
        if (latents.noise_model_inferences_) {
            return octopus::calculate_model_posterior(normal_inferences.log_evidence,
                                                      dummy_inferences.log_evidence,
                                                      latents.noise_model_inferences_->approx_log_evidence);
        } else {
            return octopus::calculate_model_posterior(normal_inferences.log_evidence,
                                                      dummy_inferences.log_evidence);
        }
    } else {
        // TODO
        return boost::none;
    }
}

auto pool_likelihood(const std::vector<SampleName>& samples,
                     const std::vector<Haplotype>& haplotypes,
                     const HaplotypeLikelihoodArray& haplotype_likelihoods)
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
        std::vector<GenotypeIndex> germline_genotype_indices {};
        latents.germline_genotypes_ = generate_all_genotypes(haplotypes, parameters_.ploidy, germline_genotype_indices);
        latents.germline_genotype_indices_ = std::move(germline_genotype_indices);
    }
}

namespace {

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
    auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end   = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}

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

template <typename T>
auto copy_greatest_probability_values(const std::vector<T>& values,
                                      const std::vector<double>& probabilities,
                                      const std::size_t n,
                                      const boost::optional<double> min_include_probability = boost::none,
                                      const boost::optional<double> max_exclude_probability = boost::none)
{
    assert(values.size() == probabilities.size());
    if (values.size() <= n) return values;
    auto value_probabilities = zip_cref(values, probabilities);
    auto last_include_itr = std::next(std::begin(value_probabilities), n);
    const auto probability_greater = [] (const auto& lhs, const auto& rhs) noexcept { return lhs.second > rhs.second; };
    std::partial_sort(std::begin(value_probabilities), last_include_itr, std::end(value_probabilities), probability_greater);
    if (min_include_probability) {
        last_include_itr = std::upper_bound(std::begin(value_probabilities), last_include_itr, *min_include_probability,
                                            [] (auto lhs, const auto& rhs) noexcept { return lhs > rhs.second; });
        if (last_include_itr == std::begin(value_probabilities)) ++last_include_itr;
    }
    if (max_exclude_probability) {
        last_include_itr = std::partition(last_include_itr, std::end(value_probabilities),
                                          [&] (const auto& p) noexcept { return p.second > *max_exclude_probability; });
    }
    std::vector<T> result {};
    result.reserve(std::distance(std::begin(value_probabilities), last_include_itr));
    std::transform(std::begin(value_probabilities), last_include_itr, std::back_inserter(result),
                   [] (const auto& p) { return p.first.get(); });
    return result;
}

template <typename G, typename I>
auto copy_greatest_probability_genotypes(const std::vector<G>& genotypes,
                                         const std::vector<I>& genotype_indices,
                                         const std::vector<double>& probabilities,
                                         const std::size_t n,
                                         const boost::optional<double> min_include_probability = boost::none,
                                         const boost::optional<double> max_exclude_probability = boost::none)
{
    assert(genotypes.size() == genotype_indices.size());
    using GenotypeReference = std::reference_wrapper<const G>;
    using GenotypeIndexReference = std::reference_wrapper<const I>;
    std::vector<std::pair<GenotypeReference, GenotypeIndexReference>> zipped {};
    zipped.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(genotype_indices), std::back_inserter(zipped),
                   [] (const auto& g, const auto& g_idx) { return std::make_pair(std::cref(g), std::cref(g_idx)); });
    auto tmp = copy_greatest_probability_values(zipped, probabilities, n, min_include_probability, max_exclude_probability);
    std::vector<G> result_genotypes {};
    result_genotypes.reserve(tmp.size());
    std::vector<I> result_indices {};
    result_indices.reserve(tmp.size());
    for (const auto& p : tmp) {
        result_genotypes.push_back(p.first.get());
        result_indices.push_back(p.second.get());
    }
    return std::make_pair(std::move(result_genotypes), std::move(result_indices));
}

auto calculate_posteriors_with_germline_likelihood_model(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                                         const std::vector<CancerGenotypeIndex>& indices,
                                                         const CancerGenotypePriorModel& prior_model,
                                                         const model::ConstantMixtureGenotypeLikelihoodModel likelihood_model,
                                                         const std::vector<SampleName>& samples)
{
    auto result = evaluate(indices, prior_model);
    GenotypeIndex flattened_index(genotypes.front().ploidy());
    for (const auto& sample : samples) {
        likelihood_model.cache().prime(sample);
        std::transform(std::cbegin(indices), std::cend(indices), std::cbegin(result), std::begin(result),
                       [&] (const auto& genotype, auto curr) {
                           auto itr = std::copy(std::cbegin(genotype.germline), std::cbegin(genotype.germline), std::begin(flattened_index));
                           std::copy(std::cbegin(genotype.somatic), std::cbegin(genotype.somatic), itr);
                           return curr + likelihood_model.evaluate(flattened_index);
                       });
    }
    maths::normalise_exp(result);
    return result;
}

void filter_with_germline_model(std::vector<CancerGenotype<Haplotype>>& genotypes,
                                std::vector<CancerGenotypeIndex>& indices,
                                const CancerGenotypePriorModel& prior_model,
                                const model::ConstantMixtureGenotypeLikelihoodModel likelihood_model,
                                const std::vector<SampleName>& samples,
                                const std::size_t n)
{
    const auto germline_model_posteriors = calculate_posteriors_with_germline_likelihood_model(genotypes, indices, prior_model, likelihood_model, samples);
    auto result = copy_greatest_probability_genotypes(genotypes, indices, germline_model_posteriors, n);
    genotypes = std::move(result.first);
    indices = std::move(result.second);
}

} // namespace

void CancerCaller::generate_cancer_genotypes(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto num_haplotypes = latents.haplotypes_.get().size();
    if (num_haplotypes == 1) return;
    const auto& germline_genotypes = latents.germline_genotypes_;
    const auto num_germline_genotypes = germline_genotypes.size();
    const auto max_possible_cancer_genotypes = num_haplotypes * num_germline_genotypes;
    const auto max_allowed_cancer_genotypes = std::max(parameters_.max_genotypes, num_germline_genotypes);
    if (max_possible_cancer_genotypes <= max_allowed_cancer_genotypes) {
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

auto calculate_max_germline_genotype_bases(const unsigned max_genotypes, const unsigned num_haplotypes,
                                           const unsigned somatic_ploidy)
{
    const auto num_somatic_genotypes = num_genotypes(num_haplotypes, somatic_ploidy);
    return std::max(max_genotypes / num_somatic_genotypes, decltype(num_somatic_genotypes) {1});
}

void CancerCaller::generate_cancer_genotypes_with_clean_normal(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto& haplotypes = latents.haplotypes_.get();
    const auto& germline_genotypes = latents.germline_genotypes_;
    const auto max_allowed_cancer_genotypes = std::max(parameters_.max_genotypes, germline_genotypes.size());
    if (!latents.cancer_genotypes_.empty()) {
        const auto max_old_cancer_genotype_bases = std::max(max_allowed_cancer_genotypes / haplotypes.size(), std::size_t {1});
        const auto& cancer_genotype_posteriors = latents.somatic_model_inferences_.posteriors.genotype_probabilities;
        if (latents.cancer_genotype_indices_) {
            std::vector<CancerGenotype<Haplotype>> old_cancer_genotype_bases {};
            std::vector<CancerGenotypeIndex> old_cancer_genotype_index_bases {};
            std::tie(old_cancer_genotype_bases, old_cancer_genotype_index_bases)
            = copy_greatest_probability_genotypes(latents.cancer_genotypes_, *latents.cancer_genotype_indices_,
                                                  cancer_genotype_posteriors, max_old_cancer_genotype_bases);
            latents.cancer_genotypes_ = extend_somatic_genotypes(old_cancer_genotype_bases, old_cancer_genotype_index_bases,
                                                                 haplotypes, *latents.cancer_genotype_indices_);
        } else {
            const auto old_cancer_genotype_bases = copy_greatest_probability_values(latents.cancer_genotypes_, cancer_genotype_posteriors, max_old_cancer_genotype_bases);
            latents.cancer_genotypes_ = extend_somatic_genotypes(old_cancer_genotype_bases, haplotypes);
        }
    } else {
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
        const auto max_germline_genotype_bases = calculate_max_germline_genotype_bases(max_allowed_cancer_genotypes, haplotypes.size(), latents.somatic_ploidy_);
        if (latents.germline_genotype_indices_) {
            std::vector<Genotype<Haplotype>> germline_bases;
            std::vector<GenotypeIndex> germline_bases_indices;
            std::tie(germline_bases, germline_bases_indices) = copy_greatest_probability_genotypes(germline_genotypes,
                                                                                                   *latents.germline_genotype_indices_,
                                                                                                   germline_normal_posteriors,
                                                                                                   max_germline_genotype_bases,
                                                                                                   1e-100, 1e-2);
            std::vector<CancerGenotypeIndex> cancer_genotype_indices {};
            latents.cancer_genotypes_ = generate_all_cancer_genotypes(germline_bases, germline_bases_indices,
                                                                      haplotypes, cancer_genotype_indices,
                                                                      latents.somatic_ploidy_);
            if (latents.cancer_genotypes_.size() > 2 * max_allowed_cancer_genotypes) {
                if (!latents.cancer_genotype_prior_model_->mutation_model().is_primed()) {
                    latents.cancer_genotype_prior_model_->mutation_model().prime(haplotypes);
                }
                const model::ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods, haplotypes};
                filter_with_germline_model(latents.cancer_genotypes_, cancer_genotype_indices, *latents.cancer_genotype_prior_model_,
                                           likelihood_model, samples_, max_allowed_cancer_genotypes);
            }
            latents.cancer_genotype_indices_ = std::move(cancer_genotype_indices);
        } else {
            auto germline_bases = copy_greatest_probability_values(germline_genotypes, germline_normal_posteriors,
                                                                   max_germline_genotype_bases, 1e-100, 1e-2);
            latents.cancer_genotypes_ = generate_all_cancer_genotypes(germline_bases, haplotypes, latents.somatic_ploidy_);
        }
    }
}

void CancerCaller::generate_cancer_genotypes_with_contaminated_normal(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    // TODO
    generate_cancer_genotypes_with_clean_normal(latents, haplotype_likelihoods);
}

namespace {

struct GenotypeReferenceEqual
{
    using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
    std::size_t operator()(const GenotypeReference& lhs, const GenotypeReference& rhs) const
    {
        return lhs.get() == rhs.get();
    }
};

template <typename BidirIt, typename T, typename Compare>
BidirIt binary_find(BidirIt first, BidirIt last, const T& value, Compare cmp)
{
    const auto itr = std::lower_bound(first, last, value, std::move(cmp));
    return (itr != last && *itr == value) ? itr : last;
}

} // namespace

void CancerCaller::generate_cancer_genotypes_with_no_normal(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto& haplotypes = latents.haplotypes_.get();
    const auto& germline_genotypes = latents.germline_genotypes_;
    const auto max_allowed_cancer_genotypes = std::max(parameters_.max_genotypes, germline_genotypes.size());
    
    if (!latents.cancer_genotypes_.empty()) {
        const auto max_old_cancer_genotype_bases = std::max(max_allowed_cancer_genotypes / haplotypes.size(), std::size_t {1});
        const auto& cancer_genotype_posteriors = latents.somatic_model_inferences_.posteriors.genotype_probabilities;
        if (latents.cancer_genotype_indices_) {
            std::vector<CancerGenotype<Haplotype>> old_cancer_genotype_bases {};
            std::vector<CancerGenotypeIndex> old_cancer_genotype_index_bases {};
            std::tie(old_cancer_genotype_bases, old_cancer_genotype_index_bases)
                = copy_greatest_probability_genotypes(latents.cancer_genotypes_, *latents.cancer_genotype_indices_,
                                                      cancer_genotype_posteriors, max_old_cancer_genotype_bases);
            latents.cancer_genotypes_ = extend_somatic_genotypes(old_cancer_genotype_bases, old_cancer_genotype_index_bases,
                                                                haplotypes, *latents.cancer_genotype_indices_);
        } else {
            const auto old_cancer_genotype_bases = copy_greatest_probability_values(latents.cancer_genotypes_, cancer_genotype_posteriors, max_old_cancer_genotype_bases);
            latents.cancer_genotypes_ = extend_somatic_genotypes(old_cancer_genotype_bases, haplotypes);
        }
    } else {
        const auto max_germline_genotype_bases = calculate_max_germline_genotype_bases(max_allowed_cancer_genotypes, haplotypes.size(), latents.somatic_ploidy_);
        const auto& germline_genotype_posteriors = latents.germline_model_inferences_.posteriors.genotype_probabilities;
        std::vector<double> germline_model_haplotype_posteriors(haplotypes.size());
        if (latents.germline_genotype_indices_) {
            GenotypeIndex buffer {};
            for (std::size_t g {0}; g < germline_genotypes.size(); ++g) {
                const auto& g_indices = (*latents.germline_genotype_indices_)[g];
                for (auto idx : g_indices) {
                    if (std::find(std::cbegin(buffer), std::cend(buffer), idx) == std::cend(buffer)) {
                        germline_model_haplotype_posteriors[idx] += germline_genotype_posteriors[g];
                    }
                }
                buffer.clear();
            }
        } else {
            std::unordered_map<HaplotypeReference, double> tmp {};
            tmp.reserve(haplotypes.size());
            for (std::size_t g {0}; g < germline_genotypes.size(); ++g) {
                for (const auto& haplotype : germline_genotypes[g].copy_unique_ref()) {
                    tmp[haplotype] += germline_genotype_posteriors[g];
                }
                std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::begin(germline_model_haplotype_posteriors),
                               [&tmp] (const auto& haplotype) { return tmp.at(haplotype); });
            }
        }
        const auto max_germline_haplotype_bases = max_num_elements(max_germline_genotype_bases, parameters_.ploidy);
        const auto top_haplotypes = copy_greatest_probability_values(haplotypes, germline_model_haplotype_posteriors,
                                                                     max_germline_haplotype_bases);
        auto germline_bases = generate_all_genotypes(top_haplotypes, parameters_.ploidy);
        if (latents.germline_genotype_indices_) {
            std::vector<GenotypeIndex> germline_bases_indices;
            germline_bases_indices.reserve(germline_bases.size());
            if (std::is_sorted(std::cbegin(germline_genotypes), std::cend(germline_genotypes), GenotypeLess {})) {
                std::sort(std::begin(germline_bases), std::end(germline_bases), GenotypeLess {});
                auto genotype_itr = std::cbegin(germline_genotypes);
                for (const auto& genotype : germline_bases) {
                    const auto match_itr = binary_find(genotype_itr, std::cend(germline_genotypes), genotype, GenotypeLess {});
                    assert(match_itr != std::cend(germline_genotypes));
                    const auto idx = std::distance(std::cbegin(germline_genotypes), match_itr);
                    germline_bases_indices.push_back((*latents.germline_genotype_indices_)[idx]);
                    genotype_itr = std::next(match_itr);
                }
            } else {
                using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
                using GenotypeReferenceIndexMap = std::unordered_map<GenotypeReference, std::size_t,
                                                                     std::hash<GenotypeReference>, GenotypeReferenceEqual>;
                GenotypeReferenceIndexMap genotype_indices {};
                genotype_indices.reserve(germline_genotypes.size());
                for (std::size_t i {0}; i < germline_genotypes.size(); ++i) {
                    genotype_indices.emplace(std::cref(germline_genotypes[i]), i);
                }
                for (const auto& genotype : germline_bases) {
                    germline_bases_indices.push_back((*latents.germline_genotype_indices_)[genotype_indices.at(genotype)]);
                }
            }
            std::vector<CancerGenotypeIndex> cancer_genotype_indices {};
            latents.cancer_genotypes_ = generate_all_cancer_genotypes(germline_bases, germline_bases_indices,
                                                                      latents.haplotypes_, cancer_genotype_indices,
                                                                      latents.somatic_ploidy_);
            if (latents.cancer_genotypes_.size() > 2 * max_allowed_cancer_genotypes) {
                if (!latents.cancer_genotype_prior_model_->mutation_model().is_primed()) {
                    latents.cancer_genotype_prior_model_->mutation_model().prime(latents.haplotypes_);
                }
                const model::ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods, latents.haplotypes_};
                filter_with_germline_model(latents.cancer_genotypes_, cancer_genotype_indices, *latents.cancer_genotype_prior_model_,
                                           likelihood_model, samples_, max_allowed_cancer_genotypes);
            }
            latents.cancer_genotype_indices_ = std::move(cancer_genotype_indices);
        } else {
            latents.cancer_genotypes_ = generate_all_cancer_genotypes(germline_bases, haplotypes, latents.somatic_ploidy_);
        }
    }
}

void CancerCaller::generate_cancer_genotypes(Latents& latents, const std::vector<Genotype<Haplotype>>& germline_genotypes) const
{
    if (latents.germline_genotype_indices_) {
        std::vector<CancerGenotypeIndex> cancer_genotype_indices {};
        latents.cancer_genotypes_ = generate_all_cancer_genotypes(germline_genotypes, *latents.germline_genotype_indices_,
                                                                  latents.haplotypes_, cancer_genotype_indices,
                                                                  latents.somatic_ploidy_);
        latents.cancer_genotype_indices_ = std::move(cancer_genotype_indices);
    } else {
        latents.cancer_genotypes_ = generate_all_cancer_genotypes(germline_genotypes, latents.haplotypes_, latents.somatic_ploidy_);
    }
}

bool CancerCaller::has_high_normal_contamination_risk(const Latents& latents) const
{
    return parameters_.normal_contamination_risk == Parameters::NormalContaminationRisk::high;
}

void CancerCaller::evaluate_germline_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(!(latents.haplotypes_.get().empty() || latents.germline_genotypes_.empty()));
    latents.germline_prior_model_ = make_germline_prior_model(latents.haplotypes_);
    latents.germline_model_ = std::make_unique<GermlineModel>(*latents.germline_prior_model_);
    const auto pooled_likelihoods = pool_likelihood(samples_,  latents.haplotypes_, haplotype_likelihoods);
    if (latents.germline_genotype_indices_) {
        latents.germline_prior_model_->prime(latents.haplotypes_);
        latents.germline_model_->prime(latents.haplotypes_);
        latents.germline_model_inferences_ = latents.germline_model_->evaluate(latents.germline_genotypes_,
                                                                               *latents.germline_genotype_indices_,
                                                                               pooled_likelihoods);
    } else {
        latents.germline_model_inferences_ = latents.germline_model_->evaluate(latents.germline_genotypes_, pooled_likelihoods);
    }
}

void CancerCaller::evaluate_cnv_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(!latents.germline_genotypes_.empty() && latents.germline_prior_model_);
    auto cnv_model_priors = get_cnv_model_priors(*latents.germline_prior_model_);
    CNVModel::AlgorithmParameters params {};
    if (parameters_.max_vb_seeds) params.max_seeds = *parameters_.max_vb_seeds;
    params.target_max_memory = this->target_max_memory();
    CNVModel cnv_model {samples_, cnv_model_priors};
    CNVModel cnv_model {samples_, cnv_model_priors, params};
    if (latents.germline_genotype_indices_) {
        cnv_model.prime(latents.haplotypes_);
        latents.cnv_model_inferences_ = cnv_model.evaluate(latents.germline_genotypes_, *latents.germline_genotype_indices_,
                                                           haplotype_likelihoods);
    } else {
        latents.cnv_model_inferences_ = cnv_model.evaluate(latents.germline_genotypes_, haplotype_likelihoods);
    }
}

void CancerCaller::evaluate_somatic_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(latents.germline_prior_model_ && !latents.cancer_genotypes_.empty());
    assert(latents.cancer_genotype_prior_model_);
    auto somatic_model_priors = get_somatic_model_priors(*latents.cancer_genotype_prior_model_, latents.somatic_ploidy_);
    SomaticModel::AlgorithmParameters params {};
    if (parameters_.max_vb_seeds) params.max_seeds = *parameters_.max_vb_seeds;
    params.target_max_memory = this->target_max_memory();
    SomaticModel model {samples_, somatic_model_priors, params};
    if (latents.cancer_genotype_indices_) {
        assert(latents.cancer_genotype_prior_model_->germline_model().is_primed());
        if (!latents.cancer_genotype_prior_model_->mutation_model().is_primed()) {
            latents.cancer_genotype_prior_model_->mutation_model().prime(latents.haplotypes_);
        }
        model.prime(latents.haplotypes_);
        latents.somatic_model_inferences_ = model.evaluate(latents.cancer_genotypes_, *latents.cancer_genotype_indices_, haplotype_likelihoods);
    } else {
        latents.somatic_model_inferences_ = model.evaluate(latents.cancer_genotypes_, haplotype_likelihoods);
    }
}

auto get_high_posterior_genotypes(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                  const model::SomaticSubcloneModel::InferredLatents& latents)
{
    return copy_greatest_probability_values(genotypes, latents.posteriors.genotype_probabilities, 10, 1e-3);
}

void CancerCaller::evaluate_noise_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    if (has_normal_sample() && !has_high_normal_contamination_risk(latents)) {
        if (!latents.normal_germline_inferences_) {
            assert(latents.germline_model_);
            haplotype_likelihoods.prime(normal_sample());
            if (latents.germline_genotype_indices_) {
                latents.normal_germline_inferences_ = latents.germline_model_->evaluate(latents.germline_genotypes_,
                                                                                        *latents.germline_genotype_indices_,
                                                                                        haplotype_likelihoods);
            } else {
                latents.normal_germline_inferences_ = latents.germline_model_->evaluate(latents.germline_genotypes_,
                                                                                        haplotype_likelihoods);
            }
        }
        assert(latents.cancer_genotype_prior_model_);
        auto noise_model_priors = get_noise_model_priors(*latents.cancer_genotype_prior_model_, latents.somatic_ploidy_);
        const SomaticModel noise_model {{*parameters_.normal_sample}, noise_model_priors};
        auto noise_genotypes = get_high_posterior_genotypes(latents.cancer_genotypes_, latents.somatic_model_inferences_);
        latents.noise_model_inferences_ = noise_model.evaluate(noise_genotypes, haplotype_likelihoods);
    }
}

void CancerCaller::set_model_priors(Latents& latents) const
{
    if (has_normal_sample()) {
        latents.model_priors_ = {.09, 0.01, 0.9};
    } else {
        latents.model_priors_ = {.09, 0.001, 0.909};
    }
}

void CancerCaller::set_model_posteriors(Latents& latents) const
{
    const auto& germline_inferences = latents.germline_model_inferences_;
    const auto& cnv_inferences      = latents.cnv_model_inferences_;
    const auto& somatic_inferences  = latents.somatic_model_inferences_;
    const auto& model_priors        = latents.model_priors_;
    if (debug_log_) {
        stream(*debug_log_) << "Germline model evidence: " << germline_inferences.log_evidence;
        stream(*debug_log_) << "CNV model evidence:      " << cnv_inferences.approx_log_evidence;
        stream(*debug_log_) << "Somatic model evidence:  " << somatic_inferences.approx_log_evidence;
    }
    const auto germline_model_jlp = std::log(model_priors.germline) + germline_inferences.log_evidence;
    const auto cnv_model_jlp      = std::log(model_priors.cnv) + cnv_inferences.approx_log_evidence;
    const auto somatic_model_jlp  = std::log(model_priors.somatic) + somatic_inferences.approx_log_evidence;
    const auto norm = maths::log_sum_exp(germline_model_jlp, cnv_model_jlp, somatic_model_jlp);
    latents.model_posteriors_.germline = std::exp(germline_model_jlp - norm);
    latents.model_posteriors_.cnv      = std::exp(cnv_model_jlp - norm);
    latents.model_posteriors_.somatic  = std::exp(somatic_model_jlp - norm);
    const auto check_sum = latents.model_posteriors_.germline + latents.model_posteriors_.cnv + latents.model_posteriors_.somatic;
    if (check_sum > 1.0) {
        latents.model_posteriors_.germline /= check_sum;
        latents.model_posteriors_.cnv /= check_sum;
        latents.model_posteriors_.somatic /= check_sum;
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
            Priors::GenotypeMixturesDirichletAlphas sample_alphas(parameters_.ploidy, parameters_.concentrations.cnv.normal);
            cnv_alphas.emplace(sample, std::move(sample_alphas));
        } else {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas(parameters_.ploidy, parameters_.concentrations.cnv.tumour);
            cnv_alphas.emplace(sample, std::move(sample_alphas));
        }
    }
    return Priors {prior_model, std::move(cnv_alphas)};
}

auto make_dirichlet_alphas(unsigned n_germline, double germline, unsigned n_somatic, double somatic)
{
    model::SomaticSubcloneModel::Priors::GenotypeMixturesDirichletAlphas result(n_germline + n_somatic);
    std::fill_n(std::begin(result), n_germline, germline);
    std::fill_n(std::rbegin(result), n_somatic, somatic);
    return result;
}

CancerCaller::SomaticModel::Priors
CancerCaller::get_somatic_model_priors(const CancerGenotypePriorModel& prior_model, const unsigned somatic_ploidy) const
{
    using Priors = SomaticModel::Priors;
    Priors::GenotypeMixturesDirichletAlphaMap alphas {};
    alphas.reserve(samples_.size());
    for (const auto& sample : samples_) {
        if (has_normal_sample() && sample == normal_sample()) {
            alphas.emplace(sample, make_dirichlet_alphas(parameters_.ploidy, parameters_.concentrations.somatic.normal_germline,
                                                         somatic_ploidy, parameters_.concentrations.somatic.normal_somatic));
        } else {
            alphas.emplace(sample, make_dirichlet_alphas(parameters_.ploidy, parameters_.concentrations.somatic.tumour_germline,
                                                         somatic_ploidy, parameters_.concentrations.somatic.tumour_somatic));
        }
    }
    return Priors {prior_model, std::move(alphas)};
}

CancerCaller::SomaticModel::Priors
CancerCaller::get_noise_model_priors(const CancerGenotypePriorModel& prior_model, const unsigned somatic_ploidy) const
{
    // The noise model is intended to capture noise that may also be present in the normal sample,
    // hence all samples have the same prior alphas.
    using Priors = SomaticModel::Priors;
    auto noise_alphas = make_dirichlet_alphas(parameters_.ploidy, parameters_.concentrations.somatic.normal_germline,
                                              somatic_ploidy, parameters_.concentrations.somatic.tumour_somatic);
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

auto compute_marginal_credible_interval(const model::SomaticSubcloneModel::Priors::GenotypeMixturesDirichletAlphas& alphas,
                                        const std::size_t k, const double mass)
{
    const auto a0 = std::accumulate(std::cbegin(alphas), std::cend(alphas), 0.0);
    return maths::beta_hdi(alphas[k], a0 - alphas[k], mass);
}

auto compute_marginal_credible_intervals(const model::SomaticSubcloneModel::Priors::GenotypeMixturesDirichletAlphas& alphas,
                                         const double mass)
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

auto compute_marginal_credible_intervals(const model::SomaticSubcloneModel::Priors::GenotypeMixturesDirichletAlphaMap& alphas,
                                         const double mass)
{
    CredibleRegionMap result {};
    result.reserve(alphas.size());
    for (const auto& p : alphas) {
        result.emplace(p.first, compute_marginal_credible_intervals(p.second, mass));
    }
    return result;
}

auto compute_credible_somatic_mass(const model::SomaticSubcloneModel::Priors::GenotypeMixturesDirichletAlphas& alphas,
                                   const unsigned somatic_ploidy, const double min_credible_somatic_frequency)
{
    if (somatic_ploidy == 1) {
        return maths::dirichlet_marginal_sf(alphas, alphas.size() - 1, min_credible_somatic_frequency);
    } else {
        double inv_result {1.0};
        for (unsigned i {1}; i <= somatic_ploidy; ++i) {
            inv_result *= maths::dirichlet_marginal_cdf(alphas, alphas.size() - i, min_credible_somatic_frequency);
        }
        return 1.0 - inv_result;
    }
}

auto compute_credible_somatic_mass(const model::SomaticSubcloneModel::Priors::GenotypeMixturesDirichletAlphaMap& alphas,
                                   const unsigned somatic_ploidy, const double min_credible_somatic_frequency)
{
    double inv_result {1.0};
    for (const auto& p : alphas) {
        if (somatic_ploidy == 1) {
            inv_result *= maths::dirichlet_marginal_cdf(p.second, p.second.size() - 1, min_credible_somatic_frequency);
        } else {
            inv_result *= 1.0 - compute_credible_somatic_mass(p.second, somatic_ploidy, min_credible_somatic_frequency);
        }
    }
    return 1.0 - inv_result;
}

auto compute_map_somatic_vaf(const model::SomaticSubcloneModel::Priors::GenotypeMixturesDirichletAlphas& alphas,
                             const unsigned somatic_ploidy)
{
    double result {0.0};
    for (unsigned i {1}; i <= somatic_ploidy; ++i) {
        result = std::max(maths::dirichlet_expectation(alphas.size() - i, alphas), result);
    }
    return result;
}

using SomaticVAFMap = std::unordered_map<SampleName, double>;

auto compute_map_somatic_vafs(const model::SomaticSubcloneModel::Priors::GenotypeMixturesDirichletAlphaMap& alphas,
                              const unsigned somatic_ploidy)
{
    SomaticVAFMap result {};
    result.reserve(alphas.size());
    for (const auto& p : alphas) {
        result.emplace(p.first, compute_map_somatic_vaf(p.second, somatic_ploidy));
    }
    return result;
}

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
    Phred<double> posterior, segregation_quality;
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
    Phred<double> posterior, segregation_quality;
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
    
    Genotype<Allele> genotype, somatic;
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
    CredibleRegionMap credible_regions;
    SomaticVAFMap somatic_map_vafs;
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
    auto p = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                             0.0, [&allele] (const auto curr, const auto& p) {
                                 return curr + (contains(p.first, allele) ? 0.0 : p.second);
                             });
    return probability_false_to_phred(p);
}

template <typename M>
VariantPosteriorVector compute_candidate_posteriors(const std::vector<Variant>& candidates, const M& genotype_posteriors)
{
    VariantPosteriorVector result {};
    result.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        result.emplace_back(candidate, marginalise(candidate.alt_allele(), genotype_posteriors));
    }
    return result;
}

// segregation probability

using BigFloat = boost::multiprecision::mpf_float_1000;

BigFloat marginalise(const Allele& allele, const std::vector<Genotype<Haplotype>>& genotypes, const std::vector<double>& probabilities)
{
    assert(genotypes.size() == probabilities.size());
    auto inv_result = std::inner_product(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(probabilities),
                                         0.0, std::plus<> {}, [&allele] (const auto& genotype, const auto probability) {
        return contains(genotype, allele) ? 0.0 : probability; });
    return BigFloat {1.0} - BigFloat {inv_result};
}

bool is_somatic(const Allele& allele, const CancerGenotype<Haplotype>& genotype)
{
    return contains(genotype.somatic(), allele) && !contains(genotype.germline(), allele);
}

BigFloat marginalise(const Allele& allele, const std::vector<CancerGenotype<Haplotype>>& genotypes, const std::vector<double>& probabilities,
                     const BigFloat somatic_mass_complement)
{
    assert(genotypes.size() == probabilities.size());
    const BigFloat contained_complement {std::inner_product(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(probabilities),
                                                  0.0, std::plus<> {}, [&] (const auto& genotype, const auto probability) {
        return contains(genotype, allele) ? 0.0 : probability; })};
    BigFloat somatic_complement {std::inner_product(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(probabilities),
                                                 0.0, std::plus<> {}, [&] (const auto& genotype, const auto probability) {
        return is_somatic(allele, genotype) ? probability : 0.0; })};
    somatic_complement *= somatic_mass_complement;
    const BigFloat result_complement {contained_complement + somatic_complement};
    assert(result_complement >= 0.0 && result_complement <= 1.0);
    return BigFloat {1.0} - result_complement;
}

Phred<double>
calculate_segregation_probability(const Allele& allele,
                                  const std::vector<Genotype<Haplotype>>& germline_genotypes,
                                  const std::vector<CancerGenotype<Haplotype>>& cancer_genotypes,
                                  const std::vector<double>& germline_genotype_probabilities,
                                  const std::vector<double>& cnv_genotype_probabilities,
                                  const std::vector<double>& cancer_genotype_probabilities,
                                  const BigFloat germline_probability,
                                  const BigFloat cnv_probability,
                                  const BigFloat somatic_probability,
                                  const BigFloat somatic_mass)
{
    auto prob_germline_segregates = marginalise(allele, germline_genotypes, germline_genotype_probabilities);
    prob_germline_segregates *= germline_probability;
    auto prob_cnv_segregates = marginalise(allele, germline_genotypes, cnv_genotype_probabilities);
    prob_cnv_segregates *= cnv_probability;
    auto prob_somatic_segregates = marginalise(allele, cancer_genotypes, cancer_genotype_probabilities, BigFloat {1.0} - somatic_mass);
    prob_somatic_segregates *= somatic_probability;
    const BigFloat prob_segregates {prob_germline_segregates + prob_cnv_segregates + prob_somatic_segregates};
    return probability_true_to_phred<double>(prob_segregates);
}

Phred<double>
calculate_segregation_probability(const Allele& allele,
                                  const std::vector<Genotype<Haplotype>>& germline_genotypes,
                                  const std::vector<CancerGenotype<Haplotype>>& cancer_genotypes,
                                  const std::vector<double>& germline_genotype_probabilities,
                                  const std::vector<double>& cnv_genotype_probabilities,
                                  const std::vector<double>& cancer_genotype_probabilities,
                                  const double germline_probability,
                                  const double cnv_probability,
                                  const double somatic_probability,
                                  const double somatic_mass)
{
    BigFloat germline_bf {germline_probability}, cnv_bf {cnv_probability}, somatic_bf {somatic_probability};
    const BigFloat norm {germline_bf + cnv_bf + somatic_bf};
    germline_bf /= norm; cnv_bf /= norm; somatic_bf /= norm;
    return calculate_segregation_probability(allele, germline_genotypes, cancer_genotypes,
                                             germline_genotype_probabilities, cnv_genotype_probabilities, cancer_genotype_probabilities,
                                             germline_bf, cnv_bf, somatic_bf, BigFloat {somatic_mass});
}

Phred<double> calculate_somatic_posterior(const double somatic_model_posterior, const double somatic_mass)
{
    BigFloat somatic_posterior {somatic_model_posterior};
    somatic_posterior *= somatic_mass;
    return probability_true_to_phred<double>(somatic_posterior);
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
        if (p.second >= min_posterior && contains_alt(genotype_call, p.first)) {
            calls.emplace_back(p.first, p.second);
        } else {
            uncalled.emplace_back(p.first);
        }
    }
    return std::make_pair(std::move(calls), std::move(uncalled));
}

// somatic variant posterior

BigFloat marginalise_somatic(const Allele& allele, const std::vector<CancerGenotype<Haplotype>>& genotypes,
                             const std::vector<double>& probabilities)
{
    BigFloat result_complement {std::inner_product(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(probabilities),
                                                   0.0, std::plus<> {}, [&allele] (const auto& genotype, auto probability) {
        return is_somatic(allele, genotype) ? 0.0 : probability; })};
    return BigFloat {1.0} - result_complement;
}

auto compute_somatic_variant_posteriors(const std::vector<VariantReference>& candidates,
                                        const std::vector<CancerGenotype<Haplotype>>& cancer_genotypes,
                                        const std::vector<double>& cancer_genotype_posteriors,
                                        const BigFloat somatic_posterior)
{
    VariantPosteriorVector result {};
    result.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        auto p = marginalise_somatic(candidate.get().alt_allele(), cancer_genotypes, cancer_genotype_posteriors);
        p *= somatic_posterior;
        result.emplace_back(candidate, probability_true_to_phred<double>(p));
    }
    return result;
}

auto compute_somatic_variant_posteriors(const std::vector<VariantReference>& candidates,
                                        const std::vector<CancerGenotype<Haplotype>>& cancer_genotypes,
                                        const std::vector<double>& cancer_genotype_posteriors,
                                        const Phred<double> somatic_posterior)
{
    return compute_somatic_variant_posteriors(candidates, cancer_genotypes, cancer_genotype_posteriors,
                                              BigFloat {somatic_posterior.probability_true().value});
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

Phred<double>
marginalise(const CancerGenotype<Allele>& genotype, const std::vector<CancerGenotype<Haplotype>>& genotypes,
            const std::vector<double>& genotype_posteriors)
{
    auto p = std::inner_product(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(genotype_posteriors),
                                0.0, std::plus<> {},  [&genotype] (const auto& g, auto probability) {
                                    return contains(g, genotype) ? 0.0 : probability; });
    return probability_false_to_phred(p);
}

auto call_somatic_genotypes(const CancerGenotype<Haplotype>& called_genotype,
                            const std::vector<GenomicRegion>& called_somatic_regions,
                            const std::vector<CancerGenotype<Haplotype>>& genotypes,
                            const std::vector<double>& genotype_posteriors,
                            const CredibleRegionMap& credible_regions,
                            const SomaticVAFMap& somatic_vafs)
{
    CancerGenotypeCalls result {};
    result.reserve(called_somatic_regions.size());
    for (const auto& region : called_somatic_regions) {
        auto genotype_chunk = copy<Allele>(called_genotype, region);
        auto posterior = marginalise(genotype_chunk, genotypes, genotype_posteriors);
        result.emplace_back(std::move(genotype_chunk), posterior);
        result.back().credible_regions = credible_regions;
        result.back().somatic_map_vafs = somatic_vafs;
    }
    return result;
}

// output

octopus::VariantCall::GenotypeCall demote(GermlineGenotypeCall call)
{
    return octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
}

std::unique_ptr<octopus::VariantCall>
transform_germline_call(GermlineVariantCall&& variant_call, GermlineGenotypeCall&& genotype_call,
                        const std::vector<SampleName>& samples, const std::vector<SampleName>& somatic_samples)
{
    std::vector<std::pair<SampleName, Call::GenotypeCall>> genotypes {};
    for (const auto& sample : samples) {
        if (std::find(std::cbegin(somatic_samples), std::cend(somatic_samples), sample) == std::cend(somatic_samples)) {
            genotypes.emplace_back(sample, demote(genotype_call));
        } else {
            auto copy = genotype_call;
            for (auto allele : genotype_call.somatic) copy.genotype.emplace(allele);
            genotypes.emplace_back(sample, demote(std::move(copy)));
        }
    }
    return std::make_unique<octopus::GermlineVariantCall>(variant_call.variant.get(), std::move(genotypes),
                                                          variant_call.segregation_quality, variant_call.posterior);
}

template <typename Container, typename T>
auto find_index(const Container& values, const T& value)
{
    const auto itr = std::find(std::cbegin(values), std::cend(values), value);
    return itr != std::cend(values) ? std::distance(std::cbegin(values), itr) : -1;
}

auto transform_somatic_calls(SomaticVariantCalls&& somatic_calls, CancerGenotypeCalls&& genotype_calls,
                             const std::vector<SampleName>& somatic_samples)
{
    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
    result.reserve(somatic_calls.size());
    std::transform(std::make_move_iterator(std::begin(somatic_calls)), std::make_move_iterator(std::end(somatic_calls)),
                   std::make_move_iterator(std::begin(genotype_calls)), std::back_inserter(result),
                   [&somatic_samples] (auto&& variant_call, auto&& genotype_call) -> std::unique_ptr<octopus::VariantCall> {
                       std::unordered_map<SampleName, SomaticCall::GenotypeCredibleRegions> credible_regions {};
                       const auto germline_ploidy = genotype_call.genotype.germline_ploidy();
                       for (const auto& p : genotype_call.credible_regions) {
                           SomaticCall::GenotypeCredibleRegions sample_credible_regions {};
                           sample_credible_regions.germline.reserve(germline_ploidy);
                           std::copy(std::cbegin(p.second), std::prev(std::cend(p.second)),
                                     std::back_inserter(sample_credible_regions.germline));
                           if (std::find(std::cbegin(somatic_samples), std::cend(somatic_samples), p.first) != std::cend(somatic_samples)) {
                               auto somatic_idx = find_index(genotype_call.genotype.somatic(), variant_call.variant.get().alt_allele());
                               sample_credible_regions.somatic = p.second[germline_ploidy + somatic_idx];
                           }
                           credible_regions.emplace(p.first, std::move(sample_credible_regions));
                       }
                       return std::make_unique<SomaticCall>(variant_call.variant.get(), std::move(genotype_call.genotype),
                                                            genotype_call.posterior, std::move(credible_regions),
                                                            genotype_call.somatic_map_vafs,
                                                            variant_call.segregation_quality, variant_call.posterior);
                   });
    return result;
}

} // namespace

Phred<double>
CancerCaller::calculate_segregation_probability(const Variant& variant, const Latents& latents, double somatic_mass) const
{
    return octopus::calculate_segregation_probability(variant.alt_allele(), latents.germline_genotypes_, latents.cancer_genotypes_,
                                                      latents.germline_model_inferences_.posteriors.genotype_probabilities,
                                                      latents.cnv_model_inferences_.posteriors.genotype_probabilities,
                                                      latents.somatic_model_inferences_.posteriors.genotype_probabilities,
                                                      latents.model_posteriors_.germline, latents.model_posteriors_.cnv,
                                                      latents.model_posteriors_.somatic, somatic_mass);
}

namespace debug {

template <typename S, typename T>
void print_variants(S&& stream, const std::vector<T>& variants)
{
    for (const auto& v : variants) stream << v.variant << " " << v.posterior << '\n';
}

} // namespace debug

std::vector<std::unique_ptr<VariantCall>>
CancerCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    // TODO: refactor this into smaller methods!
    const auto conditional_somatic_mass = calculate_somatic_mass(latents);
    const auto& model_posteriors = latents.model_posteriors_;
    log(model_posteriors);
    const auto somatic_posterior = calculate_somatic_posterior(latents.model_posteriors_.somatic, conditional_somatic_mass);
    const auto germline_genotype_posteriors = calculate_germline_genotype_posteriors(latents);
    const auto& cancer_genotype_posteriors = latents.somatic_model_inferences_.posteriors.genotype_probabilities;
    log(latents.germline_genotypes_, germline_genotype_posteriors, latents.germline_model_inferences_, latents.cnv_model_inferences_,
        latents.cancer_genotypes_, latents.somatic_model_inferences_);
    const auto germline_candidate_posteriors = compute_candidate_posteriors(candidates, germline_genotype_posteriors);
    boost::optional<Genotype<Haplotype>> called_germline_genotype {};
    boost::optional<CancerGenotype<Haplotype>> called_cancer_genotype {};
    if (model_posteriors.somatic > model_posteriors.germline && somatic_posterior >= parameters_.min_somatic_posterior) {
        if (debug_log_) *debug_log_ << "Using cancer genotype for germline genotype call";
        if (!called_cancer_genotype) {
            auto cancer_posteriors = zip_cref(latents.cancer_genotypes_, cancer_genotype_posteriors);
            called_cancer_genotype = find_map_genotype(cancer_posteriors)->first;
        }
        called_germline_genotype = called_cancer_genotype->germline();
    } else {
        called_germline_genotype = find_map_genotype(germline_genotype_posteriors)->first;
    }
    GermlineVariantCalls germline_variant_calls;
    std::vector<VariantReference> uncalled_germline_candidates;
    std::tie(germline_variant_calls, uncalled_germline_candidates) = call_candidates(germline_candidate_posteriors,
                                                                                     *called_germline_genotype,
                                                                                     parameters_.min_variant_posterior);
    
    for (auto& v : germline_variant_calls) {
        v.segregation_quality = calculate_segregation_probability(v.variant, latents, conditional_somatic_mass);
        if (v.posterior > v.segregation_quality) v.posterior = v.segregation_quality;
    }
    
    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
    Genotype<Haplotype> called_somatic_genotype {};
    std::vector<SampleName> somatic_samples {};
    if (somatic_posterior >= parameters_.min_somatic_posterior) {
        auto somatic_allele_posteriors = compute_somatic_variant_posteriors(uncalled_germline_candidates, latents.cancer_genotypes_,
                                                                            cancer_genotype_posteriors, somatic_posterior);
        if (!called_cancer_genotype) {
            auto cancer_posteriors = zip_cref(latents.cancer_genotypes_, cancer_genotype_posteriors);
            called_cancer_genotype = find_map_genotype(cancer_posteriors)->first.get();
        }
        if (called_cancer_genotype->germline() == called_germline_genotype) {
            auto somatic_variant_calls = call_somatic_variants(somatic_allele_posteriors, *called_cancer_genotype,
                                                               parameters_.min_somatic_posterior);
            const auto& somatic_alphas = latents.somatic_model_inferences_.posteriors.alphas;
            const auto credible_regions = compute_marginal_credible_intervals(somatic_alphas, parameters_.credible_mass);
            if (!somatic_variant_calls.empty()) {
                for (const auto& p : credible_regions) {
                    if (debug_log_) {
                        auto ss = stream(*debug_log_);
                        ss << p.first << " somatic credible regions: ";
                        for (auto cr : p.second) ss << '(' << cr.first << ' ' << cr.second << ") ";
                    }
                    if (std::any_of(std::next(std::cbegin(p.second), parameters_.ploidy), std::cend(p.second),
                        [this] (const auto& credible_region) { return credible_region.first >= parameters_.min_credible_somatic_frequency; })) {
                        if (has_normal_sample() && p.first == normal_sample()) {
                            somatic_samples.clear();
                            break;
                        }
                        somatic_samples.push_back(p.first);
                    }
                }
                if (latents.noise_model_inferences_ && latents.normal_germline_inferences_) {
                    const auto noise_model_evidence = latents.noise_model_inferences_->approx_log_evidence;
                    const auto germline_model_evidence = latents.normal_germline_inferences_->log_evidence;
                    if (noise_model_evidence > germline_model_evidence) {
                        // Does the normal sample contain the called somatic variant?
                        const auto& noisy_alphas = latents.noise_model_inferences_->posteriors.alphas.at(normal_sample());
                        const auto noise_mass = compute_credible_somatic_mass(noisy_alphas, latents.somatic_ploidy_, parameters_.min_expected_somatic_frequency);
                        if (noise_mass > 2 * parameters_.min_credible_somatic_frequency) {
                            somatic_samples.clear();
                        }
                    }
                }
                if (somatic_samples.empty()) {
                    somatic_variant_calls.clear();
                    somatic_variant_calls.shrink_to_fit();
                } else {
                    called_somatic_genotype = called_cancer_genotype->somatic();
                }
                for (auto& v : somatic_variant_calls) {
                    v.segregation_quality = calculate_segregation_probability(v.variant, latents, conditional_somatic_mass);
                    if (v.posterior > v.segregation_quality) v.posterior = v.segregation_quality;
                }
            }
            if (debug_log_) {
                *debug_log_ << "Called somatic variants:";
                debug::print_variants(stream(*debug_log_), somatic_variant_calls);
            }
            const auto somatic_vafs = compute_map_somatic_vafs(somatic_alphas, latents.somatic_ploidy_);
            const auto called_somatic_regions = extract_regions(somatic_variant_calls);
            auto cancer_genotype_calls = call_somatic_genotypes(*called_cancer_genotype, called_somatic_regions,
                                                                latents.cancer_genotypes_, cancer_genotype_posteriors,
                                                                credible_regions, somatic_vafs);
            result = transform_somatic_calls(std::move(somatic_variant_calls), std::move(cancer_genotype_calls), somatic_samples);
        } else if (debug_log_) {
            stream(*debug_log_) << "Conflict between called germline genotype and called cancer genotype. Not calling somatics";
        }
    }
    const auto called_germline_regions = extract_regions(germline_variant_calls);
    GermlineGenotypeCalls germline_genotype_calls {};
    germline_genotype_calls.reserve(called_germline_regions.size());
    for (const auto& region : called_germline_regions) {
        auto genotype_chunk = copy<Allele>(*called_germline_genotype, region);
        const auto inv_posterior = std::accumulate(std::cbegin(germline_genotype_posteriors),
                                                   std::cend(germline_genotype_posteriors), 0.0,
                                                   [&called_germline_genotype] (const double curr, const auto& p) {
                                                       return curr + (contains(p.first, *called_germline_genotype) ? 0.0 : p.second);
                                                   });
        if (called_somatic_genotype.ploidy() > 0) {
            germline_genotype_calls.emplace_back(std::move(genotype_chunk),
                                                 copy<Allele>(called_somatic_genotype, region),
                                                 probability_false_to_phred(inv_posterior));
        } else {
            germline_genotype_calls.emplace_back(std::move(genotype_chunk), probability_false_to_phred(inv_posterior));
        }
    }
    if (debug_log_) {
        *debug_log_ << "Called germline variants:";
        debug::print_variants(stream(*debug_log_), germline_variant_calls);
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

CancerCaller::GermlineGenotypeProbabilityMap
CancerCaller::calculate_germline_genotype_posteriors(const Latents& latents) const
{
    const auto& model_posteriors = latents.model_posteriors_;
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
    const auto& somatic_posteriors = latents.somatic_model_inferences_.posteriors.genotype_probabilities;
    for (std::size_t i {0}; i < cancer_genotypes.size(); ++i) {
        result[cancer_genotypes[i].germline()] += model_posteriors.somatic * somatic_posteriors[i];
    }
    return result;
}

double CancerCaller::calculate_somatic_mass(const CancerCaller::Latents& latents) const
{
    return compute_credible_somatic_mass(latents.somatic_model_inferences_.posteriors.alphas, latents.somatic_ploidy_,
                                         parameters_.min_expected_somatic_frequency);
}

std::vector<std::unique_ptr<ReferenceCall>>
CancerCaller::call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents, const ReadPileupMap& pileups) const
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

CancerCaller::Latents::Latents(const std::vector<Haplotype>& haplotypes, const std::vector<SampleName>& samples,
                               const CancerCaller::Parameters& parameters)
: haplotypes_ {haplotypes}
, samples_ {samples}
, parameters_ {parameters}
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

void CancerCaller::Latents::compute_haplotype_posteriors() const
{
    Latents::HaplotypeProbabilityMap result {haplotypes_.get().size()};
    for (const auto& haplotype : haplotypes_.get()) {
        result.emplace(haplotype, 0.0);
    }
    // Contribution from germline model
    for (const auto& p : zip(germline_genotypes_, germline_model_inferences_.posteriors.genotype_probabilities)) {
        for (const auto& haplotype : p.get<0>().copy_unique_ref()) {
            result.at(haplotype) += model_posteriors_.germline * p.get<1>();
        }
    }
    // Contribution from CNV model
    for (const auto& p : zip(germline_genotypes_, cnv_model_inferences_.posteriors.genotype_probabilities)) {
        for (const auto& haplotype : p.get<0>().copy_unique_ref()) {
            result.at(haplotype) += model_posteriors_.cnv * p.get<1>();
        }
    }
    if (!cancer_genotypes_.empty()) {
        const auto credible_frequency = parameters_.get().min_expected_somatic_frequency;
        const auto conditional_somatic_prob = compute_credible_somatic_mass(somatic_model_inferences_.posteriors.alphas, somatic_ploidy_, credible_frequency);
        // Contribution from somatic model
        for (const auto& p : zip(cancer_genotypes_, somatic_model_inferences_.posteriors.genotype_probabilities)) {
            for (const auto& haplotype : p.get<0>().germline().copy_unique_ref()) {
                result.at(haplotype) += model_posteriors_.somatic * p.get<1>();
            }
            for (const auto& haplotype : p.get<0>().somatic().copy_unique_ref()) {
                result.at(haplotype) += model_posteriors_.somatic * conditional_somatic_prob * p.get<1>();
            }
        }
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
                       const GermlineModel::InferredLatents& germline_inferences,
                       const CNVModel::InferredLatents& cnv_inferences,
                       const CancerGenotypeVector& cancer_genotypes,
                       const SomaticModel::InferredLatents& somatic_inferences) const
{
    if (debug_log_) {
        auto germline_posteriors = zip_cref(germline_genotypes, germline_inferences.posteriors.genotype_probabilities);
        auto map_germline = find_map_genotype(germline_posteriors);
        auto germline_log = stream(*debug_log_);
        germline_log << "MAP germline genotype: ";
        debug::print_variant_alleles(germline_log, map_germline->first);
        germline_log << ' ' << map_germline->second;
        auto cnv_posteriors = zip_cref(germline_genotypes, cnv_inferences.posteriors.genotype_probabilities);
        auto map_cnv = find_map_genotype(cnv_posteriors);
        auto cnv_log = stream(*debug_log_);
        cnv_log << "MAP CNV genotype: ";
        debug::print_variant_alleles(cnv_log, map_cnv->first);
        cnv_log << ' ' << map_cnv->second;
        auto somatic_log = stream(*debug_log_);
        auto cancer_posteriors = zip_cref(cancer_genotypes, somatic_inferences.posteriors.genotype_probabilities);
        auto map_somatic = find_map_genotype(cancer_posteriors);
        auto map_cancer_genotype = map_somatic->first.get();
        somatic_log << "MAP cancer genotype: ";
        debug::print_variant_alleles(somatic_log, map_cancer_genotype);
        somatic_log << ' ' << map_somatic->second;
        auto map_marginal_germline = find_map_genotype(germline_genotype_posteriors);
        auto marginal_germline_log = stream(*debug_log_);
        marginal_germline_log << "MAP marginal germline genotype: ";
        debug::print_variant_alleles(marginal_germline_log, map_marginal_germline->first);
        marginal_germline_log << ' ' << map_marginal_germline->second;
    }
}

} // namespace octopus
