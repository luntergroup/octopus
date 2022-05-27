// Copyright (c) 2015-2021 Daniel Cooke
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
#include <limits>

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
#include "utils/map_utils.hpp"
#include "logging/logging.hpp"
#include "core/types/calls/germline_variant_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "core/types/calls/somatic_call.hpp"
#include "core/types/calls/cnv_call.hpp"

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
    if (parameters_.max_genotypes && *parameters_.max_genotypes == 0) {
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
        std::type_index(typeid(SomaticCall)),
        std::type_index(typeid(CNVCall))
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

std::size_t CancerCaller::do_remove_duplicates(HaplotypeBlock& haplotypes) const
{
    if (parameters_.deduplicate_haplotypes_with_germline_model) {
        if (haplotypes.size() < 2) return 0;
        CoalescentModel::Parameters model_params {};
        if (parameters_.germline_prior_model_params) model_params = *parameters_.germline_prior_model_params;
        Haplotype reference {mapped_region(haplotypes), reference_.get()};
        CoalescentModel model {std::move(reference), model_params, haplotypes.size(), CoalescentModel::CachingStrategy::none};
        const CoalescentProbabilityGreater cmp {std::move(model)};
        return octopus::remove_duplicates(haplotypes, cmp);
    } else {
        return Caller::do_remove_duplicates(haplotypes);
    }
}

std::unique_ptr<CancerCaller::Caller::Latents>
CancerCaller::infer_latents(const HaplotypeBlock& haplotypes,
                            const HaplotypeLikelihoodArray& haplotype_likelihoods,
                            OptionalThreadPool workers) const
{
    // Store any intermediate results in Latents for reuse, so the order of model evaluation matters!
    auto result = std::make_unique<Latents>(haplotypes, samples_, parameters_);
    set_model_priors(*result);
    generate_germline_genotypes(*result, result->indexed_haplotypes_);
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

boost::optional<Caller::ModelPosterior>
CancerCaller::calculate_model_posterior(const HaplotypeBlock& haplotypes,
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
    latents.max_evidence_somatic_model_index_ = 0;
    latents.cancer_genotypes_.reserve(parameters_.max_somatic_haplotypes);
    latents.somatic_model_inferences_.reserve(parameters_.max_somatic_haplotypes);
    latents.somatic_model_posteriors_.reserve(latents.somatic_model_inferences_.size());
    for (unsigned somatic_ploidy {1}; somatic_ploidy <= parameters_.max_somatic_haplotypes; ++somatic_ploidy) {
        if (debug_log_) stream(*debug_log_) << "Fitting somatic model with somatic ploidy " << somatic_ploidy;
        latents.inferred_somatic_ploidy_ = somatic_ploidy;
        generate_cancer_genotypes(latents, haplotype_likelihoods);
        if (debug_log_) stream(*debug_log_) << "There are " << latents.cancer_genotypes_.back().size() << " candidate cancer genotypes";
        evaluate_somatic_model(latents, haplotype_likelihoods);
        latents.somatic_model_posteriors_.push_back(latents.somatic_model_inferences_.back().approx_log_evidence);
        if (debug_log_) stream(*debug_log_) << "Evidence for somatic model with somatic ploidy "
                             << somatic_ploidy << " is " << latents.somatic_model_posteriors_.back();
        
        if (somatic_ploidy > 1) {
            if (latents.somatic_model_inferences_.back().approx_log_evidence
              < latents.somatic_model_inferences_[somatic_ploidy - 2].approx_log_evidence) {
                  latents.inferred_somatic_ploidy_ = somatic_ploidy - 1;
                break;
            }
        } else {
            set_model_posteriors(latents);
            if (latents.model_posteriors_.somatic < std::max(latents.model_posteriors_.germline, latents.model_posteriors_.cnv)) {
                break;
            }
        }
        if (latents.haplotypes_.get().size() <= somatic_ploidy + 1) break;
    }
    maths::normalise_exp(latents.somatic_model_posteriors_);
    if (debug_log_) stream(*debug_log_) << "Best somatic model has somatic ploidy " << latents.inferred_somatic_ploidy_;
    latents.max_evidence_somatic_model_index_ = latents.inferred_somatic_ploidy_ - 1;
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

auto demote_each(const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& genotypes)
{
    MappableBlock<Genotype<IndexedHaplotype<>>> result {};
    result.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::back_inserter(result),
                   [] (const auto& genotype) { return demote(genotype); });
    return result;
}

} // namespace

boost::optional<Caller::ModelPosterior>
CancerCaller::calculate_model_posterior(const HaplotypeBlock& haplotypes,
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
        const auto dummy_genotypes = demote_each(latents.cancer_genotypes_[latents.max_evidence_somatic_model_index_]);
        const auto dummy_inferences = germline_model.evaluate(dummy_genotypes, haplotype_likelihoods);
        ModelPosterior result {};
        if (latents.noise_model_inferences_) {
            result.joint = octopus::calculate_model_posterior(normal_inferences.log_evidence,
                                                              dummy_inferences.log_evidence,
                                                              latents.noise_model_inferences_->approx_log_evidence);
        } else {
            result.joint = octopus::calculate_model_posterior(normal_inferences.log_evidence,
                                                              dummy_inferences.log_evidence);
        }
        return result;
    } else {
        // TODO
        return boost::none;
    }
}

void CancerCaller::generate_germline_genotypes(Latents& latents, const IndexedHaplotypeBlock& haplotypes) const
{
    latents.germline_genotypes_ = generate_all_genotypes(haplotypes, parameters_.ploidy);
}

namespace {

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
    auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end   = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}

template <typename Range>
auto zip_cref(const Range& values, const std::vector<double>& probabilities)
{
    using ValueReference = std::reference_wrapper<const typename Range::value_type>;
    std::vector<std::pair<ValueReference, double>> result {};
    result.reserve(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::cbegin(probabilities), std::back_inserter(result),
                   [] (const auto& v, const auto& p) noexcept { return std::make_pair(std::cref(v), p); });
    return result;
}

template <typename G>
MappableBlock<G>
copy_greatest_probability_values(const MappableBlock<G>& values,
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
    MappableBlock<G> result {mapped_region(values)};
    result.reserve(std::distance(std::begin(value_probabilities), last_include_itr));
    std::transform(std::begin(value_probabilities), last_include_itr, std::back_inserter(result),
                   [] (const auto& p) { return p.first.get(); });
    return result;
}

template <typename G>
auto copy_greatest_probability_genotypes(const MappableBlock<G>& genotypes,
                                         const std::vector<double>& probabilities,
                                         const std::size_t n,
                                         const boost::optional<double> min_include_probability = boost::none,
                                         const boost::optional<double> max_exclude_probability = boost::none)
{
    assert(genotypes.size() == probabilities.size());
    return copy_greatest_probability_values(genotypes, probabilities, n, min_include_probability, max_exclude_probability);
}

auto calculate_posteriors_with_germline_likelihood_model(const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& genotypes,
                                                         const CancerGenotypePriorModel& prior_model,
                                                         const model::ConstantMixtureGenotypeLikelihoodModel likelihood_model,
                                                         const std::vector<SampleName>& samples)
{
    auto result = evaluate(genotypes, prior_model);
    for (const auto& sample : samples) {
        likelihood_model.cache().prime(sample);
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(result), std::begin(result),
                       [&] (const auto& genotype, auto curr) {
                           return curr + likelihood_model.evaluate(demote(genotype));
                       });
    }
    maths::normalise_exp(result);
    return result;
}

void filter_with_germline_model(MappableBlock<CancerGenotype<IndexedHaplotype<>>>& genotypes,
                                const CancerGenotypePriorModel& prior_model,
                                const model::ConstantMixtureGenotypeLikelihoodModel likelihood_model,
                                const std::vector<SampleName>& samples,
                                const std::size_t n)
{
    const auto germline_model_posteriors = calculate_posteriors_with_germline_likelihood_model(genotypes, prior_model, likelihood_model, samples);
    genotypes = copy_greatest_probability_genotypes(genotypes, germline_model_posteriors, n);
}

} // namespace

void CancerCaller::generate_cancer_genotypes(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto num_haplotypes = latents.indexed_haplotypes_.size();
    if (num_haplotypes == 1) return;
    const auto& germline_genotypes = latents.germline_genotypes_;
    const auto num_germline_genotypes = germline_genotypes.size();
    const auto max_possible_cancer_genotypes = num_haplotypes * num_germline_genotypes;
    if (!parameters_.max_genotypes || max_possible_cancer_genotypes <= *parameters_.max_genotypes) {
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

namespace {

struct CancerGenotypeLess
{
    template <typename T>
    bool operator()(const CancerGenotype<T>& lhs, const CancerGenotype<T>& rhs) const
    {
        return lhs.germline() == rhs.germline() ? GenotypeLess()(lhs.somatic(), rhs.somatic()) : GenotypeLess()(lhs.germline(), rhs.germline());
    }
};

template <typename IndexType>
void erase_duplicates(MappableBlock<CancerGenotype<IndexedHaplotype<IndexType>>>& genotypes)
{
    using std::begin; using std::end;
    std::sort(begin(genotypes), end(genotypes), CancerGenotypeLess {});
    genotypes.erase(std::unique(begin(genotypes), end(genotypes)), end(genotypes));
}

} // namespace

void CancerCaller::generate_cancer_genotypes_with_clean_normal(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(parameters_.max_genotypes);
    const auto num_haplotypes = latents.indexed_haplotypes_.size();
    const auto& germline_genotypes = latents.germline_genotypes_;
    const auto max_allowed_cancer_genotypes = *parameters_.max_genotypes;
    if (!latents.cancer_genotypes_.empty()) {
        const auto max_old_cancer_genotype_bases = std::max(max_allowed_cancer_genotypes / num_haplotypes, std::size_t {1});
        const auto& cancer_genotype_posteriors = latents.somatic_model_inferences_.back().weighted_genotype_posteriors;
        const auto old_cancer_genotype_bases = copy_greatest_probability_values(latents.cancer_genotypes_.back(), cancer_genotype_posteriors, max_old_cancer_genotype_bases);
        latents.cancer_genotypes_.push_back(extend_somatic(old_cancer_genotype_bases, latents.indexed_haplotypes_));
        erase_duplicates(latents.cancer_genotypes_.back());
    } else {
        assert(latents.germline_model_);
        haplotype_likelihoods.prime(normal_sample());
        latents.normal_germline_inferences_ = latents.germline_model_->evaluate(germline_genotypes, haplotype_likelihoods);
        const auto& germline_normal_posteriors = latents.normal_germline_inferences_->posteriors.genotype_probabilities;
        const auto max_germline_genotype_bases = calculate_max_germline_genotype_bases(max_allowed_cancer_genotypes, num_haplotypes, 1);
        MappableBlock<Genotype<IndexedHaplotype<>>> germline_bases;
        germline_bases = copy_greatest_probability_genotypes(germline_genotypes, germline_normal_posteriors, max_germline_genotype_bases, 1e-100, 1e-2);
        latents.cancer_genotypes_.push_back(generate_all_cancer_genotypes(germline_bases, latents.indexed_haplotypes_, 1));
        if (latents.cancer_genotypes_.size() > 2 * max_allowed_cancer_genotypes) {
            if (!latents.cancer_genotype_prior_model_->mutation_model().is_primed()) {
                latents.cancer_genotype_prior_model_->mutation_model().prime(latents.haplotypes_.get());
            }
            const model::ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods};
            filter_with_germline_model(latents.cancer_genotypes_.back(), *latents.cancer_genotype_prior_model_,
                                       likelihood_model, samples_, max_allowed_cancer_genotypes);
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
    using GenotypeReference = std::reference_wrapper<const Genotype<IndexedHaplotype<>>>;
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
    assert(parameters_.max_genotypes);
    const auto num_haplotypes = latents.indexed_haplotypes_.size();
    const auto& germline_genotypes = latents.germline_genotypes_;
    const auto max_allowed_cancer_genotypes = *parameters_.max_genotypes;
    if (!latents.cancer_genotypes_.empty()) {
        const auto max_old_cancer_genotype_bases = std::max(max_allowed_cancer_genotypes / num_haplotypes, std::size_t {1});
        const auto& cancer_genotype_posteriors = latents.somatic_model_inferences_.back().max_evidence_params.genotype_probabilities;
        const auto old_cancer_genotype_bases = copy_greatest_probability_values(latents.cancer_genotypes_.back(), cancer_genotype_posteriors, max_old_cancer_genotype_bases);
        latents.cancer_genotypes_.push_back(extend_somatic(old_cancer_genotype_bases, latents.indexed_haplotypes_));
    } else {
        const auto max_germline_genotype_bases = calculate_max_germline_genotype_bases(max_allowed_cancer_genotypes, num_haplotypes, 1);
        const auto& germline_genotype_posteriors = latents.germline_model_inferences_.posteriors.genotype_probabilities;
        std::vector<double> germline_model_haplotype_posteriors(num_haplotypes);
        for (std::size_t g {0}; g < germline_genotypes.size(); ++g) {
            for (const auto& haplotype : collapse(germline_genotypes[g])) {
                germline_model_haplotype_posteriors[index_of(haplotype)] += germline_genotype_posteriors[g];
            }
        }
        const auto max_germline_haplotype_bases = max_num_elements(max_germline_genotype_bases, parameters_.ploidy);
        const auto top_haplotypes = copy_greatest_probability_values(latents.indexed_haplotypes_, germline_model_haplotype_posteriors,
                                                                     max_germline_haplotype_bases);
        auto germline_bases = generate_all_genotypes(top_haplotypes, parameters_.ploidy);
        latents.cancer_genotypes_.push_back(generate_all_cancer_genotypes(germline_bases, latents.indexed_haplotypes_, 1));
        if (latents.cancer_genotypes_.size() > 2 * max_allowed_cancer_genotypes) {
            if (!latents.cancer_genotype_prior_model_->mutation_model().is_primed()) {
                latents.cancer_genotype_prior_model_->mutation_model().prime(latents.haplotypes_);
            }
            const model::ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods};
            filter_with_germline_model(latents.cancer_genotypes_.back(), *latents.cancer_genotype_prior_model_,
                                       likelihood_model, samples_, max_allowed_cancer_genotypes);
        }
    }
}

void CancerCaller::generate_cancer_genotypes(Latents& latents, const MappableBlock<Genotype<IndexedHaplotype<>>>& germline_genotypes) const
{
    latents.cancer_genotypes_.push_back(generate_all_cancer_genotypes(germline_genotypes, latents.indexed_haplotypes_, latents.inferred_somatic_ploidy_));
}

bool CancerCaller::has_high_normal_contamination_risk(const Latents& latents) const
{
    return parameters_.normal_contamination_risk == Parameters::NormalContaminationRisk::high;
}

void CancerCaller::evaluate_germline_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(!(latents.indexed_haplotypes_.empty() || latents.germline_genotypes_.empty()));
    latents.germline_prior_model_ = make_germline_prior_model(latents.haplotypes_);
    latents.germline_model_ = std::make_unique<GermlineModel>(*latents.germline_prior_model_);
    const auto pooled_likelihoods = haplotype_likelihoods.merge_samples();
    latents.germline_prior_model_->prime(latents.haplotypes_);
    latents.germline_model_->prime(latents.haplotypes_);
    latents.germline_model_inferences_ = latents.germline_model_->evaluate(latents.germline_genotypes_, pooled_likelihoods);
}

void CancerCaller::evaluate_cnv_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(!latents.germline_genotypes_.empty() && latents.germline_prior_model_);
    auto cnv_model_priors = get_cnv_model_priors(*latents.germline_prior_model_);
    CNVModel::AlgorithmParameters params {};
    if (parameters_.max_vb_seeds) params.max_seeds = *parameters_.max_vb_seeds;
    params.target_max_memory = this->target_max_memory();
    CNVModel cnv_model {samples_, cnv_model_priors, params};
    cnv_model.prime(latents.haplotypes_);
    latents.cnv_model_inferences_ = cnv_model.evaluate(latents.germline_genotypes_,  haplotype_likelihoods);
}

void CancerCaller::evaluate_somatic_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(latents.germline_prior_model_ && !latents.cancer_genotypes_.empty() && !latents.cancer_genotypes_.back().empty());
    assert(latents.cancer_genotype_prior_model_);
    auto somatic_model_priors = get_somatic_model_priors(*latents.cancer_genotype_prior_model_, latents.inferred_somatic_ploidy_);
    SomaticModel::AlgorithmParameters params {};
    if (parameters_.max_vb_seeds) params.max_seeds = *parameters_.max_vb_seeds;
    params.target_max_memory = this->target_max_memory();
    SomaticModel model {samples_, somatic_model_priors, params};
    assert(latents.cancer_genotype_prior_model_->germline_model().is_primed());
    if (!latents.cancer_genotype_prior_model_->mutation_model().is_primed()) {
        latents.cancer_genotype_prior_model_->mutation_model().prime(latents.haplotypes_);
    }
    model.prime(latents.haplotypes_);
    latents.somatic_model_inferences_.push_back(model.evaluate(latents.cancer_genotypes_.back(), haplotype_likelihoods));
}

auto get_high_posterior_genotypes(const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& genotypes,
                                  const model::SomaticSubcloneModel::InferredLatents& latents)
{
    return copy_greatest_probability_values(genotypes, latents.max_evidence_params.genotype_probabilities, 10, 1e-3);
}

void CancerCaller::evaluate_noise_model(Latents& latents, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    if (has_normal_sample() && !has_high_normal_contamination_risk(latents)) {
        if (!latents.normal_germline_inferences_) {
            assert(latents.germline_model_);
            haplotype_likelihoods.prime(normal_sample());
            latents.normal_germline_inferences_ = latents.germline_model_->evaluate(latents.germline_genotypes_, haplotype_likelihoods);
        }
        assert(latents.cancer_genotype_prior_model_);
        auto noise_model_priors = get_noise_model_priors(*latents.cancer_genotype_prior_model_, latents.inferred_somatic_ploidy_);
        const SomaticModel noise_model {{*parameters_.normal_sample}, noise_model_priors};
        const auto& best_cancer_genotypes = latents.cancer_genotypes_[latents.max_evidence_somatic_model_index_];
        const auto& best_somatic_inferences = latents.somatic_model_inferences_[latents.max_evidence_somatic_model_index_];
        auto noise_genotypes = get_high_posterior_genotypes(best_cancer_genotypes, best_somatic_inferences);
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
    const auto& somatic_inferences  = latents.somatic_model_inferences_[latents.max_evidence_somatic_model_index_];
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
                            const Caller::Latents& latents,
                            OptionalThreadPool workers) const
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

struct VAFStats
{
    using CredibleRegion = std::pair<double, double>;
    CredibleRegion credible_region;
    double map, count;
};

auto compute_vaf_stats(const model::SomaticSubcloneModel::Priors::GenotypeMixturesDirichletAlphas& alphas,
                       const double credible_mass)
{
    const auto a0 = std::accumulate(std::cbegin(alphas), std::cend(alphas), 0.0);
    std::vector<VAFStats> result {};
    result.reserve(alphas.size());
    for (std::size_t i {0}; i < alphas.size(); ++i) {
        auto map_vaf = maths::dirichlet_expectation(i, alphas);
        auto vaf_cr = maths::beta_hdi(alphas[i], a0 - alphas[i], credible_mass);
        result.push_back({vaf_cr, map_vaf, alphas[i]});
    }
    return result;
}

using VAFStatsVector = std::vector<VAFStats>;
using VAFStatsMap = std::unordered_map<SampleName, VAFStatsVector>;

auto compute_vaf_stats(const model::SomaticSubcloneModel::Priors::GenotypeMixturesDirichletAlphaMap& alphas,
                       const double credible_mass)
{
    VAFStatsMap result {};
    result.reserve(alphas.size());
    for (const auto& p : alphas) {
        result.emplace(p.first, compute_vaf_stats(p.second, credible_mass));
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
    VAFStatsMap vaf_stats;
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

BigFloat marginalise(const Allele& allele, const MappableBlock<Genotype<IndexedHaplotype<>>>& genotypes, const std::vector<double>& probabilities)
{
    assert(genotypes.size() == probabilities.size());
    auto inv_result = std::inner_product(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(probabilities),
                                         0.0, std::plus<> {}, [&allele] (const auto& genotype, const auto probability) {
        return contains(genotype, allele) ? 0.0 : probability; });
    return BigFloat {1.0} - BigFloat {inv_result};
}

bool is_somatic(const Allele& allele, const CancerGenotype<IndexedHaplotype<>>& genotype)
{
    return contains(genotype.somatic(), allele) && !contains(genotype.germline(), allele);
}

BigFloat marginalise(const Allele& allele, const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& genotypes, const std::vector<double>& probabilities,
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
                                  const MappableBlock<Genotype<IndexedHaplotype<>>>& germline_genotypes,
                                  const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& cancer_genotypes,
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
                                  const MappableBlock<Genotype<IndexedHaplotype<>>>& germline_genotypes,
                                  const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& cancer_genotypes,
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

bool contains_alt(const Genotype<IndexedHaplotype<>>& genotype_call, const VariantReference& candidate)
{
    return includes(genotype_call, candidate.get().alt_allele());
}

auto call_candidates(const VariantPosteriorVector& candidate_posteriors,
                     const Genotype<IndexedHaplotype<>>& genotype_call,
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

BigFloat marginalise_somatic(const Allele& allele, const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& genotypes,
                             const std::vector<double>& probabilities)
{
    BigFloat result_complement {std::inner_product(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(probabilities),
                                                   0.0, std::plus<> {}, [&allele] (const auto& genotype, auto probability) {
        return is_somatic(allele, genotype) ? 0.0 : probability; })};
    return BigFloat {1.0} - result_complement;
}

auto compute_somatic_variant_posteriors(const std::vector<VariantReference>& candidates,
                                        const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& cancer_genotypes,
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
                                        const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& cancer_genotypes,
                                        const std::vector<double>& cancer_genotype_posteriors,
                                        const Phred<double> somatic_posterior)
{
    return compute_somatic_variant_posteriors(candidates, cancer_genotypes, cancer_genotype_posteriors,
                                              BigFloat {somatic_posterior.probability_true().value});
}

auto call_somatic_variants(const VariantPosteriorVector& somatic_variant_posteriors,
                           const CancerGenotype<IndexedHaplotype<>>& called_genotype,
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
marginalise(const CancerGenotype<Allele>& genotype, const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& genotypes,
            const std::vector<double>& genotype_posteriors)
{
    auto p = std::inner_product(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(genotype_posteriors),
                                0.0, std::plus<> {},  [&genotype] (const auto& g, auto probability) {
                                    return contains(g, genotype) ? 0.0 : probability; });
    return probability_false_to_phred(p);
}

auto call_somatic_genotypes(const CancerGenotype<IndexedHaplotype<>>& called_genotype,
                            const std::vector<GenomicRegion>& called_somatic_regions,
                            const MappableBlock<CancerGenotype<IndexedHaplotype<>>>& genotypes,
                            const std::vector<double>& genotype_posteriors,
                            const VAFStatsMap& vaf_stats)
{
    CancerGenotypeCalls result {};
    result.reserve(called_somatic_regions.size());
    for (const auto& region : called_somatic_regions) {
        auto genotype_chunk = copy<Allele>(called_genotype, region);
        auto posterior = marginalise(genotype_chunk, genotypes, genotype_posteriors);
        result.emplace_back(std::move(genotype_chunk), posterior);
        result.back().vaf_stats = vaf_stats;
    }
    return result;
}

// output

octopus::VariantCall::GenotypeCall demote(GermlineGenotypeCall call)
{
    return octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
}

std::unique_ptr<octopus::VariantCall>
transform_germline_cnv_call(GermlineVariantCall&& variant_call, GermlineGenotypeCall&& genotype_call,
                            const std::vector<SampleName>& samples, const std::vector<SampleName>& somatic_samples)
{
    std::vector<std::pair<SampleName, Call::GenotypeCall>> genotypes {};
    CNVCall::SomaticHaplotypeMap somatic_haplotypes {};
    somatic_haplotypes.reserve(samples.size());
    const auto germline_ploidy = genotype_call.genotype.ploidy();
    const auto somatic_ploidy = genotype_call.somatic.ploidy();
    for (const auto& sample : samples) {
        if (std::find(std::cbegin(somatic_samples), std::cend(somatic_samples), sample) == std::cend(somatic_samples)) {
            genotypes.emplace_back(sample, demote(genotype_call));
            somatic_haplotypes.emplace(sample, CNVCall::SomaticHaplotypeVector(germline_ploidy));
        } else {
            auto copy = genotype_call;
            for (auto allele : genotype_call.somatic) copy.genotype.emplace(allele);
            genotypes.emplace_back(sample, demote(std::move(copy)));
            CNVCall::SomaticHaplotypeVector sample_somatic_haplotypes(germline_ploidy + somatic_ploidy);
            std::fill_n(std::rbegin(sample_somatic_haplotypes), somatic_ploidy, true);
            somatic_haplotypes.emplace(sample, std::move(sample_somatic_haplotypes));
        }
    }
    return std::make_unique<CNVCall>(variant_call.variant.get(), std::move(genotypes),
                                     variant_call.segregation_quality, variant_call.posterior,
                                     std::move(somatic_haplotypes));
}

std::unique_ptr<octopus::VariantCall>
transform_germline_call(GermlineVariantCall&& variant_call, GermlineGenotypeCall&& genotype_call,
                        const std::vector<SampleName>& samples)
{
    std::vector<std::pair<SampleName, Call::GenotypeCall>> genotypes {};
    for (const auto& sample : samples) {
        genotypes.emplace_back(sample, demote(genotype_call));
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
                       SomaticCall::GenotypeStatsMap genotype_stats {};
                       genotype_stats.reserve(genotype_call.vaf_stats.size()); // num samples
                       for (const auto& p : genotype_call.vaf_stats) {
                           const VAFStatsVector& stats {p.second};
                           SomaticCall::GenotypeAlleleStats sample_stats {};
                           sample_stats.germline.reserve(genotype_call.genotype.germline_ploidy());
                           const auto convert_stats = [] (const VAFStats& stats) -> SomaticCall::AlleleStats {
                               return {stats.credible_region, stats.map, stats.count};
                           };
                           std::transform(std::cbegin(stats), std::next(std::cbegin(stats), genotype_call.genotype.germline_ploidy()),
                                          std::back_inserter(sample_stats.germline), convert_stats);
                           if (std::find(std::cbegin(somatic_samples), std::cend(somatic_samples), p.first) != std::cend(somatic_samples)) {
                               sample_stats.somatic.reserve(genotype_call.genotype.somatic_ploidy());
                               std::transform(std::next(std::cbegin(stats), genotype_call.genotype.germline_ploidy()), std::cend(stats),
                                              std::back_inserter(sample_stats.somatic), convert_stats);
                           }
                           genotype_stats.emplace(p.first, std::move(sample_stats));
                       }
                       return std::make_unique<SomaticCall>(variant_call.variant.get(), std::move(genotype_call.genotype),
                                                            variant_call.segregation_quality, genotype_call.posterior,
                                                            std::move(genotype_stats), variant_call.posterior);
                   });
    return result;
}

template <typename Map>
auto compute_posterior(const Genotype<Allele>& genotype, const Map& genotype_posteriors)
{
    auto p = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), 0.0,
                             [&genotype] (const double curr, const auto& p) {
                                 return curr + (contains(p.first, genotype) ? 0.0 : p.second);
                             });
    return probability_false_to_phred(p);
}

} // namespace

Phred<double>
CancerCaller::calculate_segregation_probability(const Variant& variant, const Latents& latents, double somatic_mass) const
{
    return octopus::calculate_segregation_probability(variant.alt_allele(),
                                                      latents.germline_genotypes_,
                                                      latents.cancer_genotypes_[latents.max_evidence_somatic_model_index_],
                                                      latents.germline_model_inferences_.posteriors.genotype_probabilities,
                                                      latents.cnv_model_inferences_.max_evidence_params.genotype_probabilities,
                                                      latents.somatic_model_inferences_[latents.max_evidence_somatic_model_index_].max_evidence_params.genotype_probabilities,
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
    const auto& best_somatic_model_inferences = latents.somatic_model_inferences_[latents.max_evidence_somatic_model_index_];
    const auto& cancer_genotype_posteriors = best_somatic_model_inferences.max_evidence_params.genotype_probabilities;
    const auto& best_cancer_genotypes_set = latents.cancer_genotypes_[latents.max_evidence_somatic_model_index_];
    log(latents.germline_genotypes_, germline_genotype_posteriors, latents.germline_model_inferences_, latents.cnv_model_inferences_,
        best_cancer_genotypes_set, best_somatic_model_inferences);
    const auto germline_candidate_posteriors = compute_candidate_posteriors(candidates, germline_genotype_posteriors);
    boost::optional<Genotype<IndexedHaplotype<>>> called_germline_genotype {};
    boost::optional<CancerGenotype<IndexedHaplotype<>>> called_cancer_genotype {};
    if (model_posteriors.somatic > model_posteriors.germline && somatic_posterior >= parameters_.min_somatic_posterior) {
        if (debug_log_) *debug_log_ << "Using cancer genotype for germline genotype call";
        if (!called_cancer_genotype) {
            auto cancer_posteriors = zip_cref(best_cancer_genotypes_set, cancer_genotype_posteriors);
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
    Genotype<IndexedHaplotype<>> called_somatic_genotype {};
    std::vector<SampleName> somatic_samples {};
    if (somatic_posterior >= parameters_.min_somatic_posterior) {
        auto somatic_allele_posteriors = compute_somatic_variant_posteriors(uncalled_germline_candidates, best_cancer_genotypes_set,
                                                                            cancer_genotype_posteriors, somatic_posterior);
        if (!called_cancer_genotype) {
            auto cancer_posteriors = zip_cref(best_cancer_genotypes_set, cancer_genotype_posteriors);
            called_cancer_genotype = find_map_genotype(cancer_posteriors)->first.get();
        }
        if (called_cancer_genotype->germline() == called_germline_genotype) {
            auto somatic_variant_calls = call_somatic_variants(somatic_allele_posteriors, *called_cancer_genotype,
                                                               parameters_.min_somatic_posterior);
            const auto& somatic_alphas = best_somatic_model_inferences.max_evidence_params.alphas;
            const auto vaf_stats = compute_vaf_stats(somatic_alphas, parameters_.credible_mass);
            if (!somatic_variant_calls.empty()) {
                for (const auto& p : vaf_stats) {
                    const auto& sample = p.first;
                    const auto& sample_vaf_stats = p.second;
                    if (debug_log_) {
                        auto ss = stream(*debug_log_);
                        ss << sample << " somatic credible regions: ";
                        for (auto stats : sample_vaf_stats) ss << '(' << stats.credible_region.first << ' ' << stats.credible_region.second << ") ";
                    }
                    if (std::any_of(std::next(std::cbegin(sample_vaf_stats), parameters_.ploidy), std::cend(sample_vaf_stats),
                        [this] (const auto& stats) { return stats.credible_region.first >= parameters_.min_credible_somatic_frequency; })) {
                        if (has_normal_sample() && sample == normal_sample()) {
                            somatic_samples.clear();
                            break;
                        }
                        somatic_samples.push_back(sample);
                    }
                }
                if (latents.noise_model_inferences_ && latents.normal_germline_inferences_) {
                    const auto noise_model_evidence = latents.noise_model_inferences_->approx_log_evidence;
                    const auto germline_model_evidence = latents.normal_germline_inferences_->log_evidence;
                    if (noise_model_evidence > germline_model_evidence) {
                        // Does the normal sample contain the called somatic variant?
                        const auto& noisy_alphas = latents.noise_model_inferences_->max_evidence_params.alphas.at(normal_sample());
                        const auto noise_mass = compute_credible_somatic_mass(noisy_alphas, latents.inferred_somatic_ploidy_, parameters_.min_expected_somatic_frequency);
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
            const auto called_somatic_regions = extract_regions(somatic_variant_calls);
            auto cancer_genotype_calls = call_somatic_genotypes(*called_cancer_genotype, called_somatic_regions,
                                                                best_cancer_genotypes_set, cancer_genotype_posteriors,
                                                                vaf_stats);
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
        const auto posterior = compute_posterior(genotype_chunk, germline_genotype_posteriors);
        if (called_somatic_genotype.ploidy() > 0) {
            germline_genotype_calls.emplace_back(std::move(genotype_chunk),
                                                 copy<Allele>(called_somatic_genotype, region),
                                                 posterior);
        } else {
            germline_genotype_calls.emplace_back(std::move(genotype_chunk), posterior);
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
                       if (somatic_samples.empty()) {
                           return transform_germline_call(std::move(variant_call), std::move(genotype_call), samples_);
                       } else {
                           return transform_germline_cnv_call(std::move(variant_call), std::move(genotype_call),
                                                              samples_, somatic_samples);
                       }
                       
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
    const auto& cnv_posteriors = latents.cnv_model_inferences_.max_evidence_params.genotype_probabilities;
    for (std::size_t i {0}; i < latents.germline_genotypes_.size(); ++i) {
        result[germline_genotypes[i]] += model_posteriors.cnv * cnv_posteriors[i];
    }
    const auto& cancer_genotypes = latents.cancer_genotypes_[latents.max_evidence_somatic_model_index_];
    const auto& somatic_posteriors = latents.somatic_model_inferences_[latents.max_evidence_somatic_model_index_].max_evidence_params.genotype_probabilities;
    for (std::size_t i {0}; i < cancer_genotypes.size(); ++i) {
        result[cancer_genotypes[i].germline()] += model_posteriors.somatic * somatic_posteriors[i];
    }
    return result;
}

double CancerCaller::calculate_somatic_mass(const CancerCaller::Latents& latents) const
{
    return compute_credible_somatic_mass(latents.somatic_model_inferences_[latents.max_evidence_somatic_model_index_].max_evidence_params.alphas,
                                         latents.inferred_somatic_ploidy_,
                                         parameters_.min_expected_somatic_frequency);
}

std::vector<std::unique_ptr<ReferenceCall>>
CancerCaller::call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents, const ReadPileupMap& pileups) const
{
    return {};
}

std::unique_ptr<GenotypePriorModel> CancerCaller::make_germline_prior_model(const HaplotypeBlock& haplotypes) const
{
    if (parameters_.germline_prior_model_params) {
        return std::make_unique<CoalescentGenotypePriorModel>(CoalescentModel {
        Haplotype {octopus::mapped_region(haplotypes), reference_},
        *parameters_.germline_prior_model_params
        });
    } else {
        return std::make_unique<UniformGenotypePriorModel>();
    }
}

// CancerCaller::Latents

CancerCaller::Latents::Latents(const HaplotypeBlock& haplotypes,
                               const std::vector<SampleName>& samples,
                               const CancerCaller::Parameters& parameters)
: haplotypes_ {haplotypes}
, indexed_haplotypes_ {index(haplotypes)}
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
    if (model_posteriors_.somatic < std::max(model_posteriors_.germline, model_posteriors_.cnv)) {
        // If no somatic variants likely then only consider the germline component
        GenotypeProbabilityMap genotype_posteriors {std::begin(germline_genotypes_), std::end(germline_genotypes_)};
        for (const auto& sample : samples_.get()) {
            insert_sample(sample, germline_model_inferences_.posteriors.genotype_probabilities, genotype_posteriors);
        }
        genotype_posteriors_ = std::make_shared<Latents::GenotypeProbabilityMap>(std::move(genotype_posteriors));
    } else {
        // If somatic variant(s) likely then consider all components
        auto total_num_genotypes = 2 * germline_genotypes_.size();
        for (const auto& genotypes : cancer_genotypes_) total_num_genotypes += genotypes.size();
        std::unordered_map<Genotype<IndexedHaplotype<>>, double> posteriors {};
        posteriors.reserve(total_num_genotypes);
        for (std::size_t g {0}; g < germline_genotypes_.size(); ++g) {
            posteriors[germline_genotypes_[g]]
                 += model_posteriors_.germline * germline_model_inferences_.posteriors.genotype_probabilities[g]
                  + model_posteriors_.cnv * cnv_model_inferences_.weighted_genotype_posteriors[g];
        }
        for (std::size_t somatic_model_idx {0}; somatic_model_idx < somatic_model_inferences_.size(); ++somatic_model_idx) {
            auto genotypes = demote_each(cancer_genotypes_[somatic_model_idx]);
            const auto& genotype_posteriors = somatic_model_inferences_[somatic_model_idx].weighted_genotype_posteriors;
            const auto somatic_model_posterior = model_posteriors_.somatic * somatic_model_posteriors_[somatic_model_idx];
            for (std::size_t g {0}; g < genotypes.size(); ++g) {
                posteriors[genotypes[g]] += somatic_model_posterior * genotype_posteriors[g];
            }
        }
        auto unique_genotypes = extract_keys(posteriors);
        auto unique_posteriors = extract_values(posteriors);
        GenotypeProbabilityMap genotype_posteriors {std::make_move_iterator(std::begin(unique_genotypes)),
                                                    std::make_move_iterator(std::end(unique_genotypes))};
        for (const auto& sample : samples_.get()) {
            insert_sample(sample, unique_posteriors, genotype_posteriors);
        }
        genotype_posteriors_ = std::make_shared<Latents::GenotypeProbabilityMap>(std::move(genotype_posteriors));
    }
}

void CancerCaller::Latents::compute_haplotype_posteriors() const
{
    Latents::HaplotypeProbabilityMap result {indexed_haplotypes_.size()};
    for (const auto& haplotype : indexed_haplotypes_) {
        result.emplace(haplotype, 0.0);
    }
    // Contribution from germline model
    for (const auto& p : zip(germline_genotypes_, germline_model_inferences_.posteriors.genotype_probabilities)) {
        for (const auto& haplotype : collapse(p.get<0>())) {
            result.at(haplotype) += model_posteriors_.germline * p.get<1>();
        }
    }
    // Contribution from CNV model
    for (const auto& p : zip(germline_genotypes_, cnv_model_inferences_.max_evidence_params.genotype_probabilities)) {
        for (const auto& haplotype : collapse(p.get<0>())) {
            result.at(haplotype) += model_posteriors_.cnv * p.get<1>();
        }
    }
    // Contribution from somatic model
    for (std::size_t somatic_model_idx {0}; somatic_model_idx < somatic_model_inferences_.size(); ++somatic_model_idx) {
        const auto& genotypes = cancer_genotypes_[somatic_model_idx];
        const auto& map_params = somatic_model_inferences_[somatic_model_idx].max_evidence_params;
        const auto conditional_somatic_prob = compute_credible_somatic_mass(map_params.alphas, somatic_model_idx + 1, parameters_.get().min_expected_somatic_frequency);
        const auto somatic_model_posterior = model_posteriors_.somatic * somatic_model_posteriors_[somatic_model_idx];
        const auto& genotype_posteriors = somatic_model_inferences_[somatic_model_idx].weighted_genotype_posteriors;
        for (const auto& p : zip(genotypes, genotype_posteriors)) {
            const auto unique_genotype = collapse(p.get<0>());
            for (const auto& haplotype : unique_genotype.germline()) {
                result.at(haplotype) += somatic_model_posterior * p.get<1>();
            }
            for (const auto& haplotype : unique_genotype.somatic()) {
                result.at(haplotype) += somatic_model_posterior * conditional_somatic_prob * p.get<1>();
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

namespace debug {

template <typename S, typename GenotypeReference>
void print_genotype_posteriors(S&& stream,
                               std::vector<std::pair<GenotypeReference, double>> genotype_posteriors,
                               const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    const auto m = std::min(n, genotype_posteriors.size());
    if (m == genotype_posteriors.size()) {
        stream << "Printing all genotype posteriors " << '\n';
    } else {
        stream << "Printing top " << m << " genotype posteriors " << '\n';
    }
    const auto mth = std::next(std::begin(genotype_posteriors), m);
    std::partial_sort(std::begin(genotype_posteriors), mth, std::end(genotype_posteriors),
                      [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
    std::for_each(std::begin(genotype_posteriors), mth,
                  [&] (const auto& p) {
                      print_variant_alleles(stream, p.first.get());
                      stream << " " << p.second << '\n';
                  });
}

} // namespace debug

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
        debug::print_variant_alleles(germline_log, map_germline->first.get());
        germline_log << ' ' << map_germline->second;
        auto cnv_posteriors = zip_cref(germline_genotypes, cnv_inferences.max_evidence_params.genotype_probabilities);
        auto map_cnv = find_map_genotype(cnv_posteriors);
        auto cnv_log = stream(*debug_log_);
        cnv_log << "MAP CNV genotype: ";
        debug::print_variant_alleles(cnv_log, map_cnv->first.get());
        cnv_log << ' ' << map_cnv->second;
        auto somatic_log = stream(*debug_log_);
        auto cancer_posteriors = zip_cref(cancer_genotypes, somatic_inferences.max_evidence_params.genotype_probabilities);
        auto map_somatic = find_map_genotype(cancer_posteriors);
        auto map_cancer_genotype = map_somatic->first.get();
        somatic_log << "MAP cancer genotype: ";
        auto weighted_cancer_posteriors = zip_cref(cancer_genotypes, somatic_inferences.weighted_genotype_posteriors);
        auto weighted_somatic_log = stream(*debug_log_);
        weighted_somatic_log << "Weighted cancer genotypes... ";
        debug::print_genotype_posteriors(weighted_somatic_log, weighted_cancer_posteriors, 10);
        debug::print_variant_alleles(somatic_log, map_cancer_genotype);
        somatic_log << ' ' << map_somatic->second;
        auto map_marginal_germline = find_map_genotype(germline_genotype_posteriors);
        auto marginal_germline_log = stream(*debug_log_);
        marginal_germline_log << "MAP marginal germline genotype: ";
        debug::print_variant_alleles(marginal_germline_log, map_marginal_germline->first);
        marginal_germline_log << ' ' << map_marginal_germline->second;
    }
    if (trace_log_) {
        auto weighted_cancer_posteriors = zip_cref(cancer_genotypes, somatic_inferences.weighted_genotype_posteriors);
        auto weighted_somatic_log = stream(*trace_log_);
        weighted_somatic_log << "Weighted cancer genotypes... ";
        debug::print_genotype_posteriors(weighted_somatic_log, weighted_cancer_posteriors);
    }
}

} // namespace octopus
