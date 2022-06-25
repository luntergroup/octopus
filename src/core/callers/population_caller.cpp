// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "population_caller.hpp"

#include <typeinfo>
#include <unordered_map>
#include <deque>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <utility>
#include <iostream>
#include <limits>

#include "basics/genomic_region.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "utils/maths.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/append.hpp"
#include "utils/select_top_k.hpp"
#include "containers/probability_matrix.hpp"
#include "core/models/genotype/individual_model.hpp"
#include "core/models/genotype/uniform_population_prior_model.hpp"
#include "core/models/genotype/coalescent_population_prior_model.hpp"
#include "core/models/genotype/uniform_genotype_prior_model.hpp"
#include "core/models/genotype/coalescent_genotype_prior_model.hpp"
#include "logging/logging.hpp"
#include "core/types/calls/germline_variant_call.hpp"
#include "core/types/calls/reference_call.hpp"

namespace octopus {

template <typename Container>
bool all_equal(const Container& c)
{
    return std::adjacent_find(std::cbegin(c), std::cend(c), std::not_equal_to<typename Container::value_type> {}) == std::cend(c);
}

PopulationCaller::PopulationCaller(Caller::Components&& components,
                                   Caller::Parameters general_parameters,
                                   Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {specific_parameters}
{
    unique_ploidies_ = parameters_.ploidies;
    std::sort(std::begin(unique_ploidies_), std::end(unique_ploidies_));
    unique_ploidies_.erase(std::unique(std::begin(unique_ploidies_), std::end(unique_ploidies_)), std::end(unique_ploidies_));
    if (unique_ploidies_.size() == 1 && unique_ploidies_[0] == 0) {
        throw std::logic_error {"PopulationCaller: at least one sample must have ploidy > 0"};
    }
}

std::string PopulationCaller::do_name() const
{
    return "population";
}

PopulationCaller::CallTypeSet PopulationCaller::do_call_types() const
{
    return {std::type_index(typeid(GermlineVariantCall))};
}

unsigned PopulationCaller::do_min_callable_ploidy() const
{
    return unique_ploidies_.front();
}

unsigned PopulationCaller::do_max_callable_ploidy() const
{
    return unique_ploidies_.back();
}

std::size_t PopulationCaller::do_remove_duplicates(HaplotypeBlock& haplotypes) const
{
    if (parameters_.deduplicate_haplotypes_with_germline_model) {
        if (haplotypes.size() < 2) return 0;
        CoalescentModel::Parameters model_params {};
        if (parameters_.prior_model_params) model_params = *parameters_.prior_model_params;
        Haplotype reference {mapped_region(haplotypes), reference_.get()};
        CoalescentModel model {std::move(reference), model_params, haplotypes.size(), CoalescentModel::CachingStrategy::none};
        const CoalescentProbabilityGreater cmp {std::move(model)};
        return octopus::remove_duplicates(haplotypes, cmp);
    } else {
        return Caller::do_remove_duplicates(haplotypes);
    }
}

// IndividualCaller::Latents public methods

namespace {

using InverseGenotypeTable = std::vector<std::vector<std::size_t>>;

auto make_inverse_genotype_table(const MappableBlock<Genotype<IndexedHaplotype<>>>& genotypes,
                                 const std::size_t num_haplotypes)
{
    InverseGenotypeTable result(num_haplotypes);
    for (auto& indices : result) indices.reserve(genotypes.size() / num_haplotypes);
    for (std::size_t genotype_idx {0}; genotype_idx < genotypes.size(); ++genotype_idx) {
        for (const auto& haplotype : genotypes[genotype_idx]) {
            result[index_of(haplotype)].push_back(genotype_idx);
        }
    }
    for (auto& indices : result) {
        std::sort(std::begin(indices), std::end(indices));
        indices.erase(std::unique(std::begin(indices), std::end(indices)), std::end(indices));
        indices.shrink_to_fit();
    }
    return result;
}

using GM = model::PopulationModel;

using GenotypeMarginalPosteriorVector = std::vector<double>;
using GenotypeMarginalPosteriorMatrix = std::vector<GenotypeMarginalPosteriorVector>;

auto calculate_haplotype_posteriors(const MappableBlock<IndexedHaplotype<>>& haplotypes,
                                    const MappableBlock<Genotype<IndexedHaplotype<>>>& genotypes,
                                    const GenotypeMarginalPosteriorMatrix& genotype_posteriors)
{
    const auto inverse_genotypes = make_inverse_genotype_table(genotypes, haplotypes.size());
    std::vector<double> result(haplotypes.size());
    auto itr = std::cbegin(inverse_genotypes);
    std::vector<std::size_t> genotype_indices(genotypes.size());
    std::iota(std::begin(genotype_indices), std::end(genotype_indices), 0);
    // noncontaining genotypes are genotypes that do not contain a particular haplotype.
    const auto num_noncontaining_genotypes = genotypes.size() - itr->size();
    std::vector<std::size_t> noncontaining_genotype_indices(num_noncontaining_genotypes);
    for (const auto& haplotype : haplotypes) {
        std::set_difference(std::cbegin(genotype_indices), std::cend(genotype_indices),
                            std::cbegin(*itr), std::cend(*itr),
                            std::begin(noncontaining_genotype_indices));
        double prob_not_observed {1};
        for (const auto& sample_genotype_posteriors : genotype_posteriors) {
            prob_not_observed *= std::accumulate(std::cbegin(noncontaining_genotype_indices), std::cend(noncontaining_genotype_indices),
                                                 0.0, [&] (const auto total, const auto i) {
                return total + sample_genotype_posteriors[i];
            });
        }
        result[index_of(haplotype)] = 1.0 - prob_not_observed;
        ++itr;
    }
    return result;
}

} // namespace

PopulationCaller::Latents::Latents(const std::vector<SampleName>& samples,
                                   const IndexedHaplotypeBlock& haplotypes,
                                   GenotypeBlock genotypes,
                                   IndependenceModelInferences&& inferences)
: genotypes_ {std::move(genotypes)}
{
    auto& genotype_marginal_posteriors = inferences.posteriors.genotype_probabilities;
    auto haplotype_posteriors = calculate_haplotype_posteriors(haplotypes, genotypes_, genotype_marginal_posteriors);
    haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>();
    haplotype_posteriors_->reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) {
        haplotype_posteriors_->emplace(haplotype, haplotype_posteriors[index_of(haplotype)]);
    }
    GenotypeProbabilityMap genotype_posteriors {std::begin(genotypes_), std::end(genotypes_)};
    for (std::size_t s {0}; s < samples.size(); ++s) {
        insert_sample(samples[s], std::move(genotype_marginal_posteriors[s]), genotype_posteriors);
    }
    genotype_posteriors_ = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
}

PopulationCaller::Latents::Latents(const std::vector<SampleName>& samples,
                                   const IndexedHaplotypeBlock& haplotypes,
                                   GenotypeBlock genotypes,
                                   ModelInferences&& inferences)
: genotypes_ {std::move(genotypes)}
, model_latents_ {std::move(inferences)}
{
    auto& genotype_marginal_posteriors = model_latents_.posteriors.marginal_genotype_probabilities;
    auto haplotype_posteriors = calculate_haplotype_posteriors(haplotypes, genotypes_, genotype_marginal_posteriors);
    haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>();
    for (const auto& haplotype : haplotypes) {
        haplotype_posteriors_->emplace(haplotype, haplotype_posteriors[index_of(haplotype)]);
    }
    GenotypeProbabilityMap genotype_posteriors {std::begin(genotypes_), std::end(genotypes_)};
    for (std::size_t s {0}; s < samples.size(); ++s) {
        insert_sample(samples[s], genotype_marginal_posteriors[s], genotype_posteriors);
    }
    genotype_posteriors_ = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
}

std::shared_ptr<PopulationCaller::Latents::HaplotypeProbabilityMap>
PopulationCaller::Latents::haplotype_posteriors() const noexcept
{
    return haplotype_posteriors_;
}

std::shared_ptr<PopulationCaller::Latents::GenotypeProbabilityMap>
PopulationCaller::Latents::genotype_posteriors() const noexcept
{
    return genotype_posteriors_;
}

using GenotypeBlock = MappableBlock<Genotype<IndexedHaplotype<>>>;
using GenotypeBlockReference = std::reference_wrapper<const GenotypeBlock>;
using GenotypesMap = std::map<unsigned, GenotypeBlock>;

std::unique_ptr<PopulationCaller::Caller::Latents>
PopulationCaller::infer_latents(const HaplotypeBlock& haplotypes,
                                const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                OptionalThreadPool workers) const
{
    if (use_independence_model()) {
        return infer_latents_with_independence_model(haplotypes, haplotype_likelihoods);
    } else {
        return infer_latents_with_joint_model(haplotypes, haplotype_likelihoods);
    }
}

boost::optional<Caller::ModelPosterior>
PopulationCaller::calculate_model_posterior(const HaplotypeBlock& haplotypes,
                                            const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                            const Caller::Latents& latents) const
{
    return calculate_model_posterior(haplotypes, haplotype_likelihoods, dynamic_cast<const Latents&>(latents));
}

static auto calculate_model_posterior(const double normal_model_log_evidence,
                                      const double dummy_model_log_evidence)
{
    constexpr double normalModelPrior {0.9999999};
    constexpr double dummyModelPrior {1.0 - normalModelPrior};
    const auto normal_model_ljp = std::log(normalModelPrior) + normal_model_log_evidence;
    const auto dummy_model_ljp  = std::log(dummyModelPrior) + dummy_model_log_evidence;
    const auto norm = maths::log_sum_exp(normal_model_ljp, dummy_model_ljp);
    return std::exp(normal_model_ljp - norm);
}

boost::optional<Caller::ModelPosterior>
PopulationCaller::calculate_model_posterior(const HaplotypeBlock& haplotypes,
                                            const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                            const Latents& latents) const
{
    const auto indexed_haplotypes = index(haplotypes);
    const auto prior_model = make_independent_prior_model(haplotypes);
    prior_model->prime(haplotypes);
    const model::IndividualModel model {*prior_model, debug_log_};
    ModelPosterior result {};
    result.samples.resize(samples_.size());
    for (std::size_t sample_idx {0}; sample_idx < samples_.size(); ++sample_idx) {
        haplotype_likelihoods.prime(samples_[sample_idx]);
        const auto genotypes = propose_model_check_genotypes(sample_idx, haplotypes, indexed_haplotypes, latents);
        const auto inferences1 = model.evaluate(genotypes.first, haplotype_likelihoods);
        const auto inferences2 = model.evaluate(genotypes.second, haplotype_likelihoods);
        result.samples[sample_idx] = octopus::calculate_model_posterior(inferences1.log_evidence, inferences2.log_evidence); 
    }
    return result;
}

//namespace debug {
//
//template <typename S>
//void print_genotype_posteriors(S&& stream, const GenotypeProbabilityMap& genotype_posteriors,
//                               std::size_t n = 5);
//void print_genotype_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
//                               std::size_t n = 5);
//template <typename S>
//void print_candidate_posteriors(S&& stream, const VariantPosteriors& candidate_posteriors,
//                                std::size_t n = 10);
//void print_candidate_posteriors(const VariantPosteriors& candidate_posteriors,
//                                std::size_t n = 10);
//void print_variant_calls(const VariantCallBlocks& calls);
//
//} // namespace debug

std::vector<std::unique_ptr<octopus::VariantCall>>
PopulationCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents, OptionalThreadPool workers) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

namespace {

using PopulationGenotypeProbabilityMap = ProbabilityMatrix<Genotype<IndexedHaplotype<>>>;
using GenotypeProbabilityMap = PopulationGenotypeProbabilityMap::InnerMap;

using VariantReference = std::reference_wrapper<const Variant>;
using VariantPosteriorVector = std::vector<std::pair<VariantReference, std::vector<Phred<double>>>>;

struct VariantCall : Mappable<VariantCall>
{
    VariantCall() = delete;
    VariantCall(const std::pair<VariantReference, std::vector<Phred<double>>>& p)
    : variant {p.first}
    , posteriors {p.second}
    {}
    VariantCall(const Variant& variant, std::vector<Phred<double>> posterior)
    : variant {variant}
    , posteriors {posterior}
    {}
    
    const GenomicRegion& mapped_region() const noexcept
    {
        return octopus::mapped_region(variant.get());
    }
    
    VariantReference variant;
    std::vector<Phred<double>> posteriors;
};

using VariantCalls = std::vector<VariantCall>;

struct GenotypeCall
{
    Genotype<Allele> genotype;
    Phred<double> posterior;
};

using GenotypeCalls = std::vector<std::vector<GenotypeCall>>;

// allele posterior calculations

using AlleleBools           = std::deque<bool>; // using std::deque because std::vector<bool> is evil
using GenotypePropertyBools = std::vector<AlleleBools>;

auto marginalise(const GenotypeProbabilityMap& genotype_posteriors,
                 const AlleleBools& contained_alleles)
{
    auto p = std::inner_product(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                std::cbegin(contained_alleles), 0.0, std::plus<> {},
                                [] (const auto& p, const bool is_contained) {
                                    return is_contained ? 0.0 : p.second;
                                });
    return probability_false_to_phred(p);
}

auto compute_sample_allele_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
                                      const GenotypePropertyBools& contained_alleles)
{
    std::vector<Phred<double>> result {};
    result.reserve(contained_alleles.size());
    for (const auto& allele : contained_alleles) {
        result.emplace_back(marginalise(genotype_posteriors, allele));
    }
    return result;
}

auto get_contained_alleles(const PopulationGenotypeProbabilityMap& genotype_posteriors,
                           const std::vector<Allele>& alleles)
{
    const auto num_genotypes = genotype_posteriors.size2();
    GenotypePropertyBools result {};
    if (num_genotypes == 0 || genotype_posteriors.empty1() || alleles.empty()) {
        return result;
    }
    result.reserve(alleles.size());
    const auto& test_sample   = genotype_posteriors.begin()->first;
    const auto genotype_begin = genotype_posteriors.begin(test_sample);
    const auto genotype_end   = genotype_posteriors.end(test_sample);
    for (const auto& allele : alleles) {
        result.emplace_back(num_genotypes);
        transform_contains(genotype_begin, genotype_end, allele, std::begin(result.back()),
                           [] (const auto& p, bool contain) { return contain; },
                           [] (const auto& p) { return p.first; });
    }
    return result;
}

auto compute_posteriors(const std::vector<SampleName>& samples,
                        const std::vector<Allele>& alleles,
                        const PopulationGenotypeProbabilityMap& genotype_posteriors)
{
    const auto contained_alleles = get_contained_alleles(genotype_posteriors, alleles);
    std::vector<std::vector<Phred<double>>> result {};
    result.reserve(genotype_posteriors.size1());
    for (const auto& sample : samples) {
        result.emplace_back(compute_sample_allele_posteriors(genotype_posteriors[sample], contained_alleles));
    }
    return result;
}

auto extract_ref_alleles(const std::vector<Variant>& variants)
{
    std::vector<Allele> result {};
    result.reserve(variants.size());
    std::transform(std::cbegin(variants), std::cend(variants), std::back_inserter(result),
                   [] (const auto& variant) { return variant.ref_allele(); });
    return result;
}

auto extract_alt_alleles(const std::vector<Variant>& variants)
{
    std::vector<Allele> result {};
    result.reserve(variants.size());
    std::transform(std::cbegin(variants), std::cend(variants), std::back_inserter(result),
                   [] (const auto& variant) { return variant.alt_allele(); });
    return result;
}

auto compute_posteriors(const std::vector<SampleName>& samples,
                        const std::vector<Variant>& variants,
                        const PopulationGenotypeProbabilityMap& genotype_posteriors)
{
    const auto allele_posteriors = compute_posteriors(samples, extract_alt_alleles(variants), genotype_posteriors);
    VariantPosteriorVector result {};
    result.reserve(variants.size());
    for (std::size_t i {0}; i < variants.size(); ++i) {
        std::vector<Phred<double>> sample_posteriors(samples.size());
        std::transform(std::cbegin(allele_posteriors), std::cend(allele_posteriors), std::begin(sample_posteriors),
                       [i] (const auto& ps) { return ps[i]; });
        result.emplace_back(variants[i], std::move(sample_posteriors));
    }
    return result;
}

// haplotype genotype calling

//auto call_genotypes(const GM::Latents& latents, const std::vector<Genotype<IndexedHaplotype<>>>& genotypes)
//{
//    const auto itr = std::max_element(std::cbegin(latents.joint_genotype_probabilities), std::cend(latents.joint_genotype_probabilities));
//    const auto& called_indices = latents.genotype_combinations[std::distance(std::cbegin(latents.joint_genotype_probabilities), itr)];
//    std::vector<Genotype<IndexedHaplotype<>>> result {};
//    result.reserve(called_indices.size());
//    std::transform(std::cbegin(called_indices), std::cend(called_indices), std::back_inserter(result),
//                   [&] (auto idx) { return genotypes[idx]; });
//    return result;
//}
//
//auto call_genotypes(const GM::Latents& latents, const std::unordered_map<unsigned, std::vector<Genotype<IndexedHaplotype<>>>>& genotypes)
//{
//    return call_genotypes(latents, std::cbegin(genotypes)->second);
//}

auto call_genotype(const PopulationGenotypeProbabilityMap::InnerMap& genotype_posteriors)
{
    return std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                            [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; })->first;
}

auto call_genotypes(const std::vector<SampleName>& samples, const PopulationGenotypeProbabilityMap& genotype_posteriors)
{
    std::vector<Genotype<IndexedHaplotype<>>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.push_back(call_genotype(genotype_posteriors[sample]));
    }
    return result;
}

// variant calling

bool has_above(const std::vector<Phred<double>>& posteriors, const Phred<double> min_posterior)
{
    return std::any_of(std::cbegin(posteriors), std::cend(posteriors), [=] (auto p) { return p >= min_posterior; });
}

bool contains_alt(const Genotype<IndexedHaplotype<>>& genotype_call, const VariantReference& candidate)
{
    return includes(genotype_call, candidate.get().alt_allele());
}

bool contains_alt(const std::vector<Genotype<IndexedHaplotype<>>>& genotype_calls, const VariantReference& candidate)
{
    return std::any_of(std::cbegin(genotype_calls), std::cend(genotype_calls),
                       [&] (const auto& genotype) { return contains_alt(genotype, candidate); });
}

VariantCalls call_candidates(const VariantPosteriorVector& candidate_posteriors,
                             const std::vector<Genotype<IndexedHaplotype<>>>& genotype_calls,
                             const Phred<double> min_posterior)
{
    VariantCalls result {};
    result.reserve(candidate_posteriors.size());
    std::copy_if(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors),
                 std::back_inserter(result),
                 [&genotype_calls, min_posterior] (const auto& p) {
                     return has_above(p.second, min_posterior) && contains_alt(genotype_calls, p.first);
                 });
    return result;
}

// polymorphism calculations

auto get_homozygous_alleles(const PopulationGenotypeProbabilityMap& genotype_posteriors,
                            const std::vector<Allele>& alleles)
{
    const auto num_genotypes = genotype_posteriors.size2();
    GenotypePropertyBools result {};
    if (num_genotypes == 0 || genotype_posteriors.empty1() || alleles.empty()) {
        return result;
    }
    result.reserve(alleles.size());
    const auto& test_sample   = genotype_posteriors.begin()->first;
    const auto genotype_begin = genotype_posteriors.begin(test_sample);
    const auto genotype_end   = genotype_posteriors.end(test_sample);
    for (const auto& allele : alleles) {
        result.emplace_back(num_genotypes);
        transform_is_homozygous(genotype_begin, genotype_end, allele, std::begin(result.back()),
                                [] (const auto& p, bool is_hom) { return is_hom; },
                                [] (const auto& p) { return p.first; });
    }
    return result;
}

auto marginalise_homozygous(const GenotypeProbabilityMap& genotype_posteriors,
                            const AlleleBools& homozygotes)
{
    return std::inner_product(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                              std::cbegin(homozygotes), 0.0, std::plus<> {},
                              [] (const auto& p, const bool is_homozygous) {
                                  return is_homozygous ? p.second : 0.0;
                              });
}

auto marginalise_homozygous(const std::vector<SampleName>& samples,
                            const PopulationGenotypeProbabilityMap& genotype_posteriors,
                            const AlleleBools& homozygotes)
{
    std::vector<double> ps(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(ps),
                   [&] (const auto& sample) {
                       return std::log(marginalise_homozygous(genotype_posteriors[sample], homozygotes));
                   });
    const auto p = std::accumulate(std::cbegin(ps), std::cend(ps), 0.0);
    return probability_false_to_phred(std::exp(p));
}

auto compute_homozygous_posteriors(const std::vector<Allele>& alleles,
                                   const std::vector<SampleName>& samples,
                                   const PopulationGenotypeProbabilityMap& genotype_posteriors)
{
    const auto homozygous_alleles = get_homozygous_alleles(genotype_posteriors, alleles);
    std::vector<Phred<double>> result(alleles.size());
    std::transform(std::cbegin(alleles), std::cend(alleles), std::cbegin(homozygous_alleles), std::begin(result),
                   [&] (const auto& allele, const auto& homozygotes) {
                       return marginalise_homozygous(samples, genotype_posteriors, homozygotes);
                   });
    return result;
}

auto compute_polymorphism_posteriors(const std::vector<Variant>& variants,
                                     const std::vector<SampleName>& samples,
                                     const PopulationGenotypeProbabilityMap& genotype_posteriors)
{
    return compute_homozygous_posteriors(extract_ref_alleles(variants), samples, genotype_posteriors);
}

// allele genotype calling

auto marginalise(const Genotype<Allele>& genotype, const GenotypeProbabilityMap& genotype_posteriors)
{
    double mass_not_contained {0};
    if (genotype.ploidy() > 0) {
        for_each_contains(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), genotype,
                          [&] (const auto& p, bool contain) { if (!contain) mass_not_contained += p.second; },
                          [] (const auto& p) { return p.first; });
    }
    return probability_false_to_phred(mass_not_contained);
}

auto call_genotypes(const std::vector<SampleName>& samples,
                    const std::vector<Genotype<IndexedHaplotype<>>>& genotype_calls,
                    const PopulationGenotypeProbabilityMap& genotype_posteriors,
                    const std::vector<GenomicRegion>& variant_regions)
{
    GenotypeCalls result {};
    result.reserve(variant_regions.size());
    for (const auto& region : variant_regions) {
        std::vector<GenotypeCall> region_calls {};
        region_calls.reserve(samples.size());
        for (std::size_t s {0}; s < samples.size(); ++s) { 
            auto genotype_chunk = copy<Allele>(genotype_calls[s], region);
            const auto posterior = marginalise(genotype_chunk, genotype_posteriors[samples[s]]);
            region_calls.push_back({std::move(genotype_chunk), posterior});
        }
        result.push_back(std::move(region_calls));
    }
    return result;
}

// output

octopus::VariantCall::GenotypeCall convert(GenotypeCall&& call)
{
    return octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
}

auto max(const std::vector<Phred<double>>& posteriors)
{
    return *std::max_element(std::cbegin(posteriors), std::cend(posteriors));
}

std::unique_ptr<octopus::VariantCall>
transform_call(const std::vector<SampleName>& samples,
               VariantCall&& variant_call,
               std::vector<GenotypeCall>&& sample_genotype_calls)
{
    std::vector<std::pair<SampleName, Call::GenotypeCall>> tmp {};
    tmp.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples),
                   std::make_move_iterator(std::begin(sample_genotype_calls)),
                   std::back_inserter(tmp),
                   [] (const auto& sample, auto&& genotype) {
                       return std::make_pair(sample, convert(std::move(genotype)));
                   });
    auto p = std::accumulate(std::cbegin(variant_call.posteriors), std::cend(variant_call.posteriors), 0.0,
                             [] (auto curr, auto x) { return curr + x.score(); });
    return std::make_unique<GermlineVariantCall>(variant_call.variant.get(), std::move(tmp), Phred<> {p});
}

auto transform_calls(const std::vector<SampleName>& samples,
                     VariantCalls&& variant_calls,
                     GenotypeCalls&& genotype_calls)
{
    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
    result.reserve(variant_calls.size());
    std::transform(std::make_move_iterator(std::begin(variant_calls)), std::make_move_iterator(std::end(variant_calls)),
                   std::make_move_iterator(std::begin(genotype_calls)), std::back_inserter(result),
                   [&samples] (auto&& variant_call, auto&& genotype_call) {
                       return transform_call(samples, std::move(variant_call), std::move(genotype_call));
                   });
    return result;
}

} // namespace

namespace debug {
void log(const PopulationGenotypeProbabilityMap& genotype_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log);
void log(const VariantPosteriorVector& candidate_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log);
void log(const Genotype<IndexedHaplotype<>>& called_genotype,
         boost::optional<logging::DebugLogger>& debug_log);
} // namespace debug

std::vector<std::unique_ptr<octopus::VariantCall>>
PopulationCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    const auto& genotype_posteriors = *latents.genotype_posteriors_;
    debug::log(genotype_posteriors, debug_log_, trace_log_);
    const auto candidate_posteriors = compute_posteriors(samples_, candidates, genotype_posteriors);
    debug::log(candidate_posteriors, debug_log_, trace_log_);
    const auto genotype_calls = call_genotypes(samples_, genotype_posteriors);
    auto variant_calls = call_candidates(candidate_posteriors, genotype_calls, parameters_.min_variant_posterior);
    const auto called_regions = extract_regions(variant_calls);
    auto allele_genotype_calls = call_genotypes(samples_, genotype_calls, genotype_posteriors, called_regions);
    return transform_calls(samples_, std::move(variant_calls), std::move(allele_genotype_calls));
}

namespace {

// reference genotype calling

struct RefCall : public Mappable<RefCall>
{
    RefCall() = default;
    
    template <typename A>
    RefCall(A&& reference_allele, double posterior)
    : reference_allele {std::forward<A>(reference_allele)}
    , posterior {posterior}
    {}
    
    const GenomicRegion& mapped_region() const noexcept { return reference_allele.mapped_region(); }
    
    Allele reference_allele;
    double posterior;
};

using RefCalls = std::vector<RefCall>;

//    double marginalise_reference_genotype(const Allele& reference_allele,
//                                          const GenotypeProbabilityMap& sample_genotype_posteriors)
//    {
//        double result {0};
//        
//        for (const auto& genotype_posterior : sample_genotype_posteriors) {
//            if (is_homozygous(genotype_posterior.first, reference_allele)) {
//                result += genotype_posterior.second;
//            }
//        }
//        
//        return result;
//    }
    
//    RefCalls call_reference(const GenotypeProbabilityMap& genotype_posteriors,
//                            const std::vector<Allele>& reference_alleles,
//                            const ReadMap::mapped_type& reads, const double min_call_posterior)
//    {
//        RefCalls result {};
//        
//        if (reference_alleles.empty()) return result;
//        
//        result.reserve(reference_alleles.size());
//        
//        for (const auto& reference_allele : reference_alleles) {
//            double posterior {0};
//            
//            if (has_coverage(reads, mapped_region(reference_allele))) {
//                posterior = marginalise_reference_genotype(reference_allele,
//                                                           genotype_posteriors);
//            }
//            
//            if (posterior >= min_call_posterior) {
//                result.emplace_back(reference_allele, posterior);
//            }
//        }
//        
//        result.shrink_to_fit();
//        
//        return result;
//    }
} // namespace

std::vector<std::unique_ptr<ReferenceCall>>
PopulationCaller::call_reference(const std::vector<Allele>& alleles,
                                 const Caller::Latents& latents,
                                 const ReadPileupMap& pileups) const
{
    return {};
}

bool PopulationCaller::use_independence_model() const noexcept
{
    return parameters_.use_independent_genotype_priors || !parameters_.prior_model_params;
}

auto generate_unique_genotypes(const MappableBlock<IndexedHaplotype<>>& haplotypes, std::vector<unsigned> ploidies)
{
    const auto region = mapped_region(haplotypes);
    std::sort(std::begin(ploidies), std::end(ploidies));
    ploidies.erase(std::unique(std::begin(ploidies), std::end(ploidies)), std::end(ploidies));
    GenotypesMap result {};
    for (auto ploidy : ploidies) {
        result[ploidy] = {generate_all_genotypes(haplotypes, ploidy), region};
    }
    return result;
}

namespace {

template <typename Container1, typename MappableType, typename Container2>
auto append(MappableBlock<MappableType, Container1>&& src, MappableBlock<MappableType, Container2>& dst)
{
    return utils::append(std::move(static_cast<Container1&>(src)), static_cast<Container2&>(dst));
}

} // namespace append

std::unique_ptr<Caller::Latents>
PopulationCaller::infer_latents_with_joint_model(const HaplotypeBlock& haplotypes,
                                                 const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto indexed_haplotypes = index(haplotypes);
    const auto prior_model = make_joint_prior_model(haplotypes);
    const model::PopulationModel model {*prior_model, {parameters_.max_genotype_combinations}, debug_log_};
    prior_model->prime(haplotypes);
    if (unique_ploidies_.size() == 1) {
        auto genotypes = generate_all_genotypes(indexed_haplotypes, parameters_.ploidies.front());
        if (debug_log_) stream(*debug_log_) << "There are " << genotypes.size() << " candidate genotypes";
        auto inferences = model.evaluate(samples_, haplotypes, genotypes, haplotype_likelihoods);
        return std::make_unique<Latents>(samples_, indexed_haplotypes, std::move(genotypes), std::move(inferences));
    } else {
        model::PopulationModel::GenotypeVector genotypes {};
        for (const auto ploidy : unique_ploidies_) {
            if (ploidy > 0) {
                append(generate_all_genotypes(indexed_haplotypes, ploidy), genotypes);
            } else {
                genotypes.push_back(Genotype<IndexedHaplotype<>> {});
            }
        }
        auto inferences = model.evaluate(samples_, parameters_.ploidies, haplotypes, genotypes, haplotype_likelihoods);
        return std::make_unique<Latents>(samples_, indexed_haplotypes, std::move(genotypes), std::move(inferences));
    }
}

std::unique_ptr<Caller::Latents>
PopulationCaller::infer_latents_with_independence_model(const HaplotypeBlock& haplotypes,
                                                        const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto indexed_haplotypes = index(haplotypes);
    const auto prior_model = make_independent_prior_model(haplotypes);
    prior_model->prime(haplotypes);
    const model::IndependentPopulationModel model {*prior_model, debug_log_};
    if (parameters_.ploidies.size() == 1) {
        auto genotypes = generate_all_genotypes(indexed_haplotypes, parameters_.ploidies.front());
        if (debug_log_) stream(*debug_log_) << "There are " << genotypes.size() << " candidate genotypes";
        auto inferences = model.evaluate(samples_, genotypes, haplotype_likelihoods);
        return std::make_unique<Latents>(samples_, indexed_haplotypes, std::move(genotypes), std::move(inferences));
    } else {
        model::PopulationModel::GenotypeVector genotypes {};
        for (const auto ploidy : unique_ploidies_) {
            if (ploidy > 0) {
                append(generate_all_genotypes(indexed_haplotypes, ploidy), genotypes);
            } else {
                genotypes.push_back(Genotype<IndexedHaplotype<>> {});
            }
        }
        auto inferences = model.evaluate(samples_, parameters_.ploidies, genotypes, haplotype_likelihoods);
        return std::make_unique<Latents>(samples_, indexed_haplotypes, std::move(genotypes), std::move(inferences));
    }
}

std::unique_ptr<PopulationPriorModel> PopulationCaller::make_joint_prior_model(const HaplotypeBlock& haplotypes) const
{
    if (parameters_.prior_model_params) {
        return std::make_unique<CoalescentPopulationPriorModel>(CoalescentModel{
        Haplotype {mapped_region(haplotypes), reference_},
        *parameters_.prior_model_params
        });
    } else {
        return std::make_unique<UniformPopulationPriorModel>();
    }
}

std::unique_ptr<GenotypePriorModel> PopulationCaller::make_independent_prior_model(const HaplotypeBlock& haplotypes) const
{
    if (parameters_.prior_model_params) {
        return std::make_unique<CoalescentGenotypePriorModel>(CoalescentModel {
        Haplotype {mapped_region(haplotypes), reference_},
        *parameters_.prior_model_params, haplotypes.size(), CoalescentModel::CachingStrategy::address
        });
    } else {
        return std::make_unique<UniformGenotypePriorModel>();
    }
}

namespace {

template <typename IndexType>
void erase_duplicates(MappableBlock<Genotype<IndexedHaplotype<IndexType>>>& genotypes)
{
    using std::begin; using std::end;
    std::sort(begin(genotypes), end(genotypes), GenotypeLess {});
    genotypes.erase(std::unique(begin(genotypes), end(genotypes)), end(genotypes));
}

} // namespace

std::pair<PopulationCaller::GenotypeBlock, PopulationCaller::GenotypeBlock>
PopulationCaller::propose_model_check_genotypes(std::size_t sample_idx,
                                                const HaplotypeBlock& haplotypes,
                                                const IndexedHaplotypeBlock& indexed_haplotypes,
                                                const Latents& latents) const
{
    constexpr std::size_t max_model_check_genotypes {5};
    const auto& genotypes = latents.genotypes_;
    const auto& marginal_genotype_posteriors = latents.model_latents_.posteriors.marginal_genotype_probabilities[sample_idx];
    const auto num_model_check_genotypes = std::min(max_model_check_genotypes, genotypes.size());
    const auto best_indices = select_top_k_indices(marginal_genotype_posteriors, num_model_check_genotypes, false);
    GenotypeBlock assumed {mapped_region(haplotypes)};
    for (const auto genotype_idx : best_indices) {
        assumed.push_back(genotypes[genotype_idx]);
    }
    auto augmented = extend(assumed, indexed_haplotypes);
    erase_duplicates(augmented);
    return {std::move(assumed), std::move(augmented)};
}

namespace debug {

namespace {

template <typename S>
void print_genotype_posteriors(S&& stream,
                               const GenotypeProbabilityMap& genotype_posteriors,
                               const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    const auto m = std::min(n, genotype_posteriors.size());
    using GenotypeReference = std::reference_wrapper<const Genotype<IndexedHaplotype<>>>;
    std::vector<std::pair<GenotypeReference, double>> v {};
    v.reserve(genotype_posteriors.size());
    std::copy(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
              std::back_inserter(v));
    const auto mth = std::next(std::begin(v), m);
    std::partial_sort(std::begin(v), mth, std::end(v),
                      [](const auto& lhs, const auto& rhs) {
                          return lhs.second > rhs.second;
                      });
    std::for_each(std::begin(v), mth,
                  [&](const auto& p) {
                      print_variant_alleles(stream, p.first.get());
                      stream << " " << p.second << '\n';
                  });
}

template <typename S>
void print_genotype_posteriors(S&& stream,
                               const PopulationGenotypeProbabilityMap& genotype_posteriors,
                               const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    for (const auto& p : genotype_posteriors) {
        stream << "Printing genotype posteriors for sample: " << p.first << '\n';
        print_genotype_posteriors(stream, p.second, n);
    }
}

void print_genotype_posteriors(const PopulationGenotypeProbabilityMap& genotype_posteriors,
                               const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    print_genotype_posteriors(std::cout, genotype_posteriors, n);
}

template <typename S>
void print_candidate_posteriors(S&& stream, const VariantPosteriorVector& candidate_posteriors,
                                const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    const auto m = std::min(n, candidate_posteriors.size());
    if (m == candidate_posteriors.size()) {
        stream << "Printing all candidate variant posteriors " << '\n';
    } else {
        stream << "Printing top " << m << " candidate variant posteriors " << '\n';
    }
    std::vector<std::pair<VariantReference, std::vector<Phred<double>>>> v {};
    v.reserve(candidate_posteriors.size());
    std::copy(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors),
              std::back_inserter(v));
    const auto mth = std::next(std::begin(v), m);
    std::partial_sort(std::begin(v), mth, std::end(v),
                      [] (const auto& lhs, const auto& rhs) {
                          return lhs.second > rhs.second;
                      });
    std::for_each(std::begin(v), mth,
                  [&] (const auto& p) {
                      stream << p.first.get() << " ";
                      for (auto x : p.second) {
                          stream << x.probability_true() << ' ';
                      }
                      stream << '\n';
                  });
}

void print_candidate_posteriors(const VariantPosteriorVector& candidate_posteriors,
                                const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    print_candidate_posteriors(std::cout, candidate_posteriors, n);
}
    
} // namespace

void log(const PopulationGenotypeProbabilityMap& genotype_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log)
{
    if (trace_log) {
        print_genotype_posteriors(stream(*trace_log), genotype_posteriors);
    }
    if (debug_log) {
        print_genotype_posteriors(stream(*debug_log), genotype_posteriors, 5);
    }
}

void log(const VariantPosteriorVector& candidate_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log)
{
    if (trace_log) {
        print_candidate_posteriors(stream(*trace_log), candidate_posteriors);
    }
    if (debug_log) {
        print_candidate_posteriors(stream(*debug_log), candidate_posteriors, 5);
    }
}
    
} // namespace debug
} // namespace octopus
