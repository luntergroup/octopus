// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "polyclone_caller.hpp"

#include <typeinfo>
#include <unordered_map>
#include <deque>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <limits>

#include "basics/genomic_region.hpp"
#include "containers/probability_matrix.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/calls/polyclone_variant_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "core/models/genotype/uniform_genotype_prior_model.hpp"
#include "core/models/genotype/coalescent_genotype_prior_model.hpp"
#include "core/models/mutation/somatic_mutation_model.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/concat.hpp"
#include "utils/select_top_k.hpp"
#include "utils/set_partition.hpp"
#include "logging/logging.hpp"

namespace octopus {

PolycloneCaller::PolycloneCaller(Caller::Components&& components,
                                 Caller::Parameters general_parameters,
                                 Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {std::move(specific_parameters)}
{
    if (parameters_.max_clones < 1) {
        throw std::logic_error {"PolycloneCaller: max_clones must be > 1"};
    }
    if (parameters_.max_clones > model::SubcloneModel::max_ploidy) {
        static std::atomic_bool warned {false};
        if (!warned) {
            warned = true;
            logging::WarningLogger log {};
            stream(log) << "Maximum supported clonality is "
                            << model::SubcloneModel::max_ploidy
                            << " but " << parameters_.max_clones << " was requested";
        }
        parameters_.max_clones = model::SubcloneModel::max_ploidy;
    }
}

std::string PolycloneCaller::do_name() const
{
    return "polyclone";
}

PolycloneCaller::CallTypeSet PolycloneCaller::do_call_types() const
{
    return {std::type_index(typeid(PolycloneVariantCall))};
}

unsigned PolycloneCaller::do_min_callable_ploidy() const
{
    return 1;
}

unsigned PolycloneCaller::do_max_callable_ploidy() const
{
    return parameters_.max_clones;
}

std::size_t PolycloneCaller::do_remove_duplicates(HaplotypeBlock& haplotypes) const
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

// PolycloneCaller::Latents public methods

PolycloneCaller::Latents::Latents(GenotypeBlock haploid_genotypes,
                                  GenotypeBlock polyploid_genotypes,
                                  HaploidModelInferences haploid_model_inferences,
                                  SubloneModelInferences subclone_model_inferences,
                                  const SampleName& sample,
                                  const std::function<double(unsigned)>& clonality_prior)
: haploid_genotypes_ {std::move(haploid_genotypes)}
, polyploid_genotypes_ {std::move(polyploid_genotypes)}
, haploid_model_inferences_ {std::move(haploid_model_inferences)}
, subclone_model_inferences_ {std::move(subclone_model_inferences)}
, model_log_posteriors_ {0, std::numeric_limits<double>::min()}
, sample_ {sample}
{
    if (!polyploid_genotypes_.empty()) {
        const auto haploid_model_prior = std::log(clonality_prior(1));
        const auto called_subclonality = polyploid_genotypes_.front().ploidy();
        const auto subclone_model_prior = std::log(clonality_prior(called_subclonality));
        const auto haploid_model_jp = haploid_model_prior + haploid_model_inferences_.log_evidence;
        const auto subclone_model_jp = subclone_model_prior + subclone_model_inferences_.approx_log_evidence;
        const auto norm = maths::log_sum_exp({haploid_model_jp, subclone_model_jp});
        model_log_posteriors_ = {haploid_model_jp - norm, subclone_model_jp - norm};
        auto log_posteriors = concat(haploid_model_inferences_.posteriors.genotype_log_probabilities,
                                     subclone_model_inferences_.max_evidence_params.genotype_log_probabilities);
        std::for_each(std::begin(log_posteriors), std::next(std::begin(log_posteriors), haploid_genotypes_.size()),
                      [&] (auto& p) { p += model_log_posteriors_.clonal; });
        std::for_each(std::next(std::begin(log_posteriors), haploid_genotypes_.size()), std::end(log_posteriors),
                      [&] (auto& p) { p += model_log_posteriors_.subclonal; });
        auto genotypes = concat(haploid_genotypes_, polyploid_genotypes_);
        genotype_log_posteriors_ = std::make_shared<GenotypeProbabilityMap>(std::make_move_iterator(std::begin(genotypes)),
                                                                            std::make_move_iterator(std::end(genotypes)));
        insert_sample(sample_, log_posteriors, *genotype_log_posteriors_);
    }
}

std::shared_ptr<PolycloneCaller::Latents::HaplotypeProbabilityMap>
PolycloneCaller::Latents::haplotype_posteriors() const noexcept
{
    if (haplotype_posteriors_ == nullptr) {
        haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>();
        for (const auto& p : (*(this->genotype_posteriors()))[sample_]) {
            for (const auto& haplotype : collapse(p.first)) {
                (*haplotype_posteriors_)[haplotype] += p.second;
            }
        }
    }
    return haplotype_posteriors_;
}

std::shared_ptr<PolycloneCaller::Latents::GenotypeProbabilityMap>
PolycloneCaller::Latents::genotype_posteriors() const noexcept
{
    if (genotype_posteriors_ == nullptr) {
        auto genotypes = concat(haploid_genotypes_, polyploid_genotypes_);
        auto posteriors = concat(haploid_model_inferences_.posteriors.genotype_probabilities,
                                 subclone_model_inferences_.weighted_genotype_posteriors);
        const ModelProbabilities model_posterior {std::exp( model_log_posteriors_.clonal), std::exp( model_log_posteriors_.subclonal)};
        std::for_each(std::begin(posteriors), std::next(std::begin(posteriors), haploid_genotypes_.size()),
                      [&] (auto& p) { p *= model_posterior.clonal; });
        std::for_each(std::next(std::begin(posteriors), haploid_genotypes_.size()), std::end(posteriors),
                      [=] (auto& p) { p *= model_posterior.subclonal; });
        genotype_posteriors_ = std::make_shared<GenotypeProbabilityMap>(std::make_move_iterator(std::begin(genotypes)),
                                                                        std::make_move_iterator(std::end(genotypes)));
        insert_sample(sample_, posteriors, *genotype_posteriors_);
    }
    return genotype_posteriors_;
}

// PolycloneCaller::Latents private methods

namespace {

auto make_sublone_model_mixture_prior_map(const SampleName& sample, const unsigned num_clones, const double alpha = 1.0)
{
    model::SubcloneModel::Priors::GenotypeMixturesDirichletAlphaMap result {};
    model::SubcloneModel::Priors::GenotypeMixturesDirichletAlphas alphas(num_clones, alpha);
    result.emplace(sample, std::move(alphas));
    return result;
}

}

std::unique_ptr<PolycloneCaller::Caller::Latents>
PolycloneCaller::infer_latents(const HaplotypeBlock& haplotypes, 
                               const HaplotypeLikelihoodArray& haplotype_likelihoods,
                               OptionalThreadPool workers) const
{
    const auto indexed_haplotypes = index(haplotypes);
    auto haploid_genotypes = generate_all_genotypes(indexed_haplotypes, 1);
    if (debug_log_) stream(*debug_log_) << "There are " << haploid_genotypes.size() << " candidate haploid genotypes";
    auto genotype_prior_model = make_prior_model(haplotypes);
    genotype_prior_model->prime(haplotypes);
    const model::IndividualModel haploid_model {*genotype_prior_model, debug_log_};
    haplotype_likelihoods.prime(sample());
    auto haploid_inferences = haploid_model.evaluate(haploid_genotypes, haplotype_likelihoods);
    if (debug_log_) stream(*debug_log_) << "Evidence for haploid model is " << haploid_inferences.log_evidence;
    GenotypeBlock polyploid_genotypes {};
    model::SubcloneModel::InferredLatents sublonal_inferences;
    fit_sublone_model(haplotypes, indexed_haplotypes, haplotype_likelihoods, *genotype_prior_model, haploid_inferences,
                      polyploid_genotypes, sublonal_inferences, workers);
    if (debug_log_) stream(*debug_log_) << "There are " << polyploid_genotypes.size() << " candidate polyploid genotypes";
    if (parameters_.haplogroup_prior && sublonal_inferences.approx_log_evidence > haploid_inferences.log_evidence) {
        fit_haplogroup_model(haplotypes, indexed_haplotypes, haplotype_likelihoods, polyploid_genotypes, *genotype_prior_model, sublonal_inferences);
    }
    using std::move;
    return std::make_unique<Latents>(move(haploid_genotypes), move(polyploid_genotypes),
                                     move(haploid_inferences), move(sublonal_inferences),
                                     sample(), parameters_.clonality_prior);
}

std::vector<std::unique_ptr<octopus::VariantCall>>
PolycloneCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents, OptionalThreadPool workers) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

namespace {

using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<IndexedHaplotype<>>>::InnerMap;
using VariantReference = std::reference_wrapper<const Variant>;
using VariantPosteriorVector = std::vector<std::pair<VariantReference, Phred<double>>>;

struct VariantCall : Mappable<VariantCall>
{
    VariantCall() = delete;
    VariantCall(const std::pair<VariantReference, Phred<double>>& p)
    : variant {p.first}
    , posterior {p.second}
    {}
    VariantCall(const Variant& variant, Phred<double> posterior)
    : variant {variant}
    , posterior {posterior}
    {}
    
    const GenomicRegion& mapped_region() const noexcept
    {
        return octopus::mapped_region(variant.get());
    }
    
    VariantReference variant;
    Phred<double> posterior;
    bool is_dummy_filtered = false;
};

using VariantCalls = std::vector<VariantCall>;

struct GenotypeCall
{
    template <typename T> GenotypeCall(T&& genotype, Phred<double> posterior)
    : genotype {std::forward<T>(genotype)}
    , posterior {posterior}
    {}
    
    Genotype<Allele> genotype;
    Phred<double> posterior;
};

using GenotypeCalls = std::vector<GenotypeCall>;

// allele posterior calculations

template <typename GenotypeOrAllele>
auto marginalise_contained(const GenotypeOrAllele& element, const GenotypeProbabilityMap& genotype_log_posteriors)
{
    thread_local std::vector<double> buffer {};
    buffer.clear();
    for_each_contains(std::cbegin(genotype_log_posteriors), std::cend(genotype_log_posteriors), element,
                      [&] (const auto& p, bool contain) { if (!contain) buffer.push_back(p.second); },
                      [] (const auto& p) { return p.first; });
    if (!buffer.empty()) {
        return log_probability_false_to_phred(std::min(maths::log_sum_exp(buffer), 0.0));
    } else {
        return Phred<> {std::numeric_limits<double>::infinity()};
    }
}

auto compute_candidate_posteriors(const std::vector<Variant>& candidates,
                                  const GenotypeProbabilityMap& genotype_log_posteriors)
{
    VariantPosteriorVector result {};
    result.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        result.emplace_back(candidate, marginalise_contained(candidate.alt_allele(), genotype_log_posteriors));
    }
    return result;
}

// variant calling

bool has_callable(const VariantPosteriorVector& variant_posteriors, const Phred<double> min_posterior) noexcept
{
    return std::any_of(std::cbegin(variant_posteriors), std::cend(variant_posteriors),
                       [=] (const auto& p) noexcept { return p.second >= min_posterior; });
}

bool contains_alt(const Genotype<IndexedHaplotype<>>& genotype_call, const VariantReference& candidate)
{
    return includes(genotype_call, candidate.get().alt_allele());
}

VariantCalls call_candidates(const VariantPosteriorVector& candidate_posteriors,
                             const Genotype<IndexedHaplotype<>>& genotype_call,
                             const Phred<double> min_posterior)
{
    VariantCalls result {};
    result.reserve(candidate_posteriors.size());
    std::copy_if(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors), std::back_inserter(result),
                 [&] (const auto& p) { return p.second >= min_posterior && contains_alt(genotype_call, p.first); });
    return result;
}

// variant genotype calling

template <typename PairIterator>
PairIterator find_map(PairIterator first, PairIterator last)
{
    return std::max_element(first, last, [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
}

auto call_genotype(const GenotypeProbabilityMap& genotype_posteriors, const bool ignore_hom_ref = false)
{
    const auto map_itr = find_map(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors));
    assert(map_itr != std::cend(genotype_posteriors));
    if (!ignore_hom_ref || !is_homozygous_reference(map_itr->first)) {
        return map_itr->first;
    } else {
        const auto lhs_map_itr = find_map(std::cbegin(genotype_posteriors), map_itr);
        const auto rhs_map_itr = find_map(std::next(map_itr), std::cend(genotype_posteriors));
        if (lhs_map_itr != map_itr) {
            if (rhs_map_itr != std::cend(genotype_posteriors)) {
                return lhs_map_itr->second < rhs_map_itr->second ? rhs_map_itr->first : lhs_map_itr->first;
            } else {
                return lhs_map_itr->first;
            }
        } else {
            return rhs_map_itr->first;
        }
    }
}

auto compute_posterior(const Genotype<Allele>& genotype, const GenotypeProbabilityMap& genotype_posteriors)
{
    auto p = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), 0.0,
                             [&] (const double curr, const auto& p) {
                                 return curr + (contains(p.first, genotype) ? 0.0 : p.second);
                             });
    return probability_false_to_phred(p);
}

GenotypeCalls call_genotypes(const Genotype<IndexedHaplotype<>>& genotype_call,
                             const GenotypeProbabilityMap& genotype_log_posteriors,
                             const std::vector<GenomicRegion>& variant_regions)
{
    GenotypeCalls result {};
    result.reserve(variant_regions.size());
    for (const auto& region : variant_regions) {
        auto genotype_chunk = copy<Allele>(genotype_call, region);
        const auto posterior = marginalise_contained(genotype_chunk, genotype_log_posteriors);
        result.emplace_back(std::move(genotype_chunk), posterior);
    }
    return result;
}

// output

using HaplotypeFrequencyStats = PolycloneVariantCall::HaplotypeFrequencyStats;
using HaplotypeFrequencyStatsVector = PolycloneVariantCall::HaplotypeFrequencyStatsVector;

auto compute_haplotype_frequency_stats(const model::SubcloneModel::Priors::GenotypeMixturesDirichletAlphas& alphas)
{
    HaplotypeFrequencyStatsVector result {};
    result.reserve(alphas.size());
    for (std::size_t i {0}; i < alphas.size(); ++i) {
        auto map_vaf = maths::dirichlet_expectation(i, alphas);
        result.push_back({alphas[i], map_vaf});
    }
    return result;
}

octopus::VariantCall::GenotypeCall convert(GenotypeCall&& call)
{
    return octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
}

std::unique_ptr<octopus::VariantCall>
transform_call(const SampleName& sample, VariantCall&& variant_call, GenotypeCall&& genotype_call, 
               const HaplotypeFrequencyStatsVector& haplotype_frequency_stats)
{
    std::vector<std::pair<SampleName, Call::GenotypeCall>> tmp {std::make_pair(sample, convert(std::move(genotype_call)))};
    std::unique_ptr<octopus::VariantCall> result {};
    if (haplotype_frequency_stats.empty()) {
        result = std::make_unique<PolycloneVariantCall>(variant_call.variant.get(), std::move(tmp), variant_call.posterior);
    } else {
        result = std::make_unique<PolycloneVariantCall>(variant_call.variant.get(), std::move(tmp), variant_call.posterior, haplotype_frequency_stats);
    }   
    return result;
}

auto transform_calls(const SampleName& sample, VariantCalls&& variant_calls, GenotypeCalls&& genotype_calls,
                    const HaplotypeFrequencyStatsVector& haplotype_frequency_stats)
{
    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
    result.reserve(variant_calls.size());
    std::transform(std::make_move_iterator(std::begin(variant_calls)), std::make_move_iterator(std::end(variant_calls)),
                   std::make_move_iterator(std::begin(genotype_calls)), std::back_inserter(result),
                   [&] (VariantCall&& variant_call, GenotypeCall&& genotype_call) {
                       return transform_call(sample, std::move(variant_call), std::move(genotype_call), haplotype_frequency_stats);
                   });
    return result;
}

} // namespace

namespace debug { namespace {

void log(const GenotypeProbabilityMap& genotype_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log);

void log(const VariantPosteriorVector& candidate_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log,
         Phred<double> min_posterior);

void log(const Genotype<IndexedHaplotype<>>& called_genotype,
         boost::optional<logging::DebugLogger>& debug_log);

} // namespace
} // namespace debug

std::vector<std::unique_ptr<octopus::VariantCall>>
PolycloneCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    log(latents);
    const auto& genotype_log_posteriors = (*latents.genotype_log_posteriors_)[sample()];
    debug::log(genotype_log_posteriors, debug_log_, trace_log_);
    const auto candidate_posteriors = compute_candidate_posteriors(candidates, genotype_log_posteriors);
    debug::log(candidate_posteriors, debug_log_, trace_log_, parameters_.min_variant_posterior);
    const bool force_call_non_ref {has_callable(candidate_posteriors, parameters_.min_variant_posterior)};
    const auto genotype_call = octopus::call_genotype(genotype_log_posteriors, force_call_non_ref);
    auto variant_calls = call_candidates(candidate_posteriors, genotype_call, parameters_.min_variant_posterior);
    const auto called_regions = extract_regions(variant_calls);
    auto genotype_calls = call_genotypes(genotype_call, genotype_log_posteriors, called_regions);
    HaplotypeFrequencyStatsVector haplotype_frequency_stats {};
    if (genotype_call.ploidy() > 1) {
        haplotype_frequency_stats = compute_haplotype_frequency_stats(latents.subclone_model_inferences_.max_evidence_params.alphas.at(sample()));
    }
    return transform_calls(sample(), std::move(variant_calls), std::move(genotype_calls), haplotype_frequency_stats);
}

std::vector<std::unique_ptr<ReferenceCall>>
PolycloneCaller::call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents, const ReadPileupMap& pileup) const
{
    return call_reference(alleles, dynamic_cast<const Latents&>(latents), pileup);
}

std::vector<std::unique_ptr<ReferenceCall>>
PolycloneCaller::call_reference(const std::vector<Allele>& alleles, const Latents& latents, const ReadPileupMap& pileup) const
{
    return {};
}

const SampleName& PolycloneCaller::sample() const noexcept
{
    return samples_.front();
}

namespace {

template <typename T>
T nth_greatest_value(std::vector<T> values, const std::size_t n)
{
    auto nth_itr = std::next(std::begin(values), n);
    std::nth_element(std::begin(values), nth_itr, std::end(values), std::greater<> {});
    return *nth_itr;
}

template <typename T>
void erase_indices(std::vector<T>& v, const std::vector<std::size_t>& indices)
{
    assert(std::is_sorted(std::cbegin(indices), std::cend(indices)));
    std::for_each(std::crbegin(indices), std::crend(indices), [&v] (auto idx) { v.erase(std::next(std::cbegin(v), idx)); });
}

template <typename GenotypeBlock>
void reduce(GenotypeBlock& genotypes,
            const std::vector<double>& genotype_probabilities,
            const std::size_t n)
{
    if (genotypes.size() <= n) return;
    const auto min_probability = nth_greatest_value(genotype_probabilities, n + 1);
    std::size_t idx {0};
    const auto is_low_probability = [&] (const auto& genotype) { return genotype_probabilities[idx++] <= min_probability; };
    genotypes.erase(std::remove_if(std::begin(genotypes), std::end(genotypes), is_low_probability), std::end(genotypes));
}

template <typename GenotypeBlock>
void reduce(GenotypeBlock& genotypes,
            const MappableBlock<Haplotype>& haplotypes,
            const GenotypePriorModel& genotype_prior_model,
            const HaplotypeLikelihoodArray& haplotype_likelihoods,
            const std::size_t n)
{
    if (genotypes.size() <= n) return;
    model::IndividualModel approx_model {genotype_prior_model};
    approx_model.prime(haplotypes);
    const auto approx_posteriors = approx_model.evaluate(genotypes, haplotype_likelihoods).posteriors.genotype_log_probabilities;
    reduce(genotypes, approx_posteriors, n);
}

auto make_hint(const std::size_t num_genotypes, const std::size_t idx, const double p = 0.9999)
{
    model::SubcloneModel::Latents::LogProbabilityVector result(num_genotypes, num_genotypes > 1 ? std::log((1 - p) / (num_genotypes - 1)) : 0);
    if (num_genotypes > 1) result[idx] = std::log(p);
    return result;
}

template <typename GenotypeBlock>
auto
propose_subclone_model_hints(const GenotypeBlock& curr_genotypes,
                             const MappableBlock<IndexedHaplotype<>>& haplotypes,
                             const model::IndividualModel::InferredLatents& haploid_latents,
                             const std::size_t max)
{
    const auto top_haploid_indices = select_top_k_indices(haploid_latents.posteriors.genotype_log_probabilities, max);
    std::vector<model::SubcloneModel::Latents::LogProbabilityVector> result {};
    result.reserve(max);
    for (std::size_t i {0}; i < haplotypes.size() && result.size() < max; ++i) {
        const auto& base = haplotypes[top_haploid_indices[i]];
        for (std::size_t j {0}; j < top_haploid_indices.size() && result.size() < max; ++j) {
            const auto& haplotype = haplotypes[top_haploid_indices[j]];
            if (i != j) {
                const Genotype<IndexedHaplotype<>> hint {base, haplotype};
                const auto hint_itr = std::find(std::cbegin(curr_genotypes), std::cend(curr_genotypes), hint);
                if (hint_itr != std::cend(curr_genotypes)) {
                    const auto hint_idx = static_cast<std::size_t>(std::distance(std::cbegin(curr_genotypes), hint_itr));
                    result.push_back(make_hint(curr_genotypes.size(), hint_idx));
                }
            }
        }
    }
    return result;
}

template <typename GenotypeBlock>
auto
propose_subclone_model_hints(const GenotypeBlock& curr_genotypes,
                             const GenotypeBlock& prev_genotypes,
                             const model::SubcloneModel::InferredLatents& sublonal_inferences,
                             const MappableBlock<IndexedHaplotype<>>& haplotypes,
                             const model::IndividualModel::InferredLatents& haploid_latents,
                             const std::size_t max)
{
    const auto top_prev_genotype_indices = select_top_k_indices(sublonal_inferences.weighted_genotype_posteriors, max);
    const auto top_haploid_indices = select_top_k_indices(haploid_latents.posteriors.genotype_log_probabilities, max);
    std::vector<model::SubcloneModel::Latents::LogProbabilityVector> result {};
    result.reserve(max);
    for (std::size_t i {0}; i < top_prev_genotype_indices.size() && result.size() < max; ++i) {
        const auto& base = prev_genotypes[top_prev_genotype_indices[i]];
        for (std::size_t j {0}; j < top_haploid_indices.size() && result.size() < max; ++j) {
            const auto& haplotype = haplotypes[top_haploid_indices[j]];
            if (!contains(base, haplotype)) {
                auto hint = base; hint.emplace(haplotype);
                const auto hint_itr = std::find(std::cbegin(curr_genotypes), std::cend(curr_genotypes), hint);
                if (hint_itr != std::cend(curr_genotypes)) {
                    const auto hint_idx = static_cast<std::size_t>(std::distance(std::cbegin(curr_genotypes), hint_itr));
                    result.push_back(make_hint(curr_genotypes.size(), hint_idx));
                }
            }
            if (i < top_prev_genotype_indices.size() - 1 && j < top_haploid_indices.size() - 1) {
                // Do we prefer to add a new haplotype to the current base or introduce a new base?
                const auto this_polyploid_prob = sublonal_inferences.weighted_genotype_posteriors[top_prev_genotype_indices[i]];
                const auto next_polyploid_prob = sublonal_inferences.weighted_genotype_posteriors[top_prev_genotype_indices[i + 1]];
                if (next_polyploid_prob > 0) {
                    const auto polyploid_ratio = this_polyploid_prob / next_polyploid_prob;
                    const auto this_haploid_prob = haploid_latents.posteriors.genotype_probabilities[top_haploid_indices[j]];
                    const auto next_haploid_prob = haploid_latents.posteriors.genotype_probabilities[top_haploid_indices[j + 1]];
                    const auto haploid_ratio = this_haploid_prob / next_haploid_prob;
                    if (next_haploid_prob == 0 || haploid_ratio < polyploid_ratio) {
                        break;
                    }
                }
            }
        }
    }
    return result;
}

} // namespace

void
PolycloneCaller::fit_sublone_model(const HaplotypeBlock& haplotypes,
                                   const IndexedHaplotypeBlock& indexed_haplotypes,
                                   const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                   GenotypePriorModel& genotype_prior_model,
                                   const model::IndividualModel::InferredLatents& haploid_latents,
                                   GenotypeBlock& prev_genotypes,
                                   model::SubcloneModel::InferredLatents& sublonal_inferences,
                                   OptionalThreadPool workers) const
{
    model::SubcloneModel::AlgorithmParameters model_params {};
    if (parameters_.max_vb_seeds) model_params.max_seeds = *parameters_.max_vb_seeds;
    model_params.target_max_memory = this->target_max_memory();
    GenotypeBlock curr_genotypes {};
    const auto haploid_log_prior = std::log(parameters_.clonality_prior(1));
    const auto max_clones = std::min(parameters_.max_clones, static_cast<unsigned>(haplotypes.size()));
    for (unsigned clonality {2}; clonality <= max_clones; ++clonality) {
        const auto clonal_model_prior = parameters_.clonality_prior(clonality);
        if (clonal_model_prior == 0.0) break;
        genotype_prior_model.unprime();
        genotype_prior_model.prime(haplotypes);
        const auto max_possible_genotypes = num_max_zygosity_genotypes_noexcept(haplotypes.size(), clonality);
        if (prev_genotypes.empty() || clonality <= 2 || !parameters_.max_genotypes ||
            (max_possible_genotypes && *max_possible_genotypes <= *parameters_.max_genotypes)) {
            curr_genotypes = generate_all_max_zygosity_genotypes(indexed_haplotypes, clonality);
        } else {
            const static auto not_included = [] (const auto& genotype, const auto& haplotype) { return !contains(genotype, haplotype); };
            if (prev_genotypes.size() * (haplotypes.size() / 2) > *parameters_.max_genotypes) {
                auto probable_prev_genotypes = prev_genotypes;
                reduce(probable_prev_genotypes, sublonal_inferences.max_evidence_params.genotype_log_probabilities,
                       *parameters_.max_genotypes / (haplotypes.size() / 2));
                curr_genotypes= extend(probable_prev_genotypes, indexed_haplotypes, not_included);
            } else {
                curr_genotypes = extend(prev_genotypes, indexed_haplotypes, not_included);
            }
        }
        if (parameters_.max_genotypes) reduce(curr_genotypes, haplotypes, genotype_prior_model, haplotype_likelihoods, *parameters_.max_genotypes);
        if (debug_log_) stream(*debug_log_) << "Generated " << curr_genotypes.size() << " genotypes with clonality " << clonality;
        if (curr_genotypes.empty()) break;
        model::SubcloneModel::Priors priors {genotype_prior_model, make_sublone_model_mixture_prior_map(sample(), clonality, parameters_.clone_mixture_prior_concentration)};
        model::SubcloneModel model {{sample()}, priors, model_params};
        model.prime(haplotypes);
        std::vector<model::SubcloneModel::Latents::LogProbabilityVector> hints {};
        if (clonality == 2) {
            hints = propose_subclone_model_hints(curr_genotypes, indexed_haplotypes, haploid_latents, model_params.max_seeds / 2);
        } else {
            hints = propose_subclone_model_hints(curr_genotypes, prev_genotypes, sublonal_inferences, indexed_haplotypes, haploid_latents, model_params.max_seeds / 2);
        }
        auto inferences = model.evaluate(curr_genotypes, haplotype_likelihoods, std::move(hints), workers);
        if (debug_log_) {
            stream(*debug_log_) << "Evidence for model with clonality " << clonality << " is " << inferences.approx_log_evidence;
            const auto max_genotype_posterior_itr = std::max_element(std::cbegin(inferences.weighted_genotype_posteriors), std::cend(inferences.weighted_genotype_posteriors));
            const std::size_t map_genotype_idx = std::distance(std::cbegin(inferences.weighted_genotype_posteriors), max_genotype_posterior_itr);
            auto debug_stream = stream(*debug_log_);
            debug_stream << "MAP genotype with clonality " << clonality << ": ";
            debug::print_variant_alleles(debug_stream, curr_genotypes[map_genotype_idx]);
            debug_stream << " " << *max_genotype_posterior_itr;
        }
        if (clonality == 2) {
            prev_genotypes = std::move(curr_genotypes);
            sublonal_inferences = std::move(inferences);
            if ((std::log(clonal_model_prior) + sublonal_inferences.approx_log_evidence)
                < (haploid_log_prior + haploid_latents.log_evidence)) break;
        } else {
            if ((std::log(clonal_model_prior) + inferences.approx_log_evidence)
                <= (std::log(parameters_.clonality_prior(clonality - 1)) + sublonal_inferences.approx_log_evidence))  break;
            prev_genotypes = std::move(curr_genotypes);
            sublonal_inferences = std::move(inferences);
        }
    }
}

void 
PolycloneCaller::fit_haplogroup_model(const HaplotypeBlock& haplotypes,
                                      const IndexedHaplotypeBlock& indexed_haplotypes,
                                      const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                      const GenotypeBlock& polyploid_genotypes,
                                      const GenotypePriorModel& genotype_prior_model,
                                      const model::SubcloneModel::InferredLatents& sublonal_inferences) const
{
    const auto called_clonality = polyploid_genotypes.front().ploidy();
    auto max_evidence = sublonal_inferences.approx_log_evidence;
    const std::size_t max_haplogrouped_genotype_candidates {10};
    const auto best_polyploid_genotype_indices = select_top_k_indices(sublonal_inferences.weighted_genotype_posteriors, max_haplogrouped_genotype_candidates);
    boost::optional<model::HaplogroupSubcloneModel::InferredLatents> called_haplogroup_latents {};
    std::vector<PartitionedGenotype<IndexedHaplotype<>>> best_haplogrouped_genotypes {};
    for (unsigned num_haplogroups {2}; num_haplogroups <= called_clonality; ++num_haplogroups) {
        std::vector<PartitionedGenotype<IndexedHaplotype<>>> haplogrouped_genotypes {};
        for (const auto& genotype_idx : best_polyploid_genotype_indices) {
            const auto& genotype = polyploid_genotypes[genotype_idx];
            set_partitions(called_clonality, num_haplogroups, [&genotype, &haplogrouped_genotypes] (const auto& partitions) {
                PartitionedGenotype<IndexedHaplotype<>> g {partitions.size()};
                for (std::size_t p {0}; p < partitions.size(); ++p) {
                    for (const auto haplotype_idx : partitions[p]) {
                        g.partition(p).emplace(genotype[haplotype_idx]);
                    }
                }
                haplogrouped_genotypes.push_back(std::move(g));
            });
        }
        if (debug_log_) stream(*debug_log_) << "Generated " << haplogrouped_genotypes.size() << " haplogrouped genotypes";
        SomaticMutationModel mutation_model {SomaticMutationModel::Parameters {1e-6, 1e-7}};
        mutation_model.prime(haplotypes);
        HaplogroupGenotypePriorModel haplogroup_prior_model {genotype_prior_model, mutation_model};
        model::HaplogroupSubcloneModel::Priors priors {haplogroup_prior_model, make_sublone_model_mixture_prior_map(sample(), called_clonality, parameters_.clone_mixture_prior_concentration)};
        model::HaplogroupSubcloneModel model {{sample()}, priors};
        auto haplogroup_latents = model.evaluate(haplogrouped_genotypes, haplotype_likelihoods, {});
        if (haplogroup_latents.approx_log_evidence <= max_evidence) break;
        called_haplogroup_latents = std::move(haplogroup_latents);
        best_haplogrouped_genotypes = std::move(haplogrouped_genotypes);
        max_evidence = haplogroup_latents.approx_log_evidence;
    }
    if (debug_log_) {
        if (called_haplogroup_latents) {
            const auto num_haplogroups = best_haplogrouped_genotypes.front().num_partitions();
            stream(*debug_log_) << num_haplogroups << " haplogroups called";
            auto top_haplogroup_indices = select_top_k_indices(called_haplogroup_latents->weighted_genotype_posteriors, 3);
            auto ss = stream(*debug_log_);
            ss << "Best haplogroups:\n";
            for (const auto& idx : top_haplogroup_indices) {
                ss << "\t"; debug::print_variant_alleles(ss, best_haplogrouped_genotypes[idx]); ss << '\n';
            }
        } else {
            *debug_log_ << "One haplogroup called";
        }
    }
}

std::unique_ptr<GenotypePriorModel> PolycloneCaller::make_prior_model(const HaplotypeBlock& haplotypes) const
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

void PolycloneCaller::log(const Latents& latents) const
{
    if (debug_log_) {
        stream(*debug_log_) << "Clonal model posterior is " << latents.model_log_posteriors_.clonal
                            << " and subclonal model posterior is " << latents.model_log_posteriors_.subclonal;
        if (latents.model_log_posteriors_.subclonal > latents.model_log_posteriors_.clonal) {
            stream(*debug_log_) << "Detected subclonality is " << latents.polyploid_genotypes_.front().ploidy();
        }
    }
}

namespace debug { namespace {

template <typename S>
void print_genotype_posteriors(S&& stream,
                               const GenotypeProbabilityMap& genotype_posteriors,
                               const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    const auto m = std::min(n, genotype_posteriors.size());
    if (m == genotype_posteriors.size()) {
        stream << "Printing all genotype posteriors " << '\n';
    } else {
        stream << "Printing top " << m << " genotype posteriors " << '\n';
    }
    using GenotypeReference = std::reference_wrapper<const Genotype<IndexedHaplotype<>>>;
    std::vector<std::pair<GenotypeReference, double>> v {};
    v.reserve(genotype_posteriors.size());
    std::copy(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), std::back_inserter(v));
    const auto mth = std::next(std::begin(v), m);
    std::partial_sort(std::begin(v), mth, std::end(v),
                      [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
    std::for_each(std::begin(v), mth, [&] (const auto& p) {
        print_variant_alleles(stream, p.first.get());
        stream << " " << p.second << '\n';
    });
}

void print_genotype_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
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
    std::vector<std::pair<VariantReference, Phred<double>>> v {};
    v.reserve(candidate_posteriors.size());
    std::copy(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors), std::back_inserter(v));
    const auto mth = std::next(std::begin(v), m);
    std::partial_sort(std::begin(v), mth, std::end(v),
                      [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
    std::for_each(std::begin(v), mth, [&] (const auto& p) {
        stream << p.first.get() << " " << p.second.probability_true() << '\n';
    });
}

void print_candidate_posteriors(const VariantPosteriorVector& candidate_posteriors,
                                const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    print_candidate_posteriors(std::cout, candidate_posteriors, n);
}

void log(const GenotypeProbabilityMap& genotype_posteriors,
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
         boost::optional<logging::TraceLogger>& trace_log,
         Phred<double> min_posterior)
{
    if (trace_log) {
        print_candidate_posteriors(stream(*trace_log), candidate_posteriors);
    }
    if (debug_log) {
        const auto n = std::count_if(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors),
                                     [=] (const auto& p) { return p.second >= min_posterior; });
        print_candidate_posteriors(stream(*debug_log), candidate_posteriors, std::max(n,  decltype(n) {5}));
    }
}

} // namespace
} // namespace debug

} // namespace octopus
