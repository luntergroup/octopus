// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "individual_caller.hpp"

#include <typeinfo>
#include <unordered_map>
#include <deque>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <limits>

#include "basics/genomic_region.hpp"
#include "containers/probability_matrix.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/calls/germline_variant_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "core/models/genotype/uniform_genotype_prior_model.hpp"
#include "core/models/genotype/coalescent_genotype_prior_model.hpp"
#include "utils/maths.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/append.hpp"
#include "utils/select_top_k.hpp"
#include "logging/logging.hpp"

namespace octopus {

IndividualCaller::IndividualCaller(Caller::Components&& components,
                                   Caller::Parameters general_parameters,
                                   Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {std::move(specific_parameters)}
{
    if (parameters_.ploidy == 0) {
        throw std::logic_error {"IndividualCaller: ploidy must be > 0"};
    }
}

std::string IndividualCaller::do_name() const
{
    return "individual";
}

IndividualCaller::CallTypeSet IndividualCaller::do_call_types() const
{
    return {std::type_index(typeid(GermlineVariantCall))};
}

unsigned IndividualCaller::do_min_callable_ploidy() const
{
    return parameters_.ploidy;
}

std::size_t IndividualCaller::do_remove_duplicates(HaplotypeBlock& haplotypes) const
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

IndividualCaller::Latents::Latents(const SampleName& sample,
                                   const IndexedHaplotypeBlock& haplotypes,
                                   IndividualCaller::GenotypeBlock genotypes,
                                   ModelInferences&& inferences)
: genotype_log_posteriors_ {}
, genotype_posteriors_ {}
, haplotype_posteriors_ {}
, haplotypes_ {std::move(haplotypes)}
, genotypes_ {std::move(genotypes)}
, latents_ {std::move(inferences)}
, samples_ {sample}
{}

std::shared_ptr<IndividualCaller::Latents::HaplotypeProbabilityMap>
IndividualCaller::Latents::haplotype_posteriors() const
{
    if (!haplotype_posteriors_) {
        haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>(calculate_haplotype_posteriors(haplotypes_));
    }
    return haplotype_posteriors_;
}

std::shared_ptr<IndividualCaller::Latents::GenotypeProbabilityMap>
IndividualCaller::Latents::genotype_posteriors() const
{
    if (!genotype_posteriors_) {
        GenotypeProbabilityMap genotype_posteriors {std::cbegin(genotypes_), std::cend(genotypes_)};
        insert_sample(samples_, latents_.posteriors.genotype_probabilities, genotype_posteriors);
        genotype_posteriors_ = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
    }
    return genotype_posteriors_;
}

std::shared_ptr<IndividualCaller::Latents::GenotypeProbabilityMap>
IndividualCaller::Latents::genotype_log_posteriors() const
{
    if (!genotype_log_posteriors_) {
        GenotypeProbabilityMap genotype_log_posteriors {std::cbegin(genotypes_), std::cend(genotypes_)};
        insert_sample(samples_, latents_.posteriors.genotype_log_probabilities, genotype_log_posteriors);
        genotype_log_posteriors_ = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_log_posteriors));
    }
    return genotype_log_posteriors_;
}

// IndividualCaller::Latents private methods

IndividualCaller::Latents::HaplotypeProbabilityMap
IndividualCaller::Latents::calculate_haplotype_posteriors(const IndexedHaplotypeBlock& haplotypes) const
{
    const auto& genotype_posteriors = *this->genotype_posteriors();
    HaplotypeProbabilityMap result {haplotypes.size()};
    for (const auto& haplotype : haplotypes) {
        result.emplace(haplotype, 0.0);
    }
    const auto& sample = std::cbegin(genotype_posteriors)->first;
    for (const auto& p : genotype_posteriors[sample]) {
        for (const auto& haplotype : collapse(p.first)) {
            result.at(haplotype) += p.second;
        }
    }
    return result;
}

std::unique_ptr<IndividualCaller::Caller::Latents>
IndividualCaller::infer_latents(const HaplotypeBlock& haplotypes,
                                const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                OptionalThreadPool workers) const
{
    const auto indexed_haplotypes = index(haplotypes);
    auto genotypes = propose_genotypes(haplotypes, indexed_haplotypes, haplotype_likelihoods);
    if (debug_log_) stream(*debug_log_) << "There are " << genotypes.size() << " candidate genotypes";
    auto prior_model = make_prior_model(haplotypes);
    prior_model->prime(haplotypes);
    model::IndividualModel model {*prior_model, debug_log_, trace_log_};
    model.prime(haplotypes);
    haplotype_likelihoods.prime(sample());
    auto inferences = model.evaluate(genotypes, haplotype_likelihoods);
    return std::make_unique<Latents>(sample(), indexed_haplotypes, std::move(genotypes), std::move(inferences));
}

boost::optional<Caller::ModelPosterior>
IndividualCaller::calculate_model_posterior(const HaplotypeBlock& haplotypes,
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
IndividualCaller::calculate_model_posterior(const HaplotypeBlock& haplotypes,
                                            const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                            const Latents& latents) const
{
    const auto indexed_haplotypes = index(haplotypes);
    const auto genotypes = propose_model_check_genotypes(haplotypes, indexed_haplotypes, latents);
    const auto prior_model = make_prior_model(haplotypes);
    prior_model->prime(haplotypes);
    const model::IndividualModel model {*prior_model, debug_log_};
    haplotype_likelihoods.prime(sample());
    const auto inferences = model.evaluate(genotypes, haplotype_likelihoods);
    ModelPosterior result {};
    result.samples = {octopus::calculate_model_posterior(latents.latents_.log_evidence, inferences.log_evidence)};
    return result;
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

using VariantList = std::vector<VariantCall>;

struct GenotypeCall : Mappable<VariantCall>
{
    template <typename T> GenotypeCall(T&& genotype, Phred<double> posterior)
    : genotype {std::forward<T>(genotype)}
    , posterior {posterior}
    {}
    
    const GenomicRegion& mapped_region() const noexcept
    {
        return octopus::mapped_region(genotype);
    }

    Genotype<Allele> genotype;
    Phred<double> posterior;
};

using GenotypeCallVector = std::vector<GenotypeCall>;

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

// variant genotype calling

template <typename PairIterator>
PairIterator find_map(PairIterator first, PairIterator last)
{
    return std::max_element(first, last, [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
}

const auto& get_map(const GenotypeProbabilityMap& genotype_log_posteriors)
{
    return find_map(std::cbegin(genotype_log_posteriors), std::cend(genotype_log_posteriors))->first;
}

GenotypeCallVector
call_genotypes(const Genotype<IndexedHaplotype<>>& genotype_call,
               const GenotypeProbabilityMap& genotype_log_posteriors,
               const std::vector<GenomicRegion>& variant_regions)
{
    GenotypeCallVector result {};
    result.reserve(variant_regions.size());
    for (const auto& region : variant_regions) {
        auto genotype_chunk = copy<Allele>(genotype_call, region);
        const auto posterior = marginalise_contained(genotype_chunk, genotype_log_posteriors);
        result.emplace_back(std::move(genotype_chunk), posterior);
    }
    return result;
}

bool matches_all_genotypes(const Genotype<IndexedHaplotype<>>& genotype, const GenotypeCallVector& genotype_calls)
{
    const auto contains_helper = [&] (const auto& g) { return contains(genotype, g.genotype); };
    return std::all_of(std::cbegin(genotype_calls), std::cend(genotype_calls), contains_helper);
}

using GenotypeProbabilityMapIterator = ProbabilityMatrix<Genotype<IndexedHaplotype<>>>::ZipIterator;
using GenotypeProbabilityIteratorVector = std::vector<GenotypeProbabilityMapIterator>;

auto get_matching_genotype_iterators(const GenotypeProbabilityMap& genotype_log_posteriors,
                                     const GenotypeCallVector& genotype_calls)
{
    GenotypeProbabilityIteratorVector result {};
    result.reserve(genotype_log_posteriors.size() / 2);
    for (auto genotype_itr = std::cbegin(genotype_log_posteriors);
         genotype_itr != std::cend(genotype_log_posteriors);
        ++genotype_itr) {
        if (matches_all_genotypes(genotype_itr->first, genotype_calls)) {
            result.push_back(genotype_itr);
        }
    }
    return result;
}

auto get_map_containing_alt(const GenotypeProbabilityIteratorVector& genotype_itrs, const VariantCall& variant)
{
    boost::optional<GenotypeProbabilityMapIterator> result {};
    for (const auto genotype_itr : genotype_itrs) {
        if ((!result || genotype_itr->second > (*result)->second) && contains_alt(genotype_itr->first, variant.variant)) {
            result = genotype_itr;
        }
    }
    return result;
}

// output

octopus::VariantCall::GenotypeCall convert(GenotypeCall&& call)
{
    return octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
}

std::unique_ptr<octopus::VariantCall>
transform_call(const SampleName& sample, VariantCall&& variant_call, GenotypeCall&& genotype_call)
{
    std::vector<std::pair<SampleName, Call::GenotypeCall>> tmp {std::make_pair(sample, convert(std::move(genotype_call)))};
    std::unique_ptr<octopus::VariantCall> result {std::make_unique<GermlineVariantCall>(variant_call.variant.get(), std::move(tmp),
                                                                                        variant_call.posterior)};
    return result;
}

auto transform_calls(const SampleName& sample, VariantList&& variant_calls,
                     GenotypeCallVector&& genotype_calls)
{
    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
    result.reserve(variant_calls.size());
    std::transform(std::make_move_iterator(std::begin(variant_calls)), std::make_move_iterator(std::end(variant_calls)),
                   std::make_move_iterator(std::begin(genotype_calls)), std::back_inserter(result),
                   [&sample] (VariantCall&& variant_call, GenotypeCall&& genotype_call) {
                       return transform_call(sample, std::move(variant_call), std::move(genotype_call));
                   });
    return result;
}

} // namespace

std::vector<std::unique_ptr<octopus::VariantCall>>
IndividualCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents, OptionalThreadPool workers) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

namespace debug {

void log(const GenotypeProbabilityMap& genotype_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log);
void log(const VariantPosteriorVector& candidate_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log,
         Phred<double> min_posterior);
void log(const Genotype<Haplotype>& called_genotype,
         boost::optional<logging::DebugLogger>& debug_log);

} // namespace debug

std::vector<std::unique_ptr<octopus::VariantCall>>
IndividualCaller::call_variants(const std::vector<Variant>& candidates,
                                const Latents& latents) const
{
    if (parameters_.ploidy == 0) return {};
    const auto& genotype_log_posteriors = (*latents.genotype_log_posteriors())[sample()];
    debug::log(genotype_log_posteriors, debug_log_, trace_log_);
    const auto candidate_posteriors = compute_candidate_posteriors(candidates, genotype_log_posteriors);
    debug::log(candidate_posteriors, debug_log_, trace_log_, parameters_.min_variant_posterior);
    const auto map_genotype = get_map(genotype_log_posteriors);
    VariantList called_variants {}, callable_variants {};
    called_variants.reserve(candidates.size());
    for (auto& candidate: candidate_posteriors) {
        if (candidate.second >= parameters_.min_variant_posterior) {
            if (contains_alt(map_genotype, candidate.first)) {
                called_variants.push_back(std::move(candidate));
            } else {
                callable_variants.push_back(std::move(candidate));
            }
        }
    }
    auto genotype_calls = call_genotypes(map_genotype, genotype_log_posteriors, extract_regions(called_variants));
    // We now try to force call candidates that meet the minimum variant posterior
    // threshold but are not part of the MAP genotype. In order to call these variants,
    // we require that they do not interfere with any of the previously called
    // variant genotypes, otherwise we could end up reporting inviable haplotypes.
    const auto has_interacting = [&] (const auto& candidate) { return has_overlapped(called_variants, candidate); };
    callable_variants.erase(std::remove_if(std::begin(callable_variants), std::end(callable_variants), has_interacting),
                            std::end(callable_variants));
    if (!callable_variants.empty()) {
        auto viable_genotypes = get_matching_genotype_iterators(genotype_log_posteriors, genotype_calls);
        if (!viable_genotypes.empty()) {
            const static auto posterior_greater = [] (const auto& lhs, const auto& rhs) {
                return lhs.posterior > rhs.posterior;
            };
            std::sort(std::begin(callable_variants), std::end(callable_variants), posterior_greater);
            unsigned num_force_calls {0};
            for (auto& candidate : callable_variants) {
                const auto is_interacting = [&] (const auto& variant) { return overlaps(variant, candidate); };
                if (std::none_of(std::crbegin(called_variants), std::next(std::crbegin(called_variants), num_force_calls), is_interacting)) {
                    const auto new_map_genotype_itr = get_map_containing_alt(viable_genotypes, candidate);
                    if (new_map_genotype_itr) {
                        auto genotype_chunk = copy<Allele>((*new_map_genotype_itr)->first, mapped_region(candidate));
                        const auto posterior = marginalise_contained(genotype_chunk, genotype_log_posteriors);
                        const auto not_still_viable = [&] (const auto& genotype_itr) {
                            return !contains(genotype_itr->first, genotype_chunk);
                        };
                        viable_genotypes.erase(std::remove_if(std::begin(viable_genotypes), std::end(viable_genotypes), not_still_viable),
                                            std::end(viable_genotypes));
                        called_variants.push_back(std::move(candidate));
                        genotype_calls.emplace_back(std::move(genotype_chunk), posterior);
                        ++num_force_calls;
                        if (viable_genotypes.empty()) break;
                    }
                }
            }
            callable_variants.clear();
            callable_variants.shrink_to_fit();
            if (num_force_calls > 0) {
                const static auto mapped_region_less = [] (const auto& lhs, const auto& rhs) { return mapped_region(lhs) < mapped_region(rhs); };
                const auto called_variants_itr = std::prev(std::end(called_variants), num_force_calls);
                std::sort(called_variants_itr, std::end(called_variants), mapped_region_less);
                std::inplace_merge(std::begin(called_variants), called_variants_itr, std::end(called_variants));
                const auto genotype_calls_itr = std::prev(std::end(genotype_calls), num_force_calls);
                std::sort(genotype_calls_itr, std::end(genotype_calls), mapped_region_less);
                std::inplace_merge(std::begin(genotype_calls), genotype_calls_itr, std::end(genotype_calls));
            }
        }
    }
    return transform_calls(sample(), std::move(called_variants), std::move(genotype_calls));
}

namespace {

// reference genotype calling

struct RefCall
{
    Allele reference_allele;
    Phred<double> posterior;
};

const GenomicRegion& mapped_region(const GenotypeProbabilityMap& genotype_posteriors)
{
    return mapped_region(std::cbegin(genotype_posteriors)->first);
}

bool contains_helper(const Haplotype& haplotype, const Allele& allele)
{
    if (!is_indel(allele)) {
        return haplotype.contains(allele);
    } else {
        return haplotype.includes(allele);
    }
}

bool has_variation(const Allele& allele, const GenotypeProbabilityMap& genotype_posteriors)
{
    return !genotype_posteriors.empty()
        && contains(mapped_region(genotype_posteriors), allele)
        && !all_of_is_homozygous(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), allele,
                                 [] (const auto& p) { return p.first; });
}

auto marginalise_homozygous(const Allele& allele, const GenotypeProbabilityMap& genotype_log_posteriors)
{
    std::deque<double> buffer {};
    for_each_is_homozygous(std::cbegin(genotype_log_posteriors), std::cend(genotype_log_posteriors), allele,
                           [&] (const auto& p, bool is_hom) { if (!is_hom) buffer.push_back(p.second); },
                           [] (const auto& p) { return p.first; });
    if (!buffer.empty()) {
        return log_probability_false_to_phred(std::min(maths::log_sum_exp(buffer), 0.0));
    } else {
        return Phred<> {std::numeric_limits<double>::infinity()};
    }
}

using ReadPileupRange = ContainedRange<ReadPileups::const_iterator>;

auto compute_homozygous_posterior(const Allele& allele,
                                  const GenotypeProbabilityMap& genotype_posteriors,
                                  const GenotypeProbabilityMap& genotype_log_posteriors,
                                  const ReadPileupRange& pileups)
{
    assert(!empty(pileups));
    if (has_variation(allele, genotype_posteriors)) {
        return marginalise_homozygous(allele, genotype_log_posteriors);
    } else {
        std::vector<AlignedRead::BaseQuality> reference_qualities {}, non_reference_qualities {};
        Allele::NucleotideSequence reference_sequence {};
        std::size_t reference_idx {0};
        for (const ReadPileup& pileup : pileups) {
            reference_sequence.assign(1, allele.sequence()[reference_idx++]);
            pileup.summaries([&] (const auto& sequence, const auto& summaries) {
                if (sequence == reference_sequence) {
                    for (const auto& summary : summaries) {
                        for (const auto base_quality : summary.base_qualities) {
                            reference_qualities.push_back(std::min(base_quality, summary.mapping_quality));
                        }
                    }
                } else if (sequence.size() != reference_sequence.size()) {
                    // indel
                    constexpr ReadPileup::BaseQuality indel_base_quality {30};
                    for (const auto& summary : summaries) {
                        non_reference_qualities.push_back(std::min(indel_base_quality, summary.mapping_quality));
                    }
                } else {
                    // snv/mnv
                    for (const auto& summary : summaries) {
                        for (const auto base_quality : summary.base_qualities) {
                            non_reference_qualities.push_back(std::min(base_quality, summary.mapping_quality));
                        }
                    }
                }
            });
        }
        const auto depth = reference_qualities.size() + non_reference_qualities.size();
        if (depth == 0) return Phred<double> {3.0};
        for (auto& q : reference_qualities) q = std::max(q, AlignedRead::BaseQuality {1});
        for (auto& q : non_reference_qualities) q = std::max(q, AlignedRead::BaseQuality {1});
        std::vector<double> reference_ln_likelihoods(depth), non_reference_ln_likelihoods(depth);
        const auto phred_to_ln = [] (auto phred) { return phred * -maths::constants::ln10Div10<>; };
        const auto phred_to_not_ln = [] (auto phred) { return std::log(1.0 - std::pow(10.0, -phred / 10.0)); };
        auto itr = std::transform(std::cbegin(reference_qualities), std::cend(reference_qualities),
                                  std::begin(reference_ln_likelihoods), phred_to_not_ln);
        std::transform(std::cbegin(non_reference_qualities), std::cend(non_reference_qualities), itr, phred_to_ln);
        itr = std::transform(std::cbegin(reference_qualities), std::cend(reference_qualities),
                             std::begin(non_reference_ln_likelihoods), phred_to_ln);
        std::transform(std::cbegin(non_reference_qualities), std::cend(non_reference_qualities), itr, phred_to_not_ln);
        auto hom_ref_ln_likelihood = std::accumulate(std::cbegin(reference_ln_likelihoods), std::cend(reference_ln_likelihoods), 0.0);
        auto het_alt_ln_likelihood = std::inner_product(std::cbegin(reference_ln_likelihoods), std::cend(reference_ln_likelihoods),
                                                        std::cbegin(non_reference_ln_likelihoods), 0.0, std::plus<> {},
                                                        [] (auto ref, auto alt) { return maths::log_sum_exp(ref, alt) - std::log(2); });
        const auto het_ln_posterior = het_alt_ln_likelihood - maths::log_sum_exp(hom_ref_ln_likelihood, het_alt_ln_likelihood);
        return log_probability_false_to_phred(het_ln_posterior);
    }
}

auto call_reference(const std::vector<Allele>& reference_alleles,
                    const GenotypeProbabilityMap& genotype_posteriors,
                    const GenotypeProbabilityMap& genotype_log_posteriors,
                    const ReadPileups& pileups)
{
    assert(std::is_sorted(std::cbegin(reference_alleles), std::cend(reference_alleles)));
    std::vector<RefCall> result {};
    result.reserve(reference_alleles.size());
    auto active_pileup_itr = std::cbegin(pileups);
    for (const auto& allele : reference_alleles) {
        const auto active_pileups = contained_range(active_pileup_itr, std::cend(pileups), contig_region(allele));
        const auto posterior = compute_homozygous_posterior(allele, genotype_posteriors, genotype_log_posteriors, active_pileups);
        result.push_back({allele, posterior});
        active_pileup_itr = active_pileups.end().base();
    }
    return result;
}

auto transform_calls(std::vector<RefCall>&& calls, const SampleName& sample, const unsigned ploidy)
{
    std::vector<std::unique_ptr<ReferenceCall>> result {};
    result.reserve(calls.size());
    std::transform(std::make_move_iterator(std::begin(calls)), std::make_move_iterator(std::end(calls)),
                   std::back_inserter(result),
                   [&] (auto&& call) {
                       std::map<SampleName, ReferenceCall::GenotypeCall> genotype {{sample, {ploidy, call.posterior}}};
                       return std::make_unique<ReferenceCall>(std::move(call.reference_allele), call.posterior, std::move(genotype));
                   });
    return result;
}

} // namespace

std::vector<std::unique_ptr<ReferenceCall>>
IndividualCaller::call_reference(const std::vector<Allele>& alleles,
                                 const Caller::Latents& latents,
                                 const ReadPileupMap& pileups) const
{
    return call_reference(alleles, dynamic_cast<const Latents&>(latents), pileups);
}

std::vector<std::unique_ptr<ReferenceCall>>
IndividualCaller::call_reference(const std::vector<Allele>& alleles,
                                 const Latents& latents,
                                 const ReadPileupMap& pileups) const
{
    const auto& genotype_posteriors = (*latents.genotype_posteriors())[sample()];
    const auto& genotype_log_posteriors = (*latents.genotype_log_posteriors())[sample()];
    auto calls = octopus::call_reference(alleles, genotype_posteriors, genotype_log_posteriors,
                                         pileups.at(sample()));
    return transform_calls(std::move(calls), sample(), parameters_.ploidy);
}

const SampleName& IndividualCaller::sample() const noexcept
{
    return samples_.front();
}

std::unique_ptr<GenotypePriorModel> IndividualCaller::make_prior_model(const HaplotypeBlock& haplotypes) const
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
template <typename Container>
void erase_complement_indices(Container& values, const std::vector<std::size_t>& indices)
{
    std::vector<bool> keep(values.size(), true);
    for (auto idx : indices) keep[idx] = false;
    if (!indices.empty()) {
        std::size_t idx {0};
        const auto is_not_index = [&] (const auto& v) { return keep[idx++]; };
        values.erase(std::remove_if(std::begin(values), std::end(values), is_not_index), std::end(values));
    } else {
        values.clear();
    }
}

template <typename Container, typename Range>
void select_top_k(Container& items, const Range& probabilities, const std::size_t k)
{
    const auto best_indices = select_top_k_indices(probabilities, k, false);
    erase_complement_indices(items, best_indices);
}

template <typename Container, typename T, typename Compare = std::less<T>>
void 
sort_by_other(Container& items, const std::vector<T>& other,
              Compare&& comp = Compare {})
{
    assert(items.size() == other.size());
    using U = typename Container::value_type;
    std::vector<std::pair<U, std::reference_wrapper<const T>>> zipped {};
    zipped.reserve(items.size());
    for (std::size_t i {0}; i < items.size(); ++i) {
        zipped.emplace_back(std::move(items[i]), other[i]);
    }
    const auto comp_helper = [&] (const auto& lhs, const auto& rhs) {
        return comp(lhs.second, rhs.second);
    };
    std::sort(std::begin(zipped), std::end(zipped), comp_helper);
    for (std::size_t i {0}; i < items.size(); ++i) {
        items[i] = std::move(zipped[i].first);
    }
}

template <typename IndexType>
void erase_duplicates(MappableBlock<Genotype<IndexedHaplotype<IndexType>>>& genotypes)
{
    using std::begin; using std::end;
    std::sort(begin(genotypes), end(genotypes), GenotypeLess {});
    genotypes.erase(std::unique(begin(genotypes), end(genotypes)), end(genotypes));
}

} // namespace

IndividualCaller::GenotypeBlock
IndividualCaller::propose_genotypes(const HaplotypeBlock& haplotypes,
                                    const IndexedHaplotypeBlock& indexed_haplotypes,
                                    const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto num_possible_genotypes = num_genotypes_noexcept(haplotypes.size(), parameters_.ploidy);
    GenotypeBlock result {mapped_region(haplotypes)};
    if (!parameters_.max_genotypes || (num_possible_genotypes && *num_possible_genotypes <= *parameters_.max_genotypes)) {
        result = generate_all_genotypes(indexed_haplotypes, parameters_.ploidy);
    } else {
        if (debug_log_) {
            if (num_possible_genotypes) {
                stream(*debug_log_) << "Applying genotype reduction as there are " << *num_possible_genotypes << " possible genotypes";
            } else {
                *debug_log_ << "Applying genotype reduction as the number of possible genotypes calculation overflowed";
            }
        }
        auto ploidy = parameters_.ploidy - 1;
        for (; ploidy > 1; --ploidy) {
            const auto num_ploidy_genotypes = num_genotypes_noexcept(haplotypes.size(), ploidy);
            if (num_ploidy_genotypes && *num_ploidy_genotypes <= *parameters_.max_genotypes) break;
        }
        if (debug_log_) stream(*debug_log_) << "Starting genotype reduction with ploidy " << ploidy;
        result = generate_all_genotypes(indexed_haplotypes, ploidy);
        auto prior_model = make_prior_model(haplotypes);
        prior_model->prime(haplotypes);
        model::IndividualModel model {*prior_model};
        model.prime(haplotypes);
        haplotype_likelihoods.prime(sample());
        for (; ploidy < parameters_.ploidy; ++ploidy) {
            if (debug_log_) stream(*debug_log_) << "Finding good genotypes with ploidy " << ploidy << " from " << result.size();
            const auto ploidy_inferences = model.evaluate(result, haplotype_likelihoods);
            sort_by_other(result, ploidy_inferences.posteriors.genotype_log_probabilities);
            GenotypeBlock next_result {};
            next_result.reserve(*parameters_.max_genotypes);
            // We dont know the right number of seed genotypes since there can be duplicates after expansion with a new haplotype.
            for (int n {0}; n < 3 && next_result.size() < *parameters_.max_genotypes; ++n) {
                const std::size_t num_seeds {std::max((*parameters_.max_genotypes - next_result.size()) / haplotypes.size(), std::size_t {1})};
                auto seed_itr = std::prev(std::end(result), num_seeds);
                const static auto is_hom_ref = [] (const auto& genotype) { return is_homozygous_reference(genotype); };
                if (std::find_if(seed_itr, std::end(result), is_hom_ref) == std::end(result)) {
                    // Ensure reference genotype is always included, helping to keep QUAL in reasonable range
                    const auto hom_ref_itr = std::find_if(std::begin(result), seed_itr, is_hom_ref);
                    std::iter_swap(hom_ref_itr, seed_itr);
                }
                extend(seed_itr, std::end(result), indexed_haplotypes, std::back_inserter(next_result));
                erase_duplicates(next_result);
                result.erase(seed_itr, std::end(result));
                result.shrink_to_fit();
            }
            result = std::move(next_result);
        }
    }
    return result;
}

IndividualCaller::GenotypeBlock 
IndividualCaller::propose_model_check_genotypes(const HaplotypeBlock& haplotypes,
                                                const IndexedHaplotypeBlock& indexed_haplotypes,
                                                const Latents& latents) const
{
    auto current_genotypes = latents.genotypes_;
    sort_by_other(current_genotypes, latents.latents_.posteriors.genotype_log_probabilities);
    const auto num_seeds = std::min(std::size_t {3}, current_genotypes.size());
    auto seed_itr = std::prev(std::end(current_genotypes), num_seeds);
    GenotypeBlock result {mapped_region(haplotypes)};
    extend(seed_itr, std::end(current_genotypes), indexed_haplotypes, std::back_inserter(result));
    erase_duplicates(result);
    return result;
}

namespace debug {

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
    std::for_each(std::begin(v), mth,
                  [&] (const auto& p) {
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
    std::for_each(std::begin(v), mth,
                  [&] (const auto& p) {
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

} // namespace debug
} // namespace octopus
