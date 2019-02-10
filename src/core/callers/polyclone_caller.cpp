// Copyright (c) 2015-2019 Daniel Cooke
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

#include "basics/genomic_region.hpp"
#include "containers/probability_matrix.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/calls/germline_variant_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "core/models/genotype/uniform_genotype_prior_model.hpp"
#include "core/models/genotype/coalescent_genotype_prior_model.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/concat.hpp"
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
}

std::string PolycloneCaller::do_name() const
{
    return "polyclone";
}

PolycloneCaller::CallTypeSet PolycloneCaller::do_call_types() const
{
    return {std::type_index(typeid(GermlineVariantCall))};
}

unsigned PolycloneCaller::do_min_callable_ploidy() const
{
    return 1;
}

unsigned PolycloneCaller::do_max_callable_ploidy() const
{
    return parameters_.max_clones;
}

std::size_t PolycloneCaller::do_remove_duplicates(std::vector<Haplotype>& haplotypes) const
{
    if (parameters_.deduplicate_haplotypes_with_germline_model) {
        if (haplotypes.size() < 2) return 0;
        CoalescentModel::Parameters model_params {};
        if (parameters_.prior_model_params) model_params = *parameters_.prior_model_params;
        Haplotype reference {mapped_region(haplotypes.front()), reference_.get()};
        CoalescentModel model {std::move(reference), model_params, haplotypes.size(), CoalescentModel::CachingStrategy::none};
        const CoalescentProbabilityGreater cmp {std::move(model)};
        return octopus::remove_duplicates(haplotypes, cmp);
    } else {
        return Caller::do_remove_duplicates(haplotypes);
    }
}

// PolycloneCaller::Latents public methods

PolycloneCaller::Latents::Latents(std::vector<Genotype<Haplotype>> haploid_genotypes, std::vector<Genotype<Haplotype>> polyploid_genotypes,
                                   HaploidModelInferences haploid_model_inferences, SubloneModelInferences subclone_model_inferences,
                                   const SampleName& sample, const std::function<double(unsigned)>& clonality_prior)
: haploid_genotypes_ {std::move(haploid_genotypes)}
, polyploid_genotypes_ {std::move(polyploid_genotypes)}
, haploid_model_inferences_ {std::move(haploid_model_inferences)}
, subclone_model_inferences_ {std::move(subclone_model_inferences)}
, model_posteriors_ {}
, sample_ {sample}
{
    if (!polyploid_genotypes_.empty()) {
        const auto haploid_model_prior = std::log(clonality_prior(1));
        const auto called_subclonality = polyploid_genotypes_.front().ploidy();
        const auto subclone_model_prior = std::log(clonality_prior(called_subclonality));
        const auto haploid_model_jp = haploid_model_prior + haploid_model_inferences_.log_evidence;
        const auto subclone_model_jp = subclone_model_prior + subclone_model_inferences_.approx_log_evidence;
        const auto norm = maths::log_sum_exp({haploid_model_jp, subclone_model_jp});
        model_posteriors_.clonal = std::exp(haploid_model_jp - norm);
        model_posteriors_.subclonal = std::exp(subclone_model_jp - norm);
    } else {
        model_posteriors_.clonal = 1.0;
        model_posteriors_.subclonal = 0.0;
    }
}

std::shared_ptr<PolycloneCaller::Latents::HaplotypeProbabilityMap>
PolycloneCaller::Latents::haplotype_posteriors() const noexcept
{
    if (haplotype_posteriors_ == nullptr) {
        haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>();
        for (const auto& p : (*(this->genotype_posteriors()))[sample_]) {
            for (const auto& haplotype : p.first.copy_unique_ref()) {
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
        const auto genotypes = concat(haploid_genotypes_, polyploid_genotypes_);
        auto posteriors = concat(haploid_model_inferences_.posteriors.genotype_probabilities,
                                 subclone_model_inferences_.posteriors.genotype_probabilities);
        std::for_each(std::begin(posteriors), std::next(std::begin(posteriors), haploid_genotypes_.size()),
                      [this] (auto& p) { p *= model_posteriors_.clonal; });
        std::for_each(std::next(std::begin(posteriors), haploid_genotypes_.size()), std::end(posteriors),
                      [this] (auto& p) { p *= model_posteriors_.subclonal; });
        genotype_posteriors_ = std::make_shared<GenotypeProbabilityMap>(std::make_move_iterator(std::begin(genotypes)),
                                                                        std::make_move_iterator(std::end(genotypes)));
        insert_sample(sample_, posteriors, *genotype_posteriors_);
    }
    return genotype_posteriors_;
}

// PolycloneCaller::Latents private methods

namespace {

auto make_sublone_model_mixture_prior_map(const SampleName& sample, const unsigned num_clones, const double alpha = 0.5)
{
    model::SubcloneModel::Priors::GenotypeMixturesDirichletAlphaMap result {};
    model::SubcloneModel::Priors::GenotypeMixturesDirichletAlphas alphas(num_clones, alpha);
    result.emplace(sample, std::move(alphas));
    return result;
}

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

void reduce(std::vector<Genotype<Haplotype>>& genotypes, const GenotypePriorModel& genotype_prior_model,
            const HaplotypeLikelihoodArray& haplotype_likelihoods, const std::size_t n)
{
    if (genotypes.size() <= n) return;
    const model::IndividualModel approx_model {genotype_prior_model};
    const auto approx_posteriors = approx_model.evaluate(genotypes, haplotype_likelihoods).posteriors.genotype_probabilities;
    const auto min_posterior = nth_greatest_value(approx_posteriors, n + 1);
    std::size_t idx {0};
    genotypes.erase(std::remove_if(std::begin(genotypes), std::end(genotypes),
                                   [&] (const auto& genotype) { return approx_posteriors[idx++] <= min_posterior; }),
                    std::end(genotypes));
}

void fit_sublone_model(const std::vector<Haplotype>& haplotypes, const HaplotypeLikelihoodArray& haplotype_likelihoods,
                       const GenotypePriorModel& genotype_prior_model, const SampleName& sample, const unsigned max_clones,
                       const double haploid_model_evidence, const std::function<double(unsigned)>& clonality_prior,
                       const std::size_t max_genotypes, std::vector<Genotype<Haplotype>>& polyploid_genotypes,
                       model::SubcloneModel::InferredLatents& sublonal_inferences,
                       boost::optional<logging::DebugLogger>& debug_log)
{
    const auto haploid_prior = std::log(clonality_prior(1));
    for (unsigned num_clones {2}; num_clones <= max_clones; ++num_clones) {
        const auto clonal_model_prior = clonality_prior(num_clones);
        if (clonal_model_prior == 0.0) break;
        auto genotypes = generate_all_max_zygosity_genotypes(haplotypes, num_clones);
        reduce(genotypes, genotype_prior_model, haplotype_likelihoods, max_genotypes);
        if (debug_log) stream(*debug_log) << "Generated " << genotypes.size() << " genotypes with clonality " << num_clones;
        if (genotypes.empty()) break;
        model::SubcloneModel::Priors subclonal_model_priors {genotype_prior_model, make_sublone_model_mixture_prior_map(sample, num_clones)};
        model::SubcloneModel subclonal_model {{sample}, subclonal_model_priors};
        auto inferences = subclonal_model.evaluate(genotypes, haplotype_likelihoods);
        if (debug_log) stream(*debug_log) << "Evidence for model with clonality " << num_clones << " is " << inferences.approx_log_evidence;
        if (num_clones == 2) {
            polyploid_genotypes = std::move(genotypes);
            sublonal_inferences = std::move(inferences);
            if ((std::log(clonal_model_prior) + sublonal_inferences.approx_log_evidence)
                < (haploid_prior + haploid_model_evidence)) break;
        } else {
            if ((std::log(clonal_model_prior) + inferences.approx_log_evidence)
                <= (std::log(clonality_prior(num_clones - 1)) + sublonal_inferences.approx_log_evidence))  break;
            polyploid_genotypes = std::move(genotypes);
            sublonal_inferences = std::move(inferences);
        }
    }
}

} // namespace

std::unique_ptr<PolycloneCaller::Caller::Latents>
PolycloneCaller::infer_latents(const std::vector<Haplotype>& haplotypes, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    auto haploid_genotypes = generate_all_genotypes(haplotypes, 1);
    if (debug_log_) stream(*debug_log_) << "There are " << haploid_genotypes.size() << " candidate haploid genotypes";
    const auto genotype_prior_model = make_prior_model(haplotypes);
    const model::IndividualModel haploid_model {*genotype_prior_model, debug_log_};
    haplotype_likelihoods.prime(sample());
    auto haploid_inferences = haploid_model.evaluate(haploid_genotypes, haplotype_likelihoods);
    if (debug_log_) stream(*debug_log_) << "Evidence for haploid model is " << haploid_inferences.log_evidence;
    std::vector<Genotype<Haplotype>> polyploid_genotypes; model::SubcloneModel::InferredLatents sublonal_inferences;
    fit_sublone_model(haplotypes, haplotype_likelihoods, *genotype_prior_model, sample(), parameters_.max_clones,
                      haploid_inferences.log_evidence, parameters_.clonality_prior, parameters_.max_genotypes, polyploid_genotypes,
                      sublonal_inferences, debug_log_);
    if (debug_log_) stream(*debug_log_) << "There are " << polyploid_genotypes.size() << " candidate polyploid genotypes";
    using std::move;
    return std::make_unique<Latents>(move(haploid_genotypes), move(polyploid_genotypes),
                                     move(haploid_inferences), move(sublonal_inferences),
                                     sample(), parameters_.clonality_prior);
}

boost::optional<double>
PolycloneCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                           const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                           const Caller::Latents& latents) const
{
    return calculate_model_posterior(haplotypes, haplotype_likelihoods, dynamic_cast<const Latents&>(latents));
}

boost::optional<double>
PolycloneCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                           const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                           const Latents& latents) const
{
    return boost::none;
}

std::vector<std::unique_ptr<octopus::VariantCall>>
PolycloneCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

namespace {

using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<Haplotype>>::InnerMap;
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

auto compute_posterior(const Allele& allele, const GenotypeProbabilityMap& genotype_posteriors)
{
    auto p = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                             0.0, [&allele] (const auto curr, const auto& p) {
        return curr + (contains(p.first, allele) ? 0.0 : p.second);
    });
    return probability_to_phred(p);
}

auto compute_candidate_posteriors(const std::vector<Variant>& candidates,
                                  const GenotypeProbabilityMap& genotype_posteriors)
{
    VariantPosteriorVector result {};
    result.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        result.emplace_back(candidate, compute_posterior(candidate.alt_allele(), genotype_posteriors));
    }
    return result;
}

// variant calling

bool has_callable(const VariantPosteriorVector& variant_posteriors, const Phred<double> min_posterior) noexcept
{
    return std::any_of(std::cbegin(variant_posteriors), std::cend(variant_posteriors),
                       [=] (const auto& p) noexcept { return p.second >= min_posterior; });
}

bool contains_alt(const Genotype<Haplotype>& genotype_call, const VariantReference& candidate)
{
    return includes(genotype_call, candidate.get().alt_allele());
}

VariantCalls call_candidates(const VariantPosteriorVector& candidate_posteriors,
                             const Genotype<Haplotype>& genotype_call,
                             const Phred<double> min_posterior)
{
    VariantCalls result {};
    result.reserve(candidate_posteriors.size());
    std::copy_if(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors), std::back_inserter(result),
                 [&genotype_call, min_posterior] (const auto& p) {
                     return p.second >= min_posterior && contains_alt(genotype_call, p.first);
                 });
    return result;
}

// variant genotype calling

template <typename PairIterator>
PairIterator find_map(PairIterator first, PairIterator last)
{
    return std::max_element(first, last, [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
}

template <typename T>
bool is_homozygous_reference(const Genotype<T>& g)
{
    return is_reference(g[0]) && g.is_homozygous();
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
                             [&genotype] (const double curr, const auto& p) {
                                 return curr + (contains(p.first, genotype) ? 0.0 : p.second);
                             });
    return probability_to_phred(p);
}

GenotypeCalls call_genotypes(const Genotype<Haplotype>& genotype_call,
                             const GenotypeProbabilityMap& genotype_posteriors,
                             const std::vector<GenomicRegion>& variant_regions)
{
    GenotypeCalls result {};
    result.reserve(variant_regions.size());
    for (const auto& region : variant_regions) {
        auto genotype_chunk = copy<Allele>(genotype_call, region);
        const auto posterior = compute_posterior(genotype_chunk, genotype_posteriors);
        result.emplace_back(std::move(genotype_chunk), posterior);
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

auto transform_calls(const SampleName& sample, VariantCalls&& variant_calls, GenotypeCalls&& genotype_calls)
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

namespace debug { namespace {

void log(const GenotypeProbabilityMap& genotype_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log);

void log(const VariantPosteriorVector& candidate_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log,
         Phred<double> min_posterior);

void log(const Genotype<Haplotype>& called_genotype,
         boost::optional<logging::DebugLogger>& debug_log);

} // namespace
} // namespace debug

std::vector<std::unique_ptr<octopus::VariantCall>>
PolycloneCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    log(latents);
    const auto& genotype_posteriors = (*latents.genotype_posteriors())[sample()];
    debug::log(genotype_posteriors, debug_log_, trace_log_);
    const auto candidate_posteriors = compute_candidate_posteriors(candidates, genotype_posteriors);
    debug::log(candidate_posteriors, debug_log_, trace_log_, parameters_.min_variant_posterior);
    const bool force_call_non_ref {has_callable(candidate_posteriors, parameters_.min_variant_posterior)};
    const auto genotype_call = octopus::call_genotype(genotype_posteriors, force_call_non_ref);
    auto variant_calls = call_candidates(candidate_posteriors, genotype_call, parameters_.min_variant_posterior);
    const auto called_regions = extract_regions(variant_calls);
    auto genotype_calls = call_genotypes(genotype_call, genotype_posteriors, called_regions);
    return transform_calls(sample(), std::move(variant_calls), std::move(genotype_calls));
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

std::unique_ptr<GenotypePriorModel> PolycloneCaller::make_prior_model(const std::vector<Haplotype>& haplotypes) const
{
    if (parameters_.prior_model_params) {
        return std::make_unique<CoalescentGenotypePriorModel>(CoalescentModel {
        Haplotype {mapped_region(haplotypes.front()), reference_},
        *parameters_.prior_model_params, haplotypes.size(), CoalescentModel::CachingStrategy::address
        });
    } else {
        return std::make_unique<UniformGenotypePriorModel>();
    }
}

void PolycloneCaller::log(const Latents& latents) const
{
    if (debug_log_) {
        stream(*debug_log_) << "Clonal model posterior is " << latents.model_posteriors_.clonal
                            << " and subclonal model posterior is " << latents.model_posteriors_.subclonal;
        if (latents.model_posteriors_.subclonal > latents.model_posteriors_.clonal) {
            stream(*debug_log_) << "Detected subclonality is " << latents.polyploid_genotypes_.front().ploidy();
        }
    }
}

namespace debug { namespace {

template <typename S>
void print_genotype_posteriors(S&& stream,
                               const GenotypeProbabilityMap& genotype_posteriors,
                               const std::size_t n)
{
    const auto m = std::min(n, genotype_posteriors.size());
    if (m == genotype_posteriors.size()) {
        stream << "Printing all genotype posteriors " << '\n';
    } else {
        stream << "Printing top " << m << " genotype posteriors " << '\n';
    }
    using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
    std::vector<std::pair<GenotypeReference, double>> v {};
    v.reserve(genotype_posteriors.size());
    std::copy(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), std::back_inserter(v));
    const auto mth = std::next(std::begin(v), m);
    std::partial_sort(std::begin(v), mth, std::end(v),
                      [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
    std::for_each(std::begin(v), mth,
                  [&] (const auto& p) {
                      print_variant_alleles(stream, p.first);
                      stream << " " << p.second << '\n';
                  });
}

void print_genotype_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
                               const std::size_t n)
{
    print_genotype_posteriors(std::cout, genotype_posteriors, n);
}

template <typename S>
void print_candidate_posteriors(S&& stream, const VariantPosteriorVector& candidate_posteriors,
                                const std::size_t n)
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
                                const std::size_t n)
{
    print_candidate_posteriors(std::cout, candidate_posteriors, n);
}

void log(const GenotypeProbabilityMap& genotype_posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log)
{
    if (trace_log) {
        print_genotype_posteriors(stream(*trace_log), genotype_posteriors, -1);
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
        print_candidate_posteriors(stream(*trace_log), candidate_posteriors, -1);
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
