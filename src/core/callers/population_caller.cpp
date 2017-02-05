// Copyright (c) 2016 Daniel Cooke
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

#include "basics/genomic_region.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "utils/maths.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "containers/probability_matrix.hpp"
#include "core/models/genotype/individual_model.hpp"
#include "core/models/genotype/uniform_population_prior_model.hpp"
#include "core/models/genotype/coalescent_population_prior_model.hpp"
#include "logging/logging.hpp"
#include "utils/germline_variant_call.hpp"
#include "utils/reference_call.hpp"

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
    if (all_equal(parameters_.ploidies)) {
        parameters_.ploidies.resize(1);
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

// IndividualCaller::Latents public methods

namespace {

using InverseGenotypeTable = std::vector<std::vector<std::size_t>>;

auto make_inverse_genotype_table(const std::vector<Haplotype>& haplotypes,
                                 const std::vector<Genotype<Haplotype>>& genotypes)
{
    assert(!haplotypes.empty() && !genotypes.empty());
    using HaplotypeReference = std::reference_wrapper<const Haplotype>;
    std::unordered_map<HaplotypeReference, std::vector<std::size_t>> result_map {haplotypes.size()};
    const auto cardinality = element_cardinality_in_genotypes(static_cast<unsigned>(haplotypes.size()),
                                                              genotypes.front().ploidy());
    for (const auto& haplotype : haplotypes) {
        auto itr = result_map.emplace(std::piecewise_construct,
                                      std::forward_as_tuple(std::cref(haplotype)),
                                      std::forward_as_tuple());
        itr.first->second.reserve(cardinality);
    }
    for (std::size_t i {0}; i < genotypes.size(); ++i) {
        for (const auto& haplotype : genotypes[i]) {
            result_map.at(haplotype).emplace_back(i);
        }
    }
    InverseGenotypeTable result {};
    result.reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) {
        auto& indices = result_map.at(haplotype);
        std::sort(std::begin(indices), std::end(indices));
        indices.erase(std::unique(std::begin(indices), std::end(indices)), std::end(indices));
        result.emplace_back(std::move(indices));
    }
    return result;
}

using GM = model::PopulationModel;

using GenotypeMarginalPosteriorVector = std::vector<double>;
using GenotypeMarginalPosteriorMatrix = std::vector<GenotypeMarginalPosteriorVector>;

auto calculate_genotype_marginal_posteriors(const GM::Latents& posteriors,
                                            const std::size_t num_genotypes,
                                            const std::size_t num_samples)
{
    GenotypeMarginalPosteriorMatrix result {num_samples, GenotypeMarginalPosteriorVector(num_genotypes, 0.0)};
    for (std::size_t i {0}; i < posteriors.genotype_combinations.size(); ++i) {
        for (std::size_t s {0}; s < num_samples; ++s) {
            result[s][posteriors.genotype_combinations[i][s]] += posteriors.joint_genotype_probabilities[i];
        }
    }
    return result;
}

auto calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                                    const std::vector<Genotype<Haplotype>>& genotypes,
                                    const GenotypeMarginalPosteriorMatrix& genotype_posteriors,
                                    const InverseGenotypeTable& inverse_genotypes)
{
    std::unordered_map<std::reference_wrapper<const Haplotype>, double> result {haplotypes.size()};
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
            prob_not_observed *= std::accumulate(std::cbegin(noncontaining_genotype_indices),
                                                 std::cend(noncontaining_genotype_indices),
                                                 0.0, [&sample_genotype_posteriors]
                                                 (const auto curr, const auto i) {
                return curr + sample_genotype_posteriors[i];
            });
        }
        result.emplace(haplotype, 1.0 - prob_not_observed);
        ++itr;
    }
    return result;
}

} // namespace

PopulationCaller::Latents::Latents(const std::vector<SampleName>& samples,
                                   const std::vector<Haplotype>& haplotypes,
                                   std::vector<Genotype<Haplotype>>&& genotypes,
                                   ModelInferences&& inferences)
: model_latents_ {std::move(inferences)}
{
    auto genotype_marginal_posteriors = calculate_genotype_marginal_posteriors(model_latents_.posteriors, genotypes.size(), samples.size());
    auto inverse_genotypes = make_inverse_genotype_table(haplotypes, genotypes);
    haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>(calculate_haplotype_posteriors(haplotypes, genotypes, genotype_marginal_posteriors, inverse_genotypes));
    GenotypeProbabilityMap genotype_posteriors {std::begin(genotypes), std::end(genotypes)};
    for (std::size_t s {0}; s < samples.size(); ++s) {
        insert_sample(samples[s], genotype_marginal_posteriors[s], genotype_posteriors);
    }
    genotype_posteriors_ = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
    genotypes_.emplace(genotypes.front().ploidy(), std::move(genotypes));
}

PopulationCaller::Latents::Latents(const std::vector<SampleName>& samples,
                                   const std::vector<Haplotype>& haplotypes,
                                   std::unordered_map<unsigned, std::vector<Genotype<Haplotype>>>&& genotypes,
                                   ModelInferences&& inferences)
: genotypes_ {std::move(genotypes)}
, model_latents_ {std::move(inferences)}
{
    
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

using GenotypesMap = std::unordered_map<unsigned, std::vector<Genotype<Haplotype>>>;

auto generate_unique_genotypes(const std::vector<Haplotype>& haplotypes, std::vector<unsigned> ploidies)
{
    std::sort(std::begin(ploidies), std::end(ploidies));
    ploidies.erase(std::unique(std::begin(ploidies), std::end(ploidies)), std::end(ploidies));
    GenotypesMap result {ploidies.size()};
    for (auto ploidy : ploidies) {
        result.emplace(ploidy, generate_all_genotypes(haplotypes, ploidy));
    }
    return result;
}

using GenotypeVectorReference = std::reference_wrapper<const std::vector<Genotype<Haplotype>>>;

auto assign_samples_to_genotypes(std::vector<unsigned> ploidies, const GenotypesMap& genotypes)
{
    std::vector<GenotypeVectorReference> result {};
    result.reserve(ploidies.size());
    for (auto ploidy : ploidies) {
        result.emplace_back(genotypes.at(ploidy));
    }
    return result;
}

std::unique_ptr<PopulationCaller::Caller::Latents>
PopulationCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                                const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    const auto prior_model = make_prior_model(haplotypes);
    const model::PopulationModel model {*prior_model, {parameters_.max_genotypes_per_sample}, debug_log_};
    if (parameters_.ploidies.size() == 1) {
        auto genotypes = generate_all_genotypes(haplotypes, parameters_.ploidies.front());
        if (debug_log_) stream(*debug_log_) << "There are " << genotypes.size() << " candidate genotypes";
        auto inferences = model.evaluate(samples_, genotypes, haplotype_likelihoods);
        return std::make_unique<Latents>(samples_, haplotypes, std::move(genotypes), std::move(inferences));
    } else {
        auto unique_genotypes = generate_unique_genotypes(haplotypes, parameters_.ploidies);
        auto sample_genotypes = assign_samples_to_genotypes(parameters_.ploidies, unique_genotypes);
        auto inferences = model.evaluate(samples_, sample_genotypes, haplotype_likelihoods);
        return std::make_unique<Latents>(samples_, haplotypes, std::move(unique_genotypes), std::move(inferences));
    }
}

//auto calculate_model_posterior(const double normal_model_log_evidence,
//                                     const double dummy_model_log_evidence)
//{
//    constexpr double normal_model_prior {0.9999999};
//    constexpr double dummy_model_prior {1.0 - normal_model_prior};
//
//    const auto normal_model_ljp = std::log(normal_model_prior) + normal_model_log_evidence;
//    const auto dummy_model_ljp  = std::log(dummy_model_prior) + dummy_model_log_evidence;
//
//    const auto norm = maths::log_sum_exp(normal_model_ljp, dummy_model_ljp);
//
//    return std::exp(dummy_model_ljp - norm);
//}

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
PopulationCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

namespace {

using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<Haplotype>>::InnerMap;

using VariantReference  = std::reference_wrapper<const Variant>;
using VariantPosteriors = std::vector<std::pair<VariantReference, std::vector<Phred<double>>>>;

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

using AlleleBools          = std::deque<bool>; // using std::deque because std::vector<bool> is evil
using GenotypeContainments = std::vector<AlleleBools>;

auto marginalise(const GenotypeProbabilityMap& genotype_posteriors,
                 const AlleleBools& contained_alleles)
{
    auto p = std::inner_product(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                std::cbegin(contained_alleles), 0.0, std::plus<void> {},
                                [] (const auto& p, const bool is_contained) {
                                    return is_contained ? 0.0 : p.second;
                                });
    return probability_to_phred(p);
}

auto compute_sample_allele_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
                                      const GenotypeContainments& contained_alleles)
{
    std::vector<Phred<double>> result {};
    result.reserve(contained_alleles.size());
    for (const auto& allele : contained_alleles) {
        result.emplace_back(marginalise(genotype_posteriors, allele));
    }
    return result;
}

auto get_contained_alleles(const ProbabilityMatrix<Genotype<Haplotype>>& genotype_posteriors,
                           const std::vector<Allele>& alleles)
{
    const auto num_genotypes = genotype_posteriors.size2();
    GenotypeContainments result {};
    if (num_genotypes == 0 || genotype_posteriors.empty1() || alleles.empty()) {
        return result;
    }
    result.reserve(alleles.size());
    const auto& test_sample   = genotype_posteriors.begin()->first;
    const auto genotype_begin = genotype_posteriors.begin(test_sample);
    const auto genotype_end   = genotype_posteriors.end(test_sample);
    for (const auto& allele : alleles) {
        result.emplace_back(num_genotypes);
        std::transform(genotype_begin, genotype_end, std::begin(result.back()),
                       [&] (const auto& p) {
                           return contains(p.first, allele);
                       });
    }
    return result;
}

auto compute_posteriors(const std::vector<SampleName>& samples,
                        const std::vector<Allele>& alleles,
                        const ProbabilityMatrix<Genotype<Haplotype>>& genotype_posteriors)
{
    const auto contained_alleles = get_contained_alleles(genotype_posteriors, alleles);
    std::vector<std::vector<Phred<double>>> result {};
    result.reserve(genotype_posteriors.size1());
    for (const auto& sample : samples) {
        result.emplace_back(compute_sample_allele_posteriors(genotype_posteriors[sample], contained_alleles));
    }
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
                        const ProbabilityMatrix<Genotype<Haplotype>>& genotype_posteriors)
{
    const auto allele_posteriors = compute_posteriors(samples, extract_alt_alleles(variants), genotype_posteriors);
    VariantPosteriors result {};
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

auto call_genotypes(const GM::Latents& latents, const std::vector<Genotype<Haplotype>>& genotypes)
{
    const auto itr = std::max_element(std::cbegin(latents.joint_genotype_probabilities), std::cend(latents.joint_genotype_probabilities));
    const auto& called_indices = latents.genotype_combinations[std::distance(std::cbegin(latents.joint_genotype_probabilities), itr)];
    std::vector<Genotype<Haplotype>> result {};
    result.reserve(called_indices.size());
    std::transform(std::cbegin(called_indices), std::cend(called_indices), std::back_inserter(result),
                   [&] (auto idx) { return genotypes[idx]; });
    return result;
}

auto call_genotypes(const GM::Latents& latents, const std::unordered_map<unsigned, std::vector<Genotype<Haplotype>>>& genotypes)
{
    return call_genotypes(latents, std::cbegin(genotypes)->second);
}

// variant calling

bool has_above(const std::vector<Phred<double>>& posteriors, const Phred<double> min_posterior)
{
    return std::any_of(std::cbegin(posteriors), std::cend(posteriors), [=] (auto p) { return p >= min_posterior; });
}

bool contains_alt(const Genotype<Haplotype>& genotype_call, const VariantReference& candidate)
{
    return includes(genotype_call, candidate.get().alt_allele());
}

bool contains_alt(const std::vector<Genotype<Haplotype>>& genotype_calls, const VariantReference& candidate)
{
    return std::any_of(std::cbegin(genotype_calls), std::cend(genotype_calls),
                       [&] (const auto& genotype) { return contains_alt(genotype, candidate); });
}

VariantCalls call_candidates(const VariantPosteriors& candidate_posteriors,
                             const std::vector<Genotype<Haplotype>>& genotype_calls,
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

// allele genotype calling

auto marginalise(const Genotype<Allele>& genotype, const GenotypeProbabilityMap& genotype_posteriors)
{
    auto p = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), 0.0,
                             [&genotype] (const double curr, const auto& p) {
                                 return curr + (contains(p.first, genotype) ? 0.0 : p.second);
                             });
    return probability_to_phred(p);
}

auto call_genotypes(const std::vector<SampleName>& samples,
                    const std::vector<Genotype<Haplotype>>& genotype_calls,
                    const ProbabilityMatrix<Genotype<Haplotype>>& genotype_posteriors,
                    const std::vector<GenomicRegion>& variant_regions)
{
    GenotypeCalls result {};
    result.reserve(variant_regions.size());
    for (const auto& region : variant_regions) {
        std::vector<GenotypeCall> region_calls {};
        region_calls.reserve(samples.size());
        for (std::size_t s {0}; s < samples.size(); ++s) {
            auto spliced_genotype = splice<Allele>(genotype_calls[s], region);
            const auto posterior = marginalise(spliced_genotype, genotype_posteriors[samples[s]]);
            region_calls.push_back({std::move(spliced_genotype), posterior});
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
    return std::make_unique<GermlineVariantCall>(variant_call.variant.get(), std::move(tmp),
                                                 max(variant_call.posteriors));
}

auto transform_calls(const std::vector<SampleName>& samples,
                     VariantCalls&& variant_calls,
                     GenotypeCalls&& genotype_calls)
{
    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
    result.reserve(variant_calls.size());
    std::transform(std::make_move_iterator(std::begin(variant_calls)),
                   std::make_move_iterator(std::end(variant_calls)),
                   std::make_move_iterator(std::begin(genotype_calls)),
                   std::back_inserter(result),
                   [&samples] (auto&& variant_call, auto&& genotype_call) {
                       return transform_call(samples, std::move(variant_call), std::move(genotype_call));
                   });
    return result;
}

} // namespace

std::vector<std::unique_ptr<octopus::VariantCall>>
PopulationCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    const auto& genotype_posteriors = (*latents.genotype_posteriors_);
    const auto candidate_posteriors = compute_posteriors(samples_, candidates, genotype_posteriors);
    const auto genotype_calls = call_genotypes(latents.model_latents_.posteriors, latents.genotypes_);
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
                                        const ReadMap& reads) const
{
    return {};
}

std::unique_ptr<PopulationPriorModel> PopulationCaller::make_prior_model(const std::vector<Haplotype>& haplotypes) const
{
    if (parameters_.prior_model_params) {
        return std::make_unique<CoalescentPopulationPriorModel>(CoalescentModel {
        Haplotype {mapped_region(haplotypes.front()), reference_},
        *parameters_.prior_model_params
        });
    } else {
        return std::make_unique<UniformPopulationPriorModel>();
    }
}

//namespace debug {
//
//template <typename S>
//void print_genotype_posteriors(S&& stream,
//                               const GenotypeProbabilityMap& genotype_posteriors,
//                               const std::size_t n)
//{
//    const auto m = std::min(n, genotype_posteriors.size());
//    if (m == genotype_posteriors.size()) {
//        stream << "Printing all genotype posteriors " << '\n';
//    } else {
//        stream << "Printing top " << m << " genotype posteriors " << '\n';
//    }
//    using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
//    std::vector<std::pair<GenotypeReference, double>> v {};
//    v.reserve(genotype_posteriors.size());
//    std::copy(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
//              std::back_inserter(v));
//    const auto mth = std::next(std::begin(v), m);
//    std::partial_sort(std::begin(v), mth, std::end(v),
//                      [] (const auto& lhs, const auto& rhs) {
//                          return lhs.second > rhs.second;
//                      });
//    std::for_each(std::begin(v), mth,
//                  [&] (const auto& p) {
//                      ::debug::print_variant_alleles(stream, p.first);
//                      stream << " " << p.second << '\n';
//                  });
//}
//
//void print_genotype_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
//                               const std::size_t n)
//{
//    print_genotype_posteriors(std::cout, genotype_posteriors, n);
//}
//
//template <typename S>
//void print_candidate_posteriors(S&& stream, const VariantPosteriors& candidate_posteriors,
//                                const std::size_t n)
//{
//    const auto m = std::min(n, candidate_posteriors.size());
//    if (m == candidate_posteriors.size()) {
//        stream << "Printing all candidate variant posteriors " << '\n';
//    } else {
//        stream << "Printing top " << m << " candidate variant posteriors " << '\n';
//    }
//    std::vector<std::pair<VariantReference, double>> v {};
//    v.reserve(candidate_posteriors.size());
//    std::copy(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors),
//              std::back_inserter(v));
//    const auto mth = std::next(std::begin(v), m);
//    std::partial_sort(std::begin(v), mth, std::end(v),
//                      [] (const auto& lhs, const auto& rhs) {
//                          return lhs.second > rhs.second;
//                      });
//    std::for_each(std::begin(v), mth,
//                  [&] (const auto& p) {
//                      stream << p.first.get() << " " << p.second << '\n';
//                  });
//}
//
//void print_candidate_posteriors(const VariantPosteriors& candidate_posteriors,
//                                const std::size_t n)
//{
//    print_candidate_posteriors(std::cout, candidate_posteriors, n);
//}
//
//} // namespace debug

} // namespace octopus
