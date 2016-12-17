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

PopulationCaller::PopulationCaller(Caller::Components&& components,
                                   Caller::Parameters general_parameters,
                                   Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {specific_parameters}
{}

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

using GenotypeMarginalPosteriorVector  = std::vector<double>;
using GenotypeMarginalPosteriorMatrix  = std::vector<GenotypeMarginalPosteriorVector>;

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
{
    auto genotype_marginal_posteriors = calculate_genotype_marginal_posteriors(inferences.posteriors, genotypes.size(), samples.size());
    auto inverse_genotypes = make_inverse_genotype_table(haplotypes, genotypes);
    haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>(calculate_haplotype_posteriors(haplotypes, genotypes, genotype_marginal_posteriors, inverse_genotypes));
    GenotypeProbabilityMap genotype_posteriors {
        std::make_move_iterator(std::begin(genotypes)),
        std::make_move_iterator(std::end(genotypes))
    };
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

std::unique_ptr<PopulationCaller::Caller::Latents>
PopulationCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                                const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    auto genotypes = generate_all_genotypes(haplotypes, parameters_.ploidy);
    if (debug_log_) stream(*debug_log_) << "There are " << genotypes.size() << " candidate genotypes";
    const auto prior_model = make_prior_model(haplotypes);
    const model::PopulationModel model {*prior_model, debug_log_};
    auto inferences = model.evaluate(samples_, genotypes, haplotype_likelihoods);
    return std::make_unique<Latents>(samples_, haplotypes, std::move(genotypes), std::move(inferences));
}

namespace {

using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<Haplotype>>::InnerMap;

using VariantReference  = std::reference_wrapper<const Variant>;
using VariantPosteriors = std::vector<std::pair<VariantReference, double>>;

struct VariantCall : Mappable<VariantCall>
{
    VariantCall() = delete;
    VariantCall(const std::pair<VariantReference, double>& p)
    : variant {p.first}
    , posterior {p.second}
    {}
    VariantCall(const Variant& variant, double posterior)
    : variant {variant}
    , posterior {posterior}
    {}
    
    const GenomicRegion& mapped_region() const noexcept
    {
        return octopus::mapped_region(variant.get());
    }
    
    VariantReference variant;
    double posterior;
};

using VariantCalls = std::vector<VariantCall>;

struct GenotypeCall
{
    GenotypeCall() = default;
    template <typename T> GenotypeCall(T&& genotype, double posterior)
    : genotype {std::forward<T>(genotype)}
    , posterior {posterior}
    {}
    
    Genotype<Allele> genotype;
    double posterior;
};

using GenotypeCalls = std::vector<GenotypeCall>;

} // namespace

//namespace debug
//{
//    template <typename S>
//    void print_genotype_posteriors(S&& stream, const GenotypeProbabilityMap& genotype_posteriors,
//                                   std::size_t n = 5);
//    void print_genotype_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
//                                   std::size_t n = 5);
//    template <typename S>
//    void print_candidate_posteriors(S&& stream, const VariantPosteriors& candidate_posteriors,
//                                    std::size_t n = 10);
//    void print_candidate_posteriors(const VariantPosteriors& candidate_posteriors,
//                                    std::size_t n = 10);
//    //        void print_variant_calls(const VariantCallBlocks& calls);
//} // namespace debug

namespace  {

//// allele posterior calculations
//
//auto marginalise(const Allele& allele, const GenotypeProbabilityMap& genotype_posteriors)
//{
//    auto result = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
//                                  0.0, [&allele] (const auto curr, const auto& p) {
//                                      return curr + (contains(p.first, allele) ? p.second : 0);
//                                  });
//    if (result > 1.0) result = 1.0; // just to account for floating point error
//    return result;
//}
//
//VariantPosteriors compute_candidate_posteriors(const std::vector<Variant>& candidates,
//                                               const GenotypeProbabilityMap& genotype_posteriors)
//{
//    VariantPosteriors result {};
//    result.reserve(candidates.size());
//    for (const auto& candidate : candidates) {
//        result.emplace_back(candidate, marginalise(candidate.alt_allele(), genotype_posteriors));
//    }
//    return result;
//}
//
//// variant calling
//
//bool contains_alt(const Genotype<Haplotype>& genotype_call, const VariantReference& candidate)
//{
//    return contains_exact(genotype_call, candidate.get().alt_allele());
//}
//
//VariantCalls call_candidates(const VariantPosteriors& candidate_posteriors,
//                             const Genotype<Haplotype>& genotype_call,
//                             const double min_posterior)
//{
//    VariantCalls result {};
//    result.reserve(candidate_posteriors.size());
//    std::copy_if(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors),
//                 std::back_inserter(result),
//                 [&genotype_call, min_posterior] (const auto& p) {
//                     return p.second >= min_posterior && contains_alt(genotype_call, p.first);
//                 });
//    return result;
//}
//
//// variant genotype calling
//
//auto call_genotype(const GenotypeProbabilityMap& genotype_posteriors)
//{
//    return std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
//                            [] (const auto& lhs, const auto& rhs) {
//                                return lhs.second < rhs.second;
//                            })->first;
//}
//
//double marginalise(const Genotype<Allele>& genotype, const GenotypeProbabilityMap& genotype_posteriors)
//{
//    return std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), 0.0,
//                           [&genotype] (const double curr, const auto& p) {
//                               return curr + ((contains(p.first, genotype)) ? p.second : 0.0);
//                           });
//}
//
//GenotypeCalls call_genotypes(const Genotype<Haplotype>& genotype_call,
//                             const GenotypeProbabilityMap& genotype_posteriors,
//                             const std::vector<GenomicRegion>& variant_regions)
//{
//    GenotypeCalls result {};
//    result.reserve(variant_regions.size());
//    for (const auto& region : variant_regions) {
//        auto spliced_genotype = splice<Allele>(genotype_call, region);
//        const auto posterior = marginalise(spliced_genotype, genotype_posteriors);
//        result.emplace_back(std::move(spliced_genotype), posterior);
//    }
//    return result;
//}
//
//// output
//
//octopus::VariantCall::GenotypeCall convert(GenotypeCall&& call)
//{
//    return octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
//}
//
//std::unique_ptr<octopus::VariantCall>
//transform_call(const SampleName& sample, VariantCall&& variant_call, GenotypeCall&& genotype_call)
//{
//    std::vector<std::pair<SampleName, Call::GenotypeCall>> tmp {
//        std::make_pair(sample, convert(std::move(genotype_call)))
//    };
//    return std::make_unique<GermlineVariantCall>(variant_call.variant.get(),
//                                                 std::move(tmp), variant_call.posterior);
//}
//
//auto transform_calls(const SampleName& sample, VariantCalls&& variant_calls,
//                     GenotypeCalls&& genotype_calls)
//{
//    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
//    result.reserve(variant_calls.size());
//    std::transform(std::make_move_iterator(std::begin(variant_calls)),
//                   std::make_move_iterator(std::end(variant_calls)),
//                   std::make_move_iterator(std::begin(genotype_calls)),
//                   std::back_inserter(result),
//                   [&sample] (VariantCall&& variant_call, GenotypeCall&& genotype_call) {
//                       return transform_call(sample, std::move(variant_call), std::move(genotype_call));
//                   });
//    return result;
//}

} // namespace

std::vector<std::unique_ptr<octopus::VariantCall>>
PopulationCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
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

std::vector<std::unique_ptr<octopus::VariantCall>>
PopulationCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
//    const auto& genotype_posteriors = (*latents.genotype_posteriors_)[sample_];
//    const auto candidate_posteriors = compute_candidate_posteriors(candidates, genotype_posteriors);
//    const auto genotype_call = call_genotype(genotype_posteriors);
//    auto variant_calls = call_candidates(candidate_posteriors, genotype_call, min_variant_posterior_);
//    const auto called_regions = extract_regions(variant_calls);
//    auto genotype_calls = call_genotypes(genotype_call, genotype_posteriors, called_regions);
//    return transform_calls(sample_, std::move(variant_calls), std::move(genotype_calls));
    return {};
}

namespace  {

// reference genotype calling

struct RefCall : public Mappable<RefCall>
{
    RefCall() = default;
    
    template <typename A>
    RefCall(A&& reference_allele, double posterior)
    :
    reference_allele {std::forward<A>(reference_allele)},
    posterior {posterior}
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

//namespace debug
//{
//    template <typename S>
//    void print_genotype_posteriors(S&& stream,
//                                   const GenotypeProbabilityMap& genotype_posteriors,
//                                   const std::size_t n)
//    {
//        const auto m = std::min(n, genotype_posteriors.size());
//        
//        if (m == genotype_posteriors.size()) {
//            stream << "Printing all genotype posteriors " << '\n';
//        } else {
//            stream << "Printing top " << m << " genotype posteriors " << '\n';
//        }
//        
//        using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
//        
//        std::vector<std::pair<GenotypeReference, double>> v {};
//        v.reserve(genotype_posteriors.size());
//        
//        std::copy(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
//                  std::back_inserter(v));
//        
//        const auto mth = std::next(std::begin(v), m);
//        
//        std::partial_sort(std::begin(v), mth, std::end(v),
//                          [] (const auto& lhs, const auto& rhs) {
//                              return lhs.second > rhs.second;
//                          });
//        
//        std::for_each(std::begin(v), mth,
//                      [&] (const auto& p) {
//                          ::debug::print_variant_alleles(stream, p.first);
//                          stream << " " << p.second << '\n';
//                      });
//    }
//    
//    void print_genotype_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
//                                   const std::size_t n)
//    {
//        print_genotype_posteriors(std::cout, genotype_posteriors, n);
//    }
//    
//    template <typename S>
//    void print_candidate_posteriors(S&& stream, const VariantPosteriors& candidate_posteriors,
//                                    const std::size_t n)
//    {
//        const auto m = std::min(n, candidate_posteriors.size());
//        
//        if (m == candidate_posteriors.size()) {
//            stream << "Printing all candidate variant posteriors " << '\n';
//        } else {
//            stream << "Printing top " << m << " candidate variant posteriors " << '\n';
//        }
//        
//        std::vector<std::pair<VariantReference, double>> v {};
//        v.reserve(candidate_posteriors.size());
//        
//        std::copy(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors),
//                  std::back_inserter(v));
//        
//        const auto mth = std::next(std::begin(v), m);
//        
//        std::partial_sort(std::begin(v), mth, std::end(v),
//                          [] (const auto& lhs, const auto& rhs) {
//                              return lhs.second > rhs.second;
//                          });
//        
//        std::for_each(std::begin(v), mth,
//                      [&] (const auto& p) {
//                          stream << p.first.get() << " " << p.second << '\n';
//                      });
//    }
//    
//    void print_candidate_posteriors(const VariantPosteriors& candidate_posteriors,
//                                    const std::size_t n)
//    {
//        print_candidate_posteriors(std::cout, candidate_posteriors, n);
//    }
//} // namespace debug
} // namespace octopus
