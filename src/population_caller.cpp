//
//  population_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_caller.hpp"

#include <unordered_map>
#include <deque>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <utility>
#include <iostream>

#include "genomic_region.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "maths.hpp"
#include "mappable_algorithms.hpp"
#include "read_utils.hpp"
#include "probability_matrix.hpp"
#include "coalescent_model.hpp"
#include "individual_genotype_model.hpp"
#include "germline_variant_call.hpp"
#include "reference_call.hpp"
#include "logging.hpp"

namespace Octopus
{
PopulationVariantCaller::CallerParameters::CallerParameters(double min_variant_posterior,
                                                            double min_refcall_posterior,
                                                            unsigned ploidy)
:
min_variant_posterior {min_variant_posterior},
min_refcall_posterior {min_refcall_posterior},
ploidy {ploidy}
{}

PopulationVariantCaller::PopulationVariantCaller(const ReferenceGenome& reference,
                                                 ReadPipe& read_pipe,
                                                 CandidateVariantGenerator&& candidate_generator,
                                                 VariantCaller::CallerParameters general_parameters,
                                                 CallerParameters specific_parameters)
:
VariantCaller {reference, read_pipe, std::move(candidate_generator), std::move(general_parameters)},
ploidy_ {specific_parameters.ploidy},
min_variant_posterior_ {specific_parameters.min_variant_posterior},
min_refcall_posterior_ {specific_parameters.min_refcall_posterior}
{}

// IndividualVariantCaller::Latents public methods

PopulationVariantCaller::Latents::Latents(const std::vector<SampleIdType>& samples,
                                          const std::vector<Haplotype>& haplotypes,
                                          std::vector<Genotype<Haplotype>>&& genotypes,
                                          ModelInferences&& inferences)
:
genotype_posteriors_ {},
haplotype_posteriors_ {},
model_log_evidence_ {inferences.log_evidence},
dummy_latents_ {}
{
//    GenotypeProbabilityMap genotype_posteriors {
//        std::make_move_iterator(std::begin(genotypes)),
//        std::make_move_iterator(std::end(genotypes))
//    };
//    
//    insert_sample(sample, inferences.posteriors.genotype_probabilities, genotype_posteriors);
//    
//    genotype_posteriors_  = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
//    haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>(calculate_haplotype_posteriors(haplotypes));
}

PopulationVariantCaller::Latents::Latents(const std::vector<SampleIdType>& samples,
                                          const std::vector<Haplotype>& haplotypes,
                                          std::vector<Genotype<Haplotype>>&& genotypes,
                                          ModelInferences&& inferences,
                                          ModelInferences&& dummy_inferences)
:
genotype_posteriors_ {},
haplotype_posteriors_ {},
model_log_evidence_ {inferences.log_evidence},
dummy_latents_ {std::move(dummy_inferences)}
{
//    GenotypeProbabilityMap genotype_posteriors {
//        std::make_move_iterator(std::begin(genotypes)),
//        std::make_move_iterator(std::end(genotypes))
//    };
//    
//    insert_sample(sample, inferences.posteriors.genotype_probabilities, genotype_posteriors);
//    
//    genotype_posteriors_  = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
//    haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>(calculate_haplotype_posteriors(haplotypes));
}

std::shared_ptr<PopulationVariantCaller::Latents::HaplotypeProbabilityMap>
PopulationVariantCaller::Latents::get_haplotype_posteriors() const noexcept
{
    return haplotype_posteriors_;
}

std::shared_ptr<PopulationVariantCaller::Latents::GenotypeProbabilityMap>
PopulationVariantCaller::Latents::get_genotype_posteriors() const noexcept
{
    return genotype_posteriors_;
}

// PopulationVariantCaller::Latents private methods

PopulationVariantCaller::Latents::HaplotypeProbabilityMap
PopulationVariantCaller::Latents::calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes)
{
    assert(genotype_posteriors_ != nullptr);
    
    HaplotypeProbabilityMap result {haplotypes.size()};
    
    for (const auto& haplotype : haplotypes) {
        result.emplace(haplotype, 0.0);
    }
    
    const auto& sample = std::cbegin(*genotype_posteriors_)->first;
    
    for (const auto& p : (*genotype_posteriors_)[sample]) {
        for (const auto& haplotype : p.first.copy_unique_ref()) {
            result.at(haplotype) += p.second;
        }
    }
    
    return result;
}

std::unique_ptr<PopulationVariantCaller::CallerLatents>
PopulationVariantCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                                       const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    CoalescentModel prior_model {Haplotype {mapped_region(haplotypes.front()), reference_}};
    
    GenotypeModel::Population model {prior_model};
    
    auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "There are " << genotypes.size() << " candidate genotypes";
    }
    
    auto inferences = model.infer_latents(samples_, genotypes, haplotype_likelihoods);
    
//    // TEST
//    GenotypeModel::Population dummy_model {ploidy_ + 1, prior_model};
//    auto dummy_genotypes = generate_all_genotypes(haplotypes, ploidy_ + 1);
//    if (debug_log_) {
//        stream(*debug_log_) << "Evaluating dummy model with " << dummy_genotypes.size() << " genotypes";
//    }
//    auto dummy_inferences = dummy_model.infer_latents(sample, dummy_genotypes, haplotype_likelihoods);
//    return std::make_unique<Latents>(sample, haplotypes, std::move(genotypes), std::move(inferences),
//                                     std::move(dummy_inferences));
//    // END TEST
    
    return std::make_unique<Latents>(samples_, haplotypes, std::move(genotypes), std::move(inferences));
}

namespace
{
    using GM = GenotypeModel::Population;
    
    using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<Haplotype>>::InnerMap;
    
    using VariantReference  = std::reference_wrapper<const Variant>;
    using VariantPosteriors = std::vector<std::pair<VariantReference, double>>;
    
    struct VariantCall : Mappable<VariantCall>
    {
        VariantCall() = delete;
        VariantCall(const std::pair<VariantReference, double>& p)
        : variant {p.first}, posterior {p.second} {}
        VariantCall(const Variant& variant, double posterior)
        : variant {variant}, posterior {posterior} {}
        
        const GenomicRegion& get_region() const noexcept { return mapped_region(variant.get()); }
        
        VariantReference variant;
        double posterior;
    };
    
    using VariantCalls = std::vector<VariantCall>;
    
    struct GenotypeCall
    {
        GenotypeCall() = default;
        template <typename T> GenotypeCall(T&& genotype, double posterior)
        : genotype {std::forward<T>(genotype)}, posterior {posterior} {}
        
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

namespace
{
    // allele posterior calculations
    
    auto marginalise(const Allele& allele, const GenotypeProbabilityMap& genotype_posteriors)
    {
        auto result = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                      0.0, [&allele] (const auto curr, const auto& p) {
                                          return curr + (contains(p.first, allele) ? p.second : 0);
                                      });
        
        if (result > 1.0) result = 1.0; // just to account for floating point error
        
        return result;
    }
    
    VariantPosteriors compute_candidate_posteriors(const std::vector<Variant>& candidates,
                                                   const GenotypeProbabilityMap& genotype_posteriors)
    {
        VariantPosteriors result {};
        result.reserve(candidates.size());
        
        for (const auto& candidate : candidates) {
            result.emplace_back(candidate, marginalise(candidate.get_alt_allele(), genotype_posteriors));
        }
        
        return result;
    }
    
    // variant calling
    
    bool contains_alt(const Genotype<Haplotype>& genotype_call, const VariantReference& candidate)
    {
        return contains_exact(genotype_call, candidate.get().get_alt_allele());
    }
    
    VariantCalls call_candidates(const VariantPosteriors& candidate_posteriors,
                                 const Genotype<Haplotype>& genotype_call,
                                 const double min_posterior)
    {
        VariantCalls result {};
        result.reserve(candidate_posteriors.size());
        
        std::copy_if(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors),
                     std::back_inserter(result),
                     [&genotype_call, min_posterior] (const auto& p) {
                         return p.second >= min_posterior && contains_alt(genotype_call, p.first);
                     });
        
        result.shrink_to_fit();
        
        return result;
    }
    
    // variant genotype calling
    
    auto call_genotype(const GenotypeProbabilityMap& genotype_posteriors)
    {
        return std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                [] (const auto& lhs, const auto& rhs) {
                                    return lhs.second < rhs.second;
                                })->first;
    }
    
    double marginalise(const Genotype<Allele>& genotype, const GenotypeProbabilityMap& genotype_posteriors)
    {
        return std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), 0.0,
                               [&genotype] (const double curr, const auto& p) {
                                   return curr + ((contains(p.first, genotype)) ? p.second : 0.0);
                               });
    }
    
    GenotypeCalls call_genotypes(const Genotype<Haplotype>& genotype_call,
                                 const GenotypeProbabilityMap& genotype_posteriors,
                                 const std::vector<GenomicRegion>& variant_regions)
    {
        GenotypeCalls result {};
        result.reserve(variant_regions.size());
        
        for (const auto& region : variant_regions) {
            auto spliced_genotype = splice<Allele>(genotype_call, region);
            
            const auto posterior = marginalise(spliced_genotype, genotype_posteriors);
            
            result.emplace_back(std::move(spliced_genotype), posterior);
        }
        
        return result;
    }
    
    // output
    
    Octopus::VariantCall::GenotypeCall convert(GenotypeCall&& call)
    {
        return Octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
    }
    
    std::unique_ptr<Octopus::VariantCall>
    transform_call(const SampleIdType& sample, VariantCall&& variant_call, GenotypeCall&& genotype_call)
    {
        std::vector<std::pair<SampleIdType, Call::GenotypeCall>> tmp {
            std::make_pair(sample, convert(std::move(genotype_call)))
        };
        return std::make_unique<GermlineVariantCall>(variant_call.variant.get(),
                                                     std::move(tmp), variant_call.posterior);
    }
    
    auto transform_calls(const SampleIdType& sample, VariantCalls&& variant_calls,
                         GenotypeCalls&& genotype_calls)
    {
        std::vector<std::unique_ptr<Octopus::VariantCall>> result {};
        result.reserve(variant_calls.size());
        
        std::transform(std::make_move_iterator(std::begin(variant_calls)),
                       std::make_move_iterator(std::end(variant_calls)),
                       std::make_move_iterator(std::begin(genotype_calls)),
                       std::back_inserter(result),
                       [&sample] (VariantCall&& variant_call, GenotypeCall&& genotype_call) {
                           return transform_call(sample, std::move(variant_call), std::move(genotype_call));
                       });
        
        return result;
    }
} // namespace

std::vector<std::unique_ptr<Octopus::VariantCall>>
PopulationVariantCaller::call_variants(const std::vector<Variant>& candidates,
                                       CallerLatents& latents) const
{
    return call_variants(candidates, dynamic_cast<Latents&>(latents));
}

//auto calculate_dummy_model_posterior(const double normal_model_log_evidence,
//                                     const double dummy_model_log_evidence)
//{
//    constexpr double normal_model_prior {0.9999999};
//    constexpr double dummy_model_prior {1.0 - normal_model_prior};
//    
//    const auto normal_model_ljp = std::log(normal_model_prior) + normal_model_log_evidence;
//    const auto dummy_model_ljp  = std::log(dummy_model_prior) + dummy_model_log_evidence;
//    
//    const auto norm = Maths::log_sum_exp(normal_model_ljp, dummy_model_ljp);
//    
//    return std::exp(dummy_model_ljp - norm);
//}

std::vector<std::unique_ptr<Octopus::VariantCall>>
PopulationVariantCaller::call_variants(const std::vector<Variant>& candidates,
                                       const Latents& latents) const
{
//    const auto& genotype_posteriors = (*latents.genotype_posteriors_)[sample_];
//    
//    if (latents.dummy_latents_) {
//        const auto dummy_model_posterior = calculate_dummy_model_posterior(latents.model_log_evidence_,
//                                                                           latents.dummy_latents_->log_evidence);
//        
//        if (debug_log_) {
//            stream(*debug_log_) << "Dummy model posterior = " << dummy_model_posterior;
//        }
//        
//        if (dummy_model_posterior > 0.5) {
//            if (debug_log_) {
//                *debug_log_ << "Skipping region due to model filter";
//            }
//            return {};
//        }
//    }
//    
//    if (TRACE_MODE) {
//        Logging::TraceLogger log {};
//        debug::print_genotype_posteriors(stream(log), genotype_posteriors, -1);
//    } else if (debug_log_) {
//        debug::print_genotype_posteriors(stream(*debug_log_), genotype_posteriors);
//    }
//    
//    const auto candidate_posteriors = compute_candidate_posteriors(candidates, genotype_posteriors);
//    
//    if (TRACE_MODE) {
//        Logging::TraceLogger log {};
//        debug::print_candidate_posteriors(stream(log), candidate_posteriors, -1);
//    } else if (debug_log_) {
//        debug::print_candidate_posteriors(stream(*debug_log_), candidate_posteriors);
//    }
//    
//    const auto genotype_call = call_genotype(genotype_posteriors);
//    
//    auto variant_calls = call_candidates(candidate_posteriors, genotype_call, min_variant_posterior_);
//    
//    const auto called_regions = extract_regions(variant_calls);
//    
//    auto genotype_calls = call_genotypes(genotype_call, genotype_posteriors, called_regions);
//    
//    return transform_calls(sample_, std::move(variant_calls), std::move(genotype_calls));
    return {};
}

namespace
{
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
        
        const GenomicRegion& get_region() const noexcept { return reference_allele.get_region(); }
        
        Allele reference_allele;
        double posterior;
    };
    
    using RefCalls = std::vector<RefCall>;
    
    double marginalise_reference_genotype(const Allele& reference_allele,
                                          const GenotypeProbabilityMap& sample_genotype_posteriors)
    {
        double result {0};
        
        for (const auto& genotype_posterior : sample_genotype_posteriors) {
            if (is_homozygous(genotype_posterior.first, reference_allele)) {
                result += genotype_posterior.second;
            }
        }
        
        return result;
    }
    
    RefCalls call_reference(const GenotypeProbabilityMap& genotype_posteriors,
                            const std::vector<Allele>& reference_alleles,
                            const ReadMap::mapped_type& reads, const double min_call_posterior)
    {
        RefCalls result {};
        
        if (reference_alleles.empty()) return result;
        
        result.reserve(reference_alleles.size());
        
        for (const auto& reference_allele : reference_alleles) {
            double posterior {0};
            
            if (has_coverage(reads, mapped_region(reference_allele))) {
                posterior = marginalise_reference_genotype(reference_allele,
                                                           genotype_posteriors);
            }
            
            if (posterior >= min_call_posterior) {
                result.emplace_back(reference_allele, posterior);
            }
        }
        
        result.shrink_to_fit();
        
        return result;
    }
} // namespace

std::vector<std::unique_ptr<ReferenceCall>>
PopulationVariantCaller::call_reference(const std::vector<Allele>& alleles,
                                        CallerLatents& latents,
                                        const ReadMap& reads) const
{
    return {};
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
} // namespace Octopus
