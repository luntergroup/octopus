//
//  cancer_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 16/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_caller.hpp"

#include <utility>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <functional>
#include <iostream>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "genomic_region.hpp"
#include "read_pipe.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "merge_transform.hpp"
#include "mappable_algorithms.hpp"
#include "maths.hpp"
#include "cancer_genotype.hpp"
#include "read_utils.hpp"
#include "probability_matrix.hpp"
#include "sequence_utils.hpp"
#include "germline_variant_call.hpp"
#include "reference_call.hpp"
#include "somatic_call.hpp"
#include "logging.hpp"

namespace Octopus
{
// public methods

CancerVariantCaller::CallerParameters::CallerParameters(double min_variant_posterior,
                                                        double min_somatic_posterior,
                                                        double min_refcall_posterior,
                                                        unsigned ploidy,
                                                        boost::optional<SampleIdType> normal_sample,
                                                        double somatic_mutation_rate,
                                                        bool call_somatics_only)
:
min_variant_posterior {min_variant_posterior},
min_somatic_posterior {min_somatic_posterior},
min_refcall_posterior {min_refcall_posterior},
ploidy {ploidy},
normal_sample {normal_sample},
somatic_mutation_rate {somatic_mutation_rate},
call_somatics_only {call_somatics_only},
max_genotypes {50'000}
{}

CancerVariantCaller::CancerVariantCaller(const ReferenceGenome& reference,
                                         ReadPipe& read_pipe,
                                         CandidateVariantGenerator&& candidate_generator,
                                         VariantCaller::CallerParameters general_parameters,
                                         CallerParameters specific_parameters)
:
VariantCaller {reference, read_pipe, std::move(candidate_generator), std::move(general_parameters)},
parameters_ {std::move(specific_parameters)}
{
    if (parameters_.min_variant_posterior == 0) {
        Logging::WarningLogger wlog {};
        wlog << "Having a minumum germline variant posterior means no somatic variants will be called";
    }
    
    if (debug_log_) {
        if (has_normal_sample()) {
            stream(*debug_log_) << "Normal sample is " << *parameters_.normal_sample;
        } else {
            *debug_log_ << "There is no normal sample";
        }
    }
}

CancerVariantCaller::Latents::Latents(std::vector<Genotype<Haplotype>>&& germline_genotypes,
                                      std::vector<CancerGenotype<Haplotype>>&& somatic_genotypes,
                                      GermlineModel::InferredLatents&& germline,
                                      CNVModel::InferredLatents&& cnv,
                                      SomaticModel::InferredLatents&& somatic)
:
germline_genotypes_ {std::move(germline_genotypes)},
somatic_genotypes_ {std::move(somatic_genotypes)},
germline_model_inferences_ {std::move(germline)},
cnv_model_inferences_ {std::move(cnv)},
somatic_model_inferences_ {std::move(somatic)},
normal_validation_inferences_ {}
{}

CancerVariantCaller::Latents::Latents(std::vector<Genotype<Haplotype>>&& germline_genotypes,
                                      std::vector<CancerGenotype<Haplotype>>&& somatic_genotypes,
                                      GermlineModel::InferredLatents&& germline,
                                      CNVModel::InferredLatents&& cnv,
                                      SomaticModel::InferredLatents&& somatic,
                                      NormalInferences&& normal_validation)
:
germline_genotypes_ {std::move(germline_genotypes)},
somatic_genotypes_ {std::move(somatic_genotypes)},
germline_model_inferences_ {std::move(germline)},
cnv_model_inferences_ {std::move(cnv)},
somatic_model_inferences_ {std::move(somatic)},
normal_validation_inferences_ {std::move(normal_validation)}
{}

std::shared_ptr<CancerVariantCaller::Latents::HaplotypeProbabilityMap>
CancerVariantCaller::Latents::get_haplotype_posteriors() const
{
    // TODO: properly
    
    Latents::HaplotypeProbabilityMap result {};
    result.reserve(500);
    
    for (const auto& p : cnv_model_inferences_.posteriors.genotype_probabilities) {
        for (const auto& haplotype : p.first) {
            result[haplotype] += p.second;
        }
    }
    
    for (const auto& p : somatic_model_inferences_.posteriors.genotype_probabilities) {
        for (const auto& haplotype : p.first.get_germline_genotype()) {
            result[haplotype] += p.second;
        }
        result[p.first.get_cancer_element()] += p.second;
    }
    
    const auto norm = Maths::sum_values(result);
    
    for (auto& p : result) {
        p.second /= norm;
    }
    
    return std::make_shared<Latents::HaplotypeProbabilityMap>(std::move(result));
}

std::shared_ptr<CancerVariantCaller::Latents::GenotypeProbabilityMap>
CancerVariantCaller::Latents::get_genotype_posteriors() const
{
    // TODO: properly
    
    GenotypeProbabilityMap genotype_posteriors {
        std::begin(germline_genotypes_),
        std::end(germline_genotypes_)
    };
    
    SampleIdType sample;
    
    // Just to extract sample names
    for (const auto& p : cnv_model_inferences_.posteriors.alphas) {
        std::tie(sample, std::ignore) = p;
        insert_sample(std::move(sample), germline_model_inferences_.posteriors.genotype_probabilities,
                      genotype_posteriors);
    }
    
    return std::make_shared<Latents::GenotypeProbabilityMap>(std::move(genotype_posteriors));
}

// private methods

bool CancerVariantCaller::has_normal_sample() const noexcept
{
    return static_cast<bool>(parameters_.normal_sample);
}

const SampleIdType& CancerVariantCaller::normal_sample() const
{
    return *parameters_.normal_sample;
}

std::unique_ptr<CancerVariantCaller::CallerLatents>
CancerVariantCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                                   const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    const auto ploidy = parameters_.ploidy;
    
    std::vector<CancerGenotype<Haplotype>> cancer_genotypes;
    std::vector<Genotype<Haplotype>> germline_genotypes;
    
    std::tie(cancer_genotypes, germline_genotypes) = generate_all_cancer_genotypes(haplotypes, ploidy);
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "There are " << germline_genotypes.size() << " candidate germline genotypes";
        stream(log) << "There are " << cancer_genotypes.size() << " candidate cancer genotypes";
    }
    
    CoalescentModel germline_prior_model {Haplotype {mapped_region(haplotypes.front()), reference_}};
    SomaticMutationModel somatic_prior_model {germline_prior_model, parameters_.somatic_mutation_rate};
    
    auto cnv_model_priors     = calculate_cnv_model_priors(germline_prior_model);
    auto somatic_model_priors = calculate_somatic_model_priors(somatic_prior_model);
    
    GermlineModel germline_model {germline_prior_model};
    CNVModel cnv_model {samples_, ploidy, std::move(cnv_model_priors)};
    SomaticModel somatic_model {samples_, ploidy, std::move(somatic_model_priors)};
    
    const auto merged_likelihoods = merge_samples(samples_, "union", haplotypes, haplotype_likelihoods);
    
    auto germline_inferences = germline_model.infer_latents("union", germline_genotypes,
                                                            merged_likelihoods);
    auto cnv_inferences = cnv_model.infer_latents(germline_genotypes, haplotype_likelihoods);
    
    filter(cancer_genotypes, germline_genotypes, germline_inferences, cnv_inferences);
    
//    Logging::DebugLogger log {};
//    log << "All cancer genotypes";
//    int i {0};
//    for (const auto& g : cancer_genotypes) {
//        auto s = stream(log); s << i++ << ": ";
//        debug::print_variant_alleles(s, g);
//    }
    
    auto somatic_inferences = somatic_model.infer_latents(cancer_genotypes, haplotype_likelihoods);
    
    if (has_normal_sample()) {
        Latents::NormalInferences normal_inferences;
        
        normal_inferences.germline = germline_model.infer_latents(normal_sample(),
                                                                  germline_genotypes,
                                                                  haplotype_likelihoods);
        
        auto dummy_genotypes = generate_all_genotypes(haplotypes, ploidy + 1);
        
        normal_inferences.dummy = germline_model.infer_latents(normal_sample(), dummy_genotypes,
                                                               haplotype_likelihoods);
        
        return std::make_unique<Latents>(std::move(germline_genotypes), std::move(cancer_genotypes),
                                         std::move(germline_inferences), std::move(cnv_inferences),
                                         std::move(somatic_inferences), std::move(normal_inferences));
    }
    
    return std::make_unique<Latents>(std::move(germline_genotypes), std::move(cancer_genotypes),
                                     std::move(germline_inferences), std::move(cnv_inferences),
                                     std::move(somatic_inferences));
}

template <typename... T>
auto zip(const T&... containers)
    -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
    auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}

auto extract_low_posterior_genotypes(const std::vector<Genotype<Haplotype>>& genotypes,
                                     const GenotypeModel::Individual::Latents& latents,
                                     const double min_posterior = 1e-30)
{
    using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
    
    auto hash = std::hash<GenotypeReference>();
    
    const auto cmp = [] (const auto& lhs, const auto& rhs) { return lhs.get() == rhs.get(); };
    
    std::unordered_set<GenotypeReference, decltype(hash), decltype(cmp)> result {
        genotypes.size(), hash, cmp
    };
    
    for (const auto& p : zip(genotypes, latents.genotype_probabilities)) {
        if (p.get<1>() < min_posterior) {
            result.emplace(p.get<0>());
        }
    }
    
    return result;
}

void CancerVariantCaller::filter(std::vector<CancerGenotype<Haplotype>>& cancer_genotypes,
                                 const std::vector<Genotype<Haplotype>>& germline_genotypes,
                                 const GermlineModel::InferredLatents& germline_inferences,
                                 const CNVModel::InferredLatents& cnv_inferences) const
{
    if (cancer_genotypes.size() <= parameters_.max_genotypes) {
        if (debug_log_) {
            *debug_log_ << "No cancer genotypes were filtered";
        }
        return;
    }
    
    if (has_normal_sample()) {
        const auto removable_germlines = extract_low_posterior_genotypes(germline_genotypes,
                                                                         germline_inferences.posteriors);
        
        const auto it = std::remove_if(std::begin(cancer_genotypes), std::end(cancer_genotypes),
                                       [&removable_germlines] (const auto& g) {
                                           return removable_germlines.count(g.get_germline_genotype()) == 1;
                                       });
        
        cancer_genotypes.erase(it, std::end(cancer_genotypes));
    } else {
        
    }
}

CancerVariantCaller::CNVModel::Priors
CancerVariantCaller::calculate_cnv_model_priors(const CoalescentModel& prior_model) const
{
    using Priors = CNVModel::Priors;
    
    Priors::GenotypeMixturesDirichletAlphaMap cnv_alphas {};
    cnv_alphas.reserve(samples_.size());
    
    for (const auto& sample : samples_) {
        if (has_normal_sample() && sample == normal_sample()) {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas {10.0, 10.0};
            cnv_alphas.emplace(sample, std::move(sample_alphas));
        } else {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas {0.75, 0.75};
            cnv_alphas.emplace(sample, std::move(sample_alphas));
        }
    }
    
    return Priors {prior_model, std::move(cnv_alphas)};
}

CancerVariantCaller::SomaticModel::Priors
CancerVariantCaller::calculate_somatic_model_priors(const SomaticMutationModel& prior_model) const
{
    using Priors = SomaticModel::Priors;
    
    Priors::GenotypeMixturesDirichletAlphaMap alphas {};
    alphas.reserve(samples_.size());
    
    for (const auto& sample : samples_) {
        if (has_normal_sample() && sample == normal_sample()) {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas {10.0, 10.0, 0.01};
            alphas.emplace(sample, std::move(sample_alphas));
        } else {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas {1.0, 1.0, 0.75};
            alphas.emplace(sample, std::move(sample_alphas));
        }
    }
    
    return Priors {prior_model, std::move(alphas)};
}

std::vector<std::unique_ptr<VariantCall>>
CancerVariantCaller::call_variants(const std::vector<Variant>& candidates, CallerLatents& latents) const
{
    return call_variants(candidates, dynamic_cast<Latents&>(latents));
}

namespace
{
    double calculate_dummy_model_posterior(const double normal_germline_model_log_evidence,
                                           const double normal_dummy_model_log_evidence)
    {
        constexpr double normal_model_prior {0.999};
        constexpr double dummy_model_prior {1.0 - normal_model_prior};
        
        const auto normal_model_ljp = std::log(normal_model_prior) + normal_germline_model_log_evidence;
        const auto dummy_model_ljp  = std::log(dummy_model_prior) + normal_dummy_model_log_evidence;
        
        const auto norm = Maths::log_sum_exp(normal_model_ljp, dummy_model_ljp);
        
        return std::exp(dummy_model_ljp - norm);
    }
}

using VariantReference  = std::reference_wrapper<const Variant>;
using VariantPosteriors = std::vector<std::pair<VariantReference, double>>;

struct TempGermlineVariantCall : Mappable<TempGermlineVariantCall>
{
    TempGermlineVariantCall() = delete;
    TempGermlineVariantCall(const std::pair<VariantReference, double>& p)
    : variant {p.first}, posterior {p.second} {}
    TempGermlineVariantCall(const Variant& variant, double posterior)
    : variant {variant}, posterior {posterior} {}
    
    const GenomicRegion& get_region() const noexcept { return mapped_region(variant.get()); }
    
    VariantReference variant;
    double posterior;
};

using GermlineVariantCalls = std::vector<TempGermlineVariantCall>;

struct GenotypeCall
{
    GenotypeCall() = default;
    template <typename T> GenotypeCall(T&& genotype, double posterior)
    : genotype {std::forward<T>(genotype)}, posterior {posterior} {}
    
    Genotype<Allele> genotype;
    double posterior;
};

using GenotypeCalls = std::vector<GenotypeCall>;

    struct CancerGenotypeCall
    {
        CancerGenotypeCall() = default;
        template <typename T>
        CancerGenotypeCall(T&& genotype, double posterior)
        : genotype {std::forward<T>(genotype)}, posterior {posterior} {}
        
        CancerGenotype<Allele> genotype;
        double posterior;
        std::unordered_map<SampleIdType, std::vector<std::pair<double, double>>> credible_regions;
    };
    
    using CancerGenotypeCalls = std::vector<CancerGenotypeCall>;
    
template <typename L>
auto find_map_genotype(const L& posteriors)
{
    return std::max_element(std::cbegin(posteriors),std::cend(posteriors),
                            [] (const auto& lhs, const auto& rhs) {
                                return lhs.second < rhs.second;
                            });
}

//namespace debug
//{
//    template <typename S>
//    void print_map_genotype(S&& stream, const GenotypeModel::Somatic::Latents& posteriors)
//    {
//        const auto& map = find_map_genotype(posteriors);
//        stream << "MAP cancer genotype is:" << '\n';
//        debug::print_variant_alleles(stream, map->first);
//        stream << " " << map->second << '\n';
//    }
//    
//    void print_map_genotype(const GenotypeModel::Somatic::Latents& posteriors)
//    {
//        print_map_genotype(std::cout, posteriors);
//    }
//} // namespace debug

    template <typename T>
    auto compute_marginal_credible_interval(const T& alphas, const double mass)
    {
        const auto a0 = std::accumulate(std::cbegin(alphas), std::cend(alphas), 0.0);
        
        std::vector<std::pair<double, double>> result {};
        result.reserve(alphas.size());
        
        for (const auto& alpha : alphas) {
            result.push_back(Maths::beta_hdi(alpha, a0 - alpha, mass));
        }
        
        return result;
    }
    
    template <typename M>
    auto compute_marginal_credible_intervals(const M& alphas, const double mass)
    {
        std::unordered_map<SampleIdType, std::vector<std::pair<double, double>>> result {};
        result.reserve(alphas.size());
        
        for (const auto& p : alphas) {
            result.emplace(p.first, compute_marginal_credible_interval(p.second, mass));
        }
        
        return result;
    }
    
std::vector<std::unique_ptr<VariantCall>>
CancerVariantCaller::call_variants(const std::vector<Variant>& candidates, Latents& latents) const
{
    if (latents.normal_validation_inferences_) {
        auto germline_evidence = latents.normal_validation_inferences_->germline.log_evidence;
        auto dummy_evidence    = latents.normal_validation_inferences_->dummy.log_evidence;
        
        auto dummy_model_posterior = calculate_dummy_model_posterior(germline_evidence, dummy_evidence);
        
        if (debug_log_) {
            stream(*debug_log_) << "Dummy model posterior = " << dummy_model_posterior;
        }
        
        if (dummy_model_posterior > 0.5) {
            if (debug_log_) {
                *debug_log_ << "Skipping region due to model filter";
            }
            return {};
        }
    }
    
    const auto model_posteriors = calculate_model_posteriors(latents);
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "Germline model posterior: " << model_posteriors.germline;
        stream(log) << "CNV model posterior:      " << model_posteriors.cnv;
        stream(log) << "Somatic model posterior:  " << model_posteriors.somatic;
        
        auto map_cnv = find_map_genotype(latents.cnv_model_inferences_.posteriors.genotype_probabilities);
        auto cnv_log = stream(log);
        cnv_log << "MAP CNV genotype is: ";
        ::debug::print_variant_alleles(cnv_log, map_cnv->first);
        auto somatic_log = stream(log);
        auto map_somatic = find_map_genotype(latents.somatic_model_inferences_.posteriors.genotype_probabilities);
        somatic_log << "MAP cancer genotype is: ";
        debug::print_variant_alleles(somatic_log, map_somatic->first);
    }
    
    const auto& posterior_alphas = latents.somatic_model_inferences_.posteriors.alphas;
    
    std::unordered_map<Genotype<Haplotype>, double> germline_genotype_posteriors {};
    
    std::transform(std::cbegin(latents.germline_genotypes_), std::cend(latents.germline_genotypes_),
                   std::cbegin(latents.germline_model_inferences_.posteriors.genotype_probabilities),
                   std::inserter(germline_genotype_posteriors, std::begin(germline_genotype_posteriors)),
                   [model_posteriors] (const auto& genotype, const auto& posterior) {
                       return std::make_pair(genotype, model_posteriors.germline * posterior);
                   });
    
    for (const auto& p : latents.cnv_model_inferences_.posteriors.genotype_probabilities) {
        germline_genotype_posteriors[p.first] += model_posteriors.cnv * p.second;
    }
    
    for (const auto& p : latents.somatic_model_inferences_.posteriors.genotype_probabilities) {
        germline_genotype_posteriors[p.first.get_germline_genotype()] += model_posteriors.somatic * p.second;
    }
    
//    std::cout << "germline genotype posteriors" << '\n';
//    for (const auto& p : germline_genotype_posteriors) {
//        ::debug::print_variant_alleles(p.first); std::cout << " " << p.second << '\n';
//    }
    
    std::vector<std::pair<VariantReference, double>> germline_candidate_posteriors {};
    germline_candidate_posteriors.reserve(candidates.size());
    
    for (const auto& candidate : candidates) {
        const auto& allele = candidate.get_alt_allele();
        
        auto result = std::accumulate(std::cbegin(germline_genotype_posteriors),
                                      std::cend(germline_genotype_posteriors),
                                      0.0, [&allele] (const auto curr, const auto& p) {
                                          return curr + (contains(p.first, allele) ? p.second : 0);
                                      });
        
        if (result > 1.0) result = 1.0;
        
        germline_candidate_posteriors.emplace_back(candidate, result);
    }
    
    const auto& called_germline_genotype = find_map_genotype(germline_genotype_posteriors)->first;
    
    GermlineVariantCalls germline_variant_calls {};
    germline_variant_calls.reserve(candidates.size());
    
    std::vector<VariantReference> uncalled_germline_candidates {};
    
    for (const auto& p : germline_candidate_posteriors) {
        if (p.second >= parameters_.min_variant_posterior) {
            if (contains_exact(called_germline_genotype, p.first.get().get_alt_allele())) {
                germline_variant_calls.emplace_back(p.first, p.second);
            }
        } else {
            uncalled_germline_candidates.emplace_back(p.first);
        }
    }
    
    std::vector<double> sample_somatic_posteriors(samples_.size());
    
    const auto ploidy = parameters_.ploidy;
    
    std::transform(std::cbegin(posterior_alphas), std::cend(posterior_alphas),
                   std::begin(sample_somatic_posteriors), [ploidy] (const auto& p) {
                       const auto a0 = std::accumulate(std::cbegin(p.second),
                                                       std::prev(std::cend(p.second)), 0.0);
                       return Maths::beta_cdf(p.second.back(), a0, 0.005);
                   });
    
    auto prob_somatic = 1 - std::accumulate(std::cbegin(sample_somatic_posteriors),
                                            std::cend(sample_somatic_posteriors),
                                            1.0, std::multiplies<> {});
    
    prob_somatic *= model_posteriors.somatic;
    
    std::vector<std::unique_ptr<Octopus::VariantCall>> result {};
    
    if (prob_somatic >= parameters_.min_somatic_posterior) {
        std::vector<Allele> candidate_somatic_alleles {};
        
        const auto& cancer_genotype_posteriors = latents.somatic_model_inferences_.posteriors.genotype_probabilities;
        
        std::vector<std::pair<CancerGenotype<Haplotype>, double>> likely_cancer_genotypes {};
        
        std::copy_if(std::cbegin(cancer_genotype_posteriors), std::cend(cancer_genotype_posteriors),
                     std::back_inserter(likely_cancer_genotypes),
                     [] (const auto& p) { return p.second > 0.001; });
        
        GermlineVariantCalls somatic_allele_posteriors {};
        
        for (const auto& candidate : uncalled_germline_candidates) {
            const auto& allele = candidate.get().get_alt_allele();
            
            const auto p = std::accumulate(std::cbegin(likely_cancer_genotypes),
                                           std::cend(likely_cancer_genotypes),
                                           0.0, [&allele] (const auto curr, const auto& p) {
                                               const auto& somatic_haplotype = p.first.get_cancer_element();
                                               
                                               if (!somatic_haplotype.contains(allele)) {
                                                   return curr + 0;
                                               }
                                               
                                               const auto& germline_genotype = p.first.get_germline_genotype();
                                               
                                               if (contains(germline_genotype, allele)) {
                                                   return curr + 0;
                                               }
                                               
                                               return curr + p.second;
                                           });
            
            double somatic_allele_posterior {1};
            
            for (const auto somatic_element_posterior : sample_somatic_posteriors) {
                somatic_allele_posterior *= 1 - (p * somatic_element_posterior);
            }
            
            somatic_allele_posteriors.emplace_back(candidate,  model_posteriors.somatic * (1 - somatic_allele_posterior));
        }
        
        const auto& called_cancer_genotype = find_map_genotype(cancer_genotype_posteriors)->first;
        
        GermlineVariantCalls somatic_variant_calls {};
        somatic_variant_calls.reserve(somatic_allele_posteriors.size());
        
        const auto min_posterior = parameters_.min_somatic_posterior;
        
        std::copy_if(std::begin(somatic_allele_posteriors), std::end(somatic_allele_posteriors),
                     std::back_inserter(somatic_variant_calls),
                     [min_posterior, &called_cancer_genotype] (const auto& p) {
                         return p.posterior >= min_posterior
                         && contains_exact(called_cancer_genotype, p.variant.get().get_alt_allele());
                     });
        
        const auto called_somatic_regions = extract_regions(somatic_variant_calls);
        
        CancerGenotypeCalls cancer_genotype_calls {};
        cancer_genotype_calls.reserve(called_somatic_regions.size());
        
        auto credible_regions = compute_marginal_credible_intervals(latents.somatic_model_inferences_.posteriors.alphas, 0.99);
        
        for (const auto& region : called_somatic_regions) {
            auto spliced_genotype = splice<Allele>(called_cancer_genotype, region);
            
            const auto posterior = std::accumulate(std::cbegin(likely_cancer_genotypes),
                                                   std::cend(likely_cancer_genotypes), 0.0,
                                                   [&spliced_genotype] (const double curr, const auto& p) {
                                                       return curr + ((contains(p.first, spliced_genotype)) ? p.second : 0.0);
                                                   });
            
            cancer_genotype_calls.emplace_back(std::move(spliced_genotype), posterior);
            cancer_genotype_calls.back().credible_regions = credible_regions;
        }
        
        std::transform(std::make_move_iterator(std::begin(somatic_variant_calls)),
                       std::make_move_iterator(std::end(somatic_variant_calls)),
                       std::make_move_iterator(std::begin(cancer_genotype_calls)),
                       std::back_inserter(result),
                       [] (auto&& variant_call, auto&& genotype_call) {
                           std::unordered_map<SampleIdType, SomaticCall::GenotypeCredibleRegions> credible_regions {};
                           
                           for (const auto& p : genotype_call.credible_regions) {
                               SomaticCall::GenotypeCredibleRegions sample_credible_regions {};
                               
                               sample_credible_regions.germline_credible_regions.reserve(p.second.size() - 1);
                               std::copy(std::cbegin(p.second), std::prev(std::cend(p.second)),
                                         std::back_inserter(sample_credible_regions.germline_credible_regions));
                               
                               sample_credible_regions.somatic_credible_region = p.second.back();
                               
                               credible_regions.emplace(p.first, std::move(sample_credible_regions));
                           }
                           
                           return std::make_unique<SomaticCall>(variant_call.variant.get(),
                                                                std::move(genotype_call.genotype),
                                                                genotype_call.posterior,
                                                                std::move(credible_regions),
                                                                variant_call.posterior);
                       });
    }
    
    const auto called_germline_regions = extract_regions(germline_variant_calls);
    
    GenotypeCalls germline_genotype_calls {};
    germline_genotype_calls.reserve(called_germline_regions.size());
    
    for (const auto& region : called_germline_regions) {
        auto spliced_genotype = splice<Allele>(called_germline_genotype, region);
        
        const auto posterior = std::accumulate(std::cbegin(germline_genotype_posteriors),
                                               std::cend(germline_genotype_posteriors), 0.0,
                                               [&called_germline_genotype] (const double curr, const auto& p) {
                                                   return curr + ((contains(p.first, called_germline_genotype)) ? p.second : 0.0);
                                               });
        
        germline_genotype_calls.emplace_back(std::move(spliced_genotype), posterior);
    }
    
    return result;
}

CancerVariantCaller::ModelPosteriors::ModelPosteriors(double germline, double cnv, double somatic)
: germline {germline}, cnv {cnv}, somatic {somatic} {}

CancerVariantCaller::ModelPosteriors
CancerVariantCaller::calculate_model_posteriors(const Latents& inferences) const
{
    const auto& germline_inferences = inferences.germline_model_inferences_;
    const auto& cnv_inferences      = inferences.cnv_model_inferences_;
    const auto& somatic_inferences  = inferences.somatic_model_inferences_;
    
    const double germline_model_prior {0.89999};
    const double cnv_model_prior      {0.1};
    const double somatic_model_prior  {0.01 * 0.001};
    
    const auto germline_model_jlp = std::log(germline_model_prior) + germline_inferences.log_evidence;
    const auto cnv_model_jlp      = std::log(cnv_model_prior) + cnv_inferences.approx_log_evidence;
    const auto somatic_model_jlp  = std::log(somatic_model_prior) + somatic_inferences.approx_log_evidence;
    
    const auto norm = Maths::log_sum_exp(germline_model_jlp, cnv_model_jlp, somatic_model_jlp);
    
    auto germline_model_posterior = std::exp(germline_model_jlp - norm);
    auto cnv_model_posterior      = std::exp(cnv_model_jlp - norm);
    auto somatic_model_posterior  = std::exp(somatic_model_jlp - norm);
    
    return ModelPosteriors {germline_model_posterior, cnv_model_posterior, somatic_model_posterior};
}

std::vector<std::unique_ptr<VariantCall>>
CancerVariantCaller::call_germline_variants(const std::vector<Variant>& candidates,
                                            const GenotypeModel::Individual::Latents& posteriors) const
{
    return {};
}

std::vector<std::unique_ptr<VariantCall>>
CancerVariantCaller::call_cnv_variants(const std::vector<Variant>& candidates,
                                       const GenotypeModel::CNV::Latents& posteriors) const
{
    const auto credible_regions = compute_marginal_credible_intervals(posteriors.alphas, 0.99);
    
    if (debug_log_) {
        auto slog = stream(*debug_log_);
        for (const auto& p : credible_regions) {
            slog << p.first << '\n';
            for (const auto& r : p.second) {
                slog << "    " << r.first << " " << r.second << '\n';
            }
        }
    }
    
    return {};
}

using SomaticModelLatents = GenotypeModel::Somatic::Latents;

bool overlaps(const std::pair<double, double>& lhs, const std::pair<double, double>& rhs)
{
    return std::min(lhs.second, rhs.second) - std::max(lhs.first, rhs.first) > 0;
}

std::vector<std::unique_ptr<VariantCall>>
CancerVariantCaller::call_somatic_variants(const std::vector<Variant>& candidates,
                                           const GenotypeModel::Somatic::Latents& posteriors) const
{
//    if (DEBUG_MODE) {
//        Logging::DebugLogger log {};
//        debug::print_map_genotype(stream(log), posteriors);
//    }
    
    const auto credible_regions = compute_marginal_credible_intervals(posteriors.alphas, 0.99);
    
    if (has_normal_sample()) {
        const auto& normal_somatic_credible_region = credible_regions.at(normal_sample())[parameters_.ploidy];
        
        if (debug_log_) {
            auto slog = stream(*debug_log_);
            
            for (const auto& p : credible_regions) {
                const auto& somatic_credible_region = p.second[parameters_.ploidy];
                
                if (!overlaps(normal_somatic_credible_region, somatic_credible_region)) {
                    slog << "SOMATIC!!!" << '\n';
                }
                
                slog << p.first << '\n';
                for (const auto& r : p.second) {
                    slog << "    " << r.first << " " << r.second << '\n';
                }
            }
        }
    }
    
    
    return {};
}

std::vector<std::unique_ptr<ReferenceCall>>
CancerVariantCaller::call_reference(const std::vector<Allele>& alleles,
                                    CallerLatents& latents,
                                    const ReadMap& reads) const
{
    return {};
}

} // namespace Octopus
