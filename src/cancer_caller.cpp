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
//    if (debug_log_) {
//        if (has_normal_sample()) {
//            stream(*debug_log_) << "Normal sample is \"" << *parameters_.normal_sample << "\"";
//        } else {
//            *debug_log_ << "There is no normal sample";
//        }
//    }
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
    
    GermlineModel germline_model {ploidy, germline_prior_model};
    CNVModel cnv_model {samples_, ploidy, std::move(cnv_model_priors)};
    SomaticModel somatic_model {samples_, ploidy, std::move(somatic_model_priors)};
    
    const auto merged_likelihoods = merge_samples(samples_, "union", haplotypes, haplotype_likelihoods);
    
    auto germline_inferences = germline_model.infer_latents("union", germline_genotypes,
                                                            merged_likelihoods);
    auto cnv_inferences = cnv_model.infer_latents(germline_genotypes, haplotype_likelihoods);
    
    filter(cancer_genotypes, germline_genotypes, germline_inferences, cnv_inferences);
    
    auto somatic_inferences = somatic_model.infer_latents(cancer_genotypes, haplotype_likelihoods);
    
    if (has_normal_sample()) {
        Latents::NormalInferences normal_inferences;
        
        normal_inferences.germline = germline_model.infer_latents(normal_sample(),
                                                                  germline_genotypes,
                                                                  haplotype_likelihoods);
        
        GenotypeModel::Individual dummy_model {ploidy + 1, germline_prior_model};
        
        auto dummy_genotypes = generate_all_genotypes(haplotypes, ploidy + 1);
        
        normal_inferences.dummy = dummy_model.infer_latents(normal_sample(), dummy_genotypes,
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
    
    auto cmp = [] (const auto& lhs, const auto& rhs) {
        return lhs.get() == rhs.get();
    };
    
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
            Priors::GenotypeMixturesDirichletAlphas sample_alphas {0.5, 0.5};
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
            Priors::GenotypeMixturesDirichletAlphas sample_alphas {5.0, 5.0, 0.01};
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
        constexpr double normal_model_prior {0.9999999};
        constexpr double dummy_model_prior {1.0 - normal_model_prior};
        
        const auto normal_model_ljp = std::log(normal_model_prior) + normal_germline_model_log_evidence;
        const auto dummy_model_ljp  = std::log(dummy_model_prior) + normal_dummy_model_log_evidence;
        
        const auto norm = Maths::log_sum_exp(normal_model_ljp, dummy_model_ljp);
        
        return std::exp(dummy_model_ljp - norm);
    }
}

std::vector<std::unique_ptr<VariantCall>>
CancerVariantCaller::call_variants(const std::vector<Variant>& candidates, Latents& latents) const
{
    if (latents.normal_validation_inferences_) {
        auto germline_evidence = latents.normal_validation_inferences_->germline.log_evidence;
        auto dummy_evidence    = latents.normal_validation_inferences_->dummy.log_evidence;
        
        auto dummy_model_posterior = calculate_dummy_model_posterior(germline_evidence, dummy_evidence);
        
        if (dummy_model_posterior > 0.5) {
            return {};
        }
    }
    
    const auto model_posteriors = calculate_model_posteriors(latents);
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "Germline model posterior: " << model_posteriors.germline;
        stream(log) << "CNV model posterior:      " << model_posteriors.cnv;
        stream(log) << "Somatic model posterior:  " << model_posteriors.somatic;
    }
    
    if (model_posteriors.somatic > model_posteriors.germline) {
        if (model_posteriors.somatic > model_posteriors.cnv) {
            return call_somatic_variants(candidates, latents.somatic_model_inferences_.posteriors);
        } else {
            return call_cnv_variants(candidates, latents.cnv_model_inferences_.posteriors);
        }
    } else if (model_posteriors.cnv > model_posteriors.germline) {
        return call_cnv_variants(candidates, latents.cnv_model_inferences_.posteriors);
    } else {
        return call_germline_variants(candidates, latents.germline_model_inferences_.posteriors);
    }
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

template <typename L>
auto find_map_genotype(const L& posteriors)
{
    return *std::max_element(std::cbegin(posteriors.genotype_probabilities),
                             std::cend(posteriors.genotype_probabilities),
                             [] (const auto& lhs, const auto& rhs) {
                                 return lhs.second < rhs.second;
                             });
}

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

namespace debug
{
    template <typename S>
    void print_map_genotype(S&& stream, const GenotypeModel::Somatic::Latents& posteriors)
    {
        const auto& map = find_map_genotype(posteriors);
        stream << "MAP cancer genotype is:" << '\n';
        debug::print_variant_alleles(stream, map.first);
        stream << " " << map.second << '\n';
    }
    
    void print_map_genotype(const GenotypeModel::Somatic::Latents& posteriors)
    {
        print_map_genotype(std::cout, posteriors);
    }
} // namespace debug

using SomaticModelLatents = GenotypeModel::Somatic::Latents;

bool overlaps(const std::pair<double, double>& lhs, const std::pair<double, double>& rhs)
{
    return std::min(lhs.second, rhs.second) - std::max(lhs.first, rhs.first) > 0;
}

std::vector<std::unique_ptr<VariantCall>>
CancerVariantCaller::call_somatic_variants(const std::vector<Variant>& candidates,
                                           const GenotypeModel::Somatic::Latents& posteriors) const
{
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        debug::print_map_genotype(stream(log), posteriors);
    }
    
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
