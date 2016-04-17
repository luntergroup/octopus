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
#include <iostream>

#include "genomic_region.hpp"
#include "read_pipe.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "merge_transform.hpp"
#include "vcf_record.hpp"
#include "mappable_algorithms.hpp"
#include "maths.hpp"
#include "cancer_genotype.hpp"
#include "read_utils.hpp"
#include "string_utils.hpp"
#include "probability_matrix.hpp"
#include "sequence_utils.hpp"

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
max_genotypes {500'000}
{}

CancerVariantCaller::CancerVariantCaller(const ReferenceGenome& reference,
                                         ReadPipe& read_pipe,
                                         CandidateVariantGenerator&& candidate_generator,
                                         VariantCaller::CallerParameters general_parameters,
                                         CallerParameters specific_parameters)
:
VariantCaller {reference, read_pipe, std::move(candidate_generator), std::move(general_parameters)},
parameters_ {std::move(specific_parameters)}
{}

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
somatic_model_inferences_ {std::move(somatic)}
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
    
    filter(cancer_genotypes);
    
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
    auto somatic_inferences = somatic_model.infer_latents(cancer_genotypes, haplotype_likelihoods);
    
    return std::make_unique<Latents>(std::move(germline_genotypes), std::move(cancer_genotypes),
                                     std::move(germline_inferences), std::move(cnv_inferences),
                                     std::move(somatic_inferences));
}

void CancerVariantCaller::filter(std::vector<CancerGenotype<Haplotype>>& genotypes) const
{
    
}

CancerVariantCaller::CNVModel::Priors
CancerVariantCaller::calculate_cnv_model_priors(const CoalescentModel& prior_model) const
{
    using Priors = CNVModel::Priors;
    
    Priors::GenotypeMixturesDirichletAlphaMap cnv_alphas {};
    cnv_alphas.reserve(samples_.size());
    
    for (const auto& sample : samples_) {
        if (parameters_.normal_sample && sample == *parameters_.normal_sample) {
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
        if (parameters_.normal_sample && sample == *parameters_.normal_sample) {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas {5.0, 5.0, 0.01};
            alphas.emplace(sample, std::move(sample_alphas));
        } else {
            Priors::GenotypeMixturesDirichletAlphas sample_alphas {1.0, 1.0, 0.75};
            alphas.emplace(sample, std::move(sample_alphas));
        }
    }
    
    return Priors {prior_model, std::move(alphas)};
}

namespace
{
using GM = GenotypeModel::Somatic;

struct VariantCall
{
    VariantCall() = default;
    template <typename T>
    VariantCall(T&& variants, double posterior)
    : variants {std::forward<T>(variants)}, posterior {posterior} {}
    
    std::vector<Variant> variants;
    double posterior;
};

using VariantCalls = std::vector<VariantCall>;

struct GermlineGenotypeCall
{
    GermlineGenotypeCall() = default;
    template <typename T>
    GermlineGenotypeCall(T&& genotype, double posterior) : genotype {std::forward<T>(genotype)}, posterior {posterior} {}
    
    Genotype<Allele> genotype;
    double posterior;
};

using GermlineGenotypeCalls = std::vector<GermlineGenotypeCall>;

struct SomaticCall
{
    SomaticCall() = default;
    template <typename T>
    SomaticCall(T&& allele, double posterior) : allele {std::forward<T>(allele)}, posterior {posterior} {}
    
    Allele allele;
    double posterior;
};

using SomaticCalls = std::vector<SomaticCall>;

struct RefCall
{
    RefCall() = default;
    template <typename A, typename T>
    RefCall(A&& reference_allele, double posterior, T&& sample_posteriors)
    :
    reference_allele {std::forward<A>(reference_allele)},
    posterior {posterior},
    sample_posteriors {std::forward<T>(sample_posteriors)}
    {}
    
    Allele reference_allele;
    double posterior;
    std::vector<std::pair<SampleIdType, double>> sample_posteriors;
};

using RefCalls = std::vector<RefCall>;

static std::vector<VcfRecord::SequenceType> to_vcf_genotype(const Genotype<Allele>& genotype)
{
    std::vector<VcfRecord::SequenceType> result {};
    result.reserve(genotype.ploidy());
    for (const auto& allele : genotype) result.push_back(allele.get_sequence());
    return result;
}

VcfRecord::Builder output_germline_variant_call(const Allele& reference_allele,
                                       const std::vector<Variant>& variants,
                                       const GermlineGenotypeCall& genotype_call,
                                       double posterior,
                                       const ReferenceGenome& reference,
                                       const ReadMap& reads,
                                       const GenomicRegion& phase_region)
{
    using std::to_string;
    
    auto result = VcfRecord::Builder();
    
    auto phred_quality = Maths::probability_to_phred(posterior);
    
    const auto& region = mapped_region(reference_allele);
    
    result.set_chromosome(contig_name(region));
    result.set_position(region_begin(region));
    result.set_ref_allele(reference_allele.get_sequence());
    result.set_alt_alleles(extract_alt_allele_sequences(variants));
    result.set_quality(phred_quality);
    
    result.add_info("NS",  to_string(count_samples_with_coverage(reads, region)));
    result.add_info("DP",  to_string(sum_max_coverages(reads, region)));
    result.add_info("SB",  Octopus::to_string(strand_bias(reads, region), 2));
    result.add_info("BQ",  to_string(static_cast<unsigned>(rmq_base_quality(reads, region))));
    result.add_info("MQ",  to_string(static_cast<unsigned>(rmq_mapping_quality(reads, region))));
    result.add_info("MQ0", to_string(count_mapq_zero(reads, region)));
    
    result.set_format({"GT", "FT", "GP", "PS", "PQ", "DP", "BQ", "MQ"});
    
    auto vcf_genotype = to_vcf_genotype(genotype_call.genotype);
    
    for (const auto& sample_reads : reads) {
        const auto& sample = sample_reads.first;
        result.add_genotype(sample, vcf_genotype, VcfRecord::Builder::Phasing::Phased);
        
        result.add_genotype_field(sample, "FT", "."); // TODO
        result.add_genotype_field(sample, "GP", to_string(Maths::probability_to_phred(genotype_call.posterior)));
        result.add_genotype_field(sample, "PS", to_string(region_begin(phase_region) + 1));
        result.add_genotype_field(sample, "PQ", "60"); // TODO
        result.add_genotype_field(sample, "DP", to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result;
}

VcfRecord::Builder output_somatic_variant_call(const Allele& somatic_mutation,
                                               double posterior,
                                               const ReferenceGenome& reference,
                                               const ReadMap& reads)
{
    using std::to_string;
    
    auto result = VcfRecord::Builder();
    
    auto phred_quality = Maths::probability_to_phred(posterior);
    
    const auto& region = mapped_region(somatic_mutation);
    
    result.set_chromosome(contig_name(region));
    result.set_position(region_begin(region));
    result.set_ref_allele(reference.get_sequence(region));
    result.set_alt_allele(somatic_mutation.get_sequence());
    result.set_quality(phred_quality);
    
    result.add_info("SOMATIC", {});
    result.add_info("NS",  to_string(count_samples_with_coverage(reads, region)));
    result.add_info("DP",  to_string(sum_max_coverages(reads, region)));
    result.add_info("SB",  Octopus::to_string(strand_bias(reads, region), 2));
    result.add_info("BQ",  to_string(static_cast<unsigned>(rmq_base_quality(reads, region))));
    result.add_info("MQ",  to_string(static_cast<unsigned>(rmq_mapping_quality(reads, region))));
    result.add_info("MQ0", to_string(count_mapq_zero(reads, region)));
    
    result.set_format({"DP", "BQ", "MQ"});
    
    for (const auto& sample_reads : reads) {
        const auto& sample = sample_reads.first;
        result.add_genotype_field(sample, "DP", to_string(max_coverage(sample_reads.second, region)));
        result.add_genotype_field(sample, "BQ", to_string(static_cast<unsigned>(rmq_base_quality(sample_reads.second, region))));
        result.add_genotype_field(sample, "MQ", to_string(static_cast<unsigned>(rmq_mapping_quality(sample_reads.second, region))));
    }
    
    return result;
}

VcfRecord::Builder output_reference_call(RefCall call, ReferenceGenome& reference, const ReadMap& reads)
{
    using std::to_string;
    
    auto result = VcfRecord::Builder();
    
    auto phred_quality = Maths::probability_to_phred(call.posterior);
    
    const auto& region = mapped_region(call.reference_allele);
    
    result.set_chromosome(contig_name(region));
    result.set_position(region_begin(region));
    result.set_ref_allele(call.reference_allele.get_sequence().front());
    
    result.set_quality(phred_quality);
    
    result.set_filters({"REFCALL"});
    if (region_size(region) > 1) {
        result.add_info("END", to_string(region_end(region))); // - 1 as VCF uses closed intervals
    }
    
    result.add_info("NS",  to_string(count_samples_with_coverage(reads, region)));
    result.add_info("DP",  to_string(sum_max_coverages(reads, region)));
    result.add_info("SB",  Octopus::to_string(strand_bias(reads, region), 2));
    result.add_info("BQ",  to_string(static_cast<unsigned>(rmq_base_quality(reads, region))));
    result.add_info("MQ",  to_string(static_cast<unsigned>(rmq_mapping_quality(reads, region))));
    result.add_info("MQ0", to_string(count_mapq_zero(reads, region)));
    
    result.set_format({"GT", "GP", "DP", "BQ", "MQ"});
    
    for (const auto& sample_posteior : call.sample_posteriors) {
        const auto& sample = sample_posteior.first;
        result.add_homozygous_ref_genotype(sample, 2);
        result.add_genotype_field(sample, "GP", to_string(Maths::probability_to_phred(sample_posteior.second)));
        result.add_genotype_field(sample, "DP", to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result;
}
} // namespace
    
std::vector<VcfRecord::Builder>
CancerVariantCaller::call_variants(const std::vector<Variant>& candidates,
                                   const std::vector<Allele>& callable_alleles,
                                   CallerLatents* latents,
                                   const Phaser::PhaseSet& phase_set,
                                   const ReadMap& reads) const
{
    const auto dlatents = dynamic_cast<Latents*>(latents);
    
    const auto model_posteriors = calculate_model_posteriors(*dlatents);
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "Germline model posterior: " << model_posteriors.germline;
        stream(log) << "CNV model posterior:      " << model_posteriors.cnv;
        stream(log) << "Somatic model posterior:  " << model_posteriors.somatic;
    }
    
    if (model_posteriors.somatic > model_posteriors.germline) {
        if (model_posteriors.somatic > model_posteriors.cnv) {
            return call_somatic_variants(candidates, callable_alleles,
                                         dlatents->somatic_model_inferences_.posteriors,
                                         phase_set, reads);
        } else {
            return call_cnv_variants(candidates, callable_alleles,
                                     dlatents->cnv_model_inferences_.posteriors,
                                     phase_set, reads);
        }
    } else if (model_posteriors.cnv > model_posteriors.germline) {
        return call_cnv_variants(candidates, callable_alleles,
                                 dlatents->cnv_model_inferences_.posteriors,
                                 phase_set, reads);
    } else {
        return call_germline_variants(candidates, callable_alleles,
                                      dlatents->germline_model_inferences_.posteriors,
                                      phase_set, reads);
    }
}

CancerVariantCaller::ModelPosteriors::ModelPosteriors(double germline, double cnv, double somatic)
: germline {germline}, cnv {cnv}, somatic {somatic} {}

CancerVariantCaller::ModelPosteriors
CancerVariantCaller::calculate_model_posteriors(const Latents& inferences) const
{
    const auto& germline_inferences = inferences.germline_model_inferences_;
    const auto& cnv_inferences = inferences.cnv_model_inferences_;
    const auto& somatic_inferences = inferences.somatic_model_inferences_;
    
    double germline_model_prior {0.89999};
    double cnv_model_prior {0.1};
    double somatic_model_prior {0.01 * 0.001};
    
    auto germline_model_jlp = std::log(germline_model_prior) + germline_inferences.log_evidence;
    auto cnv_model_jlp = std::log(cnv_model_prior) + cnv_inferences.approx_log_evidence;
    auto somatic_model_jlp = std::log(somatic_model_prior) + somatic_inferences.approx_log_evidence;
    
    auto norm = Maths::log_sum_exp(germline_model_jlp, cnv_model_jlp, somatic_model_jlp);
    
    auto germline_model_posterior = std::exp(germline_model_jlp - norm);
    auto cnv_model_posterior = std::exp(cnv_model_jlp - norm);
    auto somatic_model_posterior = std::exp(somatic_model_jlp - norm);
    
    return ModelPosteriors {germline_model_posterior, cnv_model_posterior, somatic_model_posterior};
}

std::vector<VcfRecord::Builder>
CancerVariantCaller::call_germline_variants(const std::vector<Variant>& candidates,
                                            const std::vector<Allele>& callable_alleles,
                                            const GenotypeModel::Individual::Latents& posteriors,
                                            const Phaser::PhaseSet& phase_set,
                                            const ReadMap& reads) const
{
    std::vector<VcfRecord::Builder> result {};
    
    if (parameters_.call_somatics_only) {
        return result;
    }
    
    return result;
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

std::vector<VcfRecord::Builder>
CancerVariantCaller::call_cnv_variants(const std::vector<Variant>& candidates,
                                       const std::vector<Allele>& callable_alleles,
                                       const GenotypeModel::CNV::Latents& posteriors,
                                       const Phaser::PhaseSet& phase_set,
                                       const ReadMap& reads) const
{
    const auto credible_regions = compute_marginal_credible_intervals(posteriors.alphas, 0.99);
    
    for (const auto& p : credible_regions) {
        std::cout << p.first << std::endl;
        for (const auto& r : p.second) {
            std::cout << "\t" << r.first << " " << r.second << '\n';
        }
    }
    
    std::vector<VcfRecord::Builder> result {};
    
    return result;
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

std::vector<VcfRecord::Builder>
CancerVariantCaller::call_somatic_variants(const std::vector<Variant>& candidates,
                                           const std::vector<Allele>& callable_alleles,
                                           const GenotypeModel::Somatic::Latents& posteriors,
                                           const Phaser::PhaseSet& phase_set,
                                           const ReadMap& reads) const
{
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        debug::print_map_genotype(stream(log), posteriors);
    }
    
    const auto credible_regions = compute_marginal_credible_intervals(posteriors.alphas, 0.99);
    
    for (const auto& p : credible_regions) {
        std::cout << p.first << std::endl;
        for (const auto& r : p.second) {
            std::cout << "\t" << r.first << " " << r.second << '\n';
        }
    }
    
    std::vector<VcfRecord::Builder> result {};
    
    return result;
}
} // namespace Octopus
