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
#include "coalescent_model.hpp"
#include "somatic_model.hpp"
#include "individual_genotype_model.hpp"
#include "cnv_genotype_model.hpp"

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
call_somatics_only {call_somatics_only}
{}

namespace
{
    using GM = GenotypeModel::Cancer;
}

CancerVariantCaller::CancerVariantCaller(const ReferenceGenome& reference,
                                         ReadPipe& read_pipe,
                                         CandidateVariantGenerator&& candidate_generator,
                                         VariantCaller::CallerParameters general_parameters,
                                         CallerParameters specific_parameters)
:
VariantCaller {reference, read_pipe, std::move(candidate_generator), std::move(general_parameters)},
parameters_ {std::move(specific_parameters)}
{}

CancerVariantCaller::Latents::Latents(ModelLatents&& model_latents)
:
model_latents_ {std::move(model_latents)}
{}

std::shared_ptr<CancerVariantCaller::Latents::HaplotypeProbabilityMap>
CancerVariantCaller::Latents::get_haplotype_posteriors() const
{
    return nullptr;
}

std::shared_ptr<CancerVariantCaller::Latents::GenotypeProbabilityMap>
CancerVariantCaller::Latents::get_genotype_posteriors() const
{
    return nullptr;
}
    
// private methods

namespace
{

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
    
    //result.set_filters({"PASS"}); // TODO
    
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

std::unique_ptr<CancerVariantCaller::CallerLatents>
CancerVariantCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                                   const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    Haplotype reference_haplotype {mapped_region(haplotypes.front()), reference_};
    
    CoalescentModel germline_prior_model {reference_haplotype};
    
    GenotypeModel::Individual germline_model {parameters_.ploidy, germline_prior_model};
    
    SomaticModel cancer_prior_model {germline_prior_model, parameters_.somatic_mutation_rate};
    
    GenotypeModel::CNV::Priors::GenotypeMixturesDirichletAlphaMap cnv_alphas {};
    cnv_alphas.reserve(samples_.size());
    for (const auto& sample : samples_) {
        if (parameters_.normal_sample && sample == *parameters_.normal_sample) {
            GenotypeModel::CNV::Priors::GenotypeMixturesDirichletAlphas sample_alphas {10.0, 10.0};
            cnv_alphas.emplace(sample, std::move(sample_alphas));
        } else {
            GenotypeModel::CNV::Priors::GenotypeMixturesDirichletAlphas sample_alphas {0.5, 0.5};
            cnv_alphas.emplace(sample, std::move(sample_alphas));
        }
    }
    GenotypeModel::CNV::Priors cnv_priors {germline_prior_model, std::move(cnv_alphas)};
    GenotypeModel::CNV cnv_model {samples_, parameters_.ploidy, std::move(cnv_priors)};
    
    GM::Priors::GenotypeMixturesDirichletAlphaMap alphas {};
    alphas.reserve(samples_.size());
    
    for (const auto& sample : samples_) {
        if (parameters_.normal_sample && sample == *parameters_.normal_sample) {
            GM::Priors::GenotypeMixturesDirichletAlphas sample_alphas {5.0, 5.0, 0.01};
            alphas.emplace(sample, std::move(sample_alphas));
        } else {
            GM::Priors::GenotypeMixturesDirichletAlphas sample_alphas {1.0, 1.0, 0.75};
            alphas.emplace(sample, std::move(sample_alphas));
        }
    }
    
    GM::Priors priors {cancer_prior_model, std::move(alphas)};
    
    GM::Cancer somatic_model {samples_, parameters_.ploidy, std::move(priors)};
    
    std::vector<CancerGenotype<Haplotype>> cancer_genotypes;
    std::vector<Genotype<Haplotype>> germline_genotypes;
    
    std::tie(cancer_genotypes, germline_genotypes) = generate_all_cancer_genotypes(haplotypes, parameters_.ploidy);
    
    std::cout << "there are " << haplotypes.size() << " haplotypes" << std::endl;
    std::cout << "there are " << germline_genotypes.size() << " germline genotypes" << std::endl;
    std::cout << "there are " << cancer_genotypes.size() << " cancer genotypes" << std::endl;
    
    auto merged_likelihoods = merge_samples(samples_, "union", haplotypes, haplotype_likelihoods);
    
    auto germline_inferences = germline_model.infer_latents("union", germline_genotypes,
                                                            merged_likelihoods);
    
    
    auto cnv_inferences = cnv_model.infer_latents(germline_genotypes, haplotype_likelihoods);
    auto somatic_inferences = somatic_model.infer_latents(cancer_genotypes, haplotype_likelihoods);
    
    std::cout << "Log evidence for germline model " << germline_inferences.log_evidence << std::endl;
    std::cout << "Log evidence for CNV model " << cnv_inferences.approx_log_evidence << std::endl;
    std::cout << "Log evidence for cancer model " << somatic_inferences.approx_log_evidence << std::endl;
    
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
    
    std::cout << "germline model posterior = " << germline_model_posterior << std::endl;
    std::cout << "CNV model posterior = " << cnv_model_posterior << std::endl;
    std::cout << "somatic model posterior = " << somatic_model_posterior << std::endl;
    
    const double cutoff {0.001};
    
    std::cout << "CNV genotypes with posterior > " << cutoff << ":" << std::endl;
    for (const auto& p : cnv_inferences.posteriors.genotype_probabilities) {
        if (p.second > cutoff) {
            ::debug::print_variant_alleles(p.first);
            std::cout << ' ' << p.second << std::endl;
        }
    }
    
    for (const auto& sample : samples_) {
        auto a = cnv_inferences.posteriors.alphas.at(sample);
        std::cout << sample << " " << a[0] << " " << a[1] << std::endl;
        auto a0 = std::accumulate(std::cbegin(a), std::cend(a), 0.0);
        const double m {0.99};
        auto p0 = Maths::beta_hdi(a[0], a0 - a[0], m);
        std::cout << "a[0] " <<  100 * m << "% credible region = (" << p0.first << ", " << p0.second << ")" << std::endl;
        auto p1 = Maths::beta_hdi(a[1], a0 - a[1], m);
        std::cout << "a[1] " <<  100 * m << "% credible region = (" << p1.first << ", " << p1.second << ")" << std::endl;
    }
    
    auto it = std::max_element(std::cbegin(somatic_inferences.posteriors.genotype_probabilities),
                               std::cend(somatic_inferences.posteriors.genotype_probabilities),
                               [] (const auto& lhs, const auto& rhs) {
                                   return lhs.second < rhs.second;
                               });
    
    std::cout << "MAP cancer genotype is ";
    debug::print_variant_alleles(it->first);
    std::cout << ' ' << it->second << std::endl;
    
    if (haplotypes.size() <= parameters_.ploidy && parameters_.normal_sample) {
        const auto& map_cancer_genotype = it->first;
        
        GenotypeModel::Individual individual_model {parameters_.ploidy, germline_prior_model};
        
        auto map_germline_inferences = individual_model.infer_latents(*parameters_.normal_sample,
                                                                      germline_genotypes,
                                                                      haplotype_likelihoods);
        
        auto it2 = std::find(std::cbegin(germline_genotypes), std::cend(germline_genotypes),
                             map_cancer_genotype.get_germline_genotype());
        
        auto p = map_germline_inferences.posteriors.genotype_probabilities[std::distance(std::cbegin(germline_genotypes), it2)];
        
        std::cout << "Posterior of map cancer genotype germline for individual model = "
        << p << std::endl;
        
        std::cout << "Germline model log evidence for normal sample = " << map_germline_inferences.log_evidence << std::endl;
    }
    
    std::cout << "Cancer genotypes with posterior > " << cutoff << ":" << std::endl;
    for (const auto& p : somatic_inferences.posteriors.genotype_probabilities) {
        if (p.second > cutoff) {
            debug::print_variant_alleles(p.first);
            std::cout << ' ' << p.second << std::endl;
        }
    }
    
    for (const auto& sample : samples_) {
        std::cout << sample << "";
        if (sample == parameters_.normal_sample) {
            std::cout << " (normal)";
        }
        std::cout << std::endl;
        auto posterior_alphas = somatic_inferences.posteriors.alphas.at(sample);
        std::cout << std::setprecision(3);
        std::cout << "\talphas: " << posterior_alphas[0] << " " << posterior_alphas[1] << " " << posterior_alphas[2] << std::endl;
        auto a0 = std::accumulate(std::cbegin(posterior_alphas), std::cend(posterior_alphas), 0.0);
        const double m {0.99};
        auto p = Maths::beta_hdi(posterior_alphas[2], a0 - posterior_alphas[2], m);
        std::cout << "\tsomatic probability " << 100 * m << "% credible region = ("
        << p.first << ", " << p.second << ")" << std::endl;
    }
    
    throw std::runtime_error {"whoops"};
    
    return nullptr;
    //return std::make_unique<Latents>(std::move(inferred_latents));
}

std::vector<VcfRecord::Builder>
CancerVariantCaller::call_variants(const std::vector<Variant>& candidates,
                                   const std::vector<Allele>& callable_alleles,
                                   CallerLatents* latents,
                                   const Phaser::PhaseSet& phase_set,
                                   const ReadMap& reads) const
{
    const auto dlatents = dynamic_cast<Latents*>(latents);
    
    std::vector<VcfRecord::Builder> result {};
    
    return result;
}
} // namespace Octopus
