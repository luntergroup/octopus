//
//  cancer_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 16/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_caller.hpp"

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <iterator>
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
#include "haplotype_prior_model.hpp"
#include "probability_matrix.hpp"

namespace Octopus
{
// public methods

CancerVariantCaller::CancerVariantCaller(const ReferenceGenome& reference,
                                         ReadPipe& read_pipe,
                                         CandidateVariantGenerator&& candidate_generator,
                                         RefCallType refcall_type, double min_variant_posterior,
                                         double min_somatic_posterior, double min_refcall_posterior,
                                         const SampleIdType& normal_sample, bool call_somatics_only)
:
VariantCaller {reference, read_pipe, std::move(candidate_generator), refcall_type},
genotype_model_ {normal_sample},
normal_sample_ {normal_sample},
min_variant_posterior_ {min_variant_posterior},
min_somatic_mutation_posterior_ {min_somatic_posterior},
min_refcall_posterior_ {min_refcall_posterior},
call_somatics_only_ {call_somatics_only}
{}

namespace
{
    using GM = GenotypeModel::Cancer;
}

// private methods

namespace
{
using GermlineGenotypeMarginalMap = std::unordered_map<Genotype<Haplotype>, double>;
using CancerHaplotypeMarginalMap  = std::unordered_map<Haplotype, double>;

using MarginalAllelePosteriorMap = std::unordered_map<Allele, double>;
using AllelePosteriorMap         = std::unordered_map<SampleIdType, MarginalAllelePosteriorMap>;

struct VariantCall
{
    VariantCall() = default;
    template <typename T>
    VariantCall(T&& variants, double posterior) : variants {std::forward<T>(variants)}, posterior {posterior} {}
    
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

std::ostream& operator<<(std::ostream& os, const std::array<double, 3>& arr)
{
    os << arr[0] << " " << arr[1] << " " << arr[2];
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::unordered_map<SampleIdType, std::array<double, 3>>& m)
{
    for (const auto& p : m) os << p.first << ": " << p.second << "\n";
    return os;
}

std::vector<GenomicRegion> get_segment_regions(const std::vector<std::vector<Variant>>& segments)
{
    std::vector<GenomicRegion> result {};
    result.reserve(segments.size());
    for (const auto& segment : segments) result.push_back(segment.front().get_region());
    return result;
}

auto find_map_genotype(const GM::GenotypeProbabilityMap& genotype_posteriors)
{
    return *std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                             [] (const auto& lhs, const auto& rhs) {
                                 return lhs.second < rhs.second;
                             });
}

template <typename Map>
auto find_map_genotype(const Map& map)
{
    return *std::max_element(std::cbegin(map), std::cend(map),
                             [] (const auto& lhs, const auto& rhs) {
                                 return lhs.second < rhs.second;
                             });
}

GermlineGenotypeMarginalMap
marginalise_germline_genotypes(const GM::GenotypeProbabilityMap& genotype_posteriors,
                               unsigned num_haplotypes)
{
    GermlineGenotypeMarginalMap result {};
    result.reserve(num_genotypes(num_haplotypes, 2));
    
    for (const auto& genotype_posterior : genotype_posteriors) {
        result[genotype_posterior.first.get_germline_genotype()] += genotype_posterior.second;
    }
    
    return result;
}

CancerHaplotypeMarginalMap
marginalise_cancer_haplotypes(const GM::GenotypeProbabilityMap& genotype_posteriors,
                              unsigned num_haplotypes)
{
    CancerHaplotypeMarginalMap result {};
    result.reserve(num_haplotypes);
    
    for (const auto& genotype_posterior : genotype_posteriors) {
        result[genotype_posterior.first.get_cancer_element()] += genotype_posterior.second;
    }
    
    return result;
}

static double marginalise(const Allele& allele, const GermlineGenotypeMarginalMap& germline_genotype_posteriors)
{
    double result {0.0};
    
    for (const auto& genotype_posterior : germline_genotype_posteriors) {
        if (contains(genotype_posterior.first, allele)) result += genotype_posterior.second;
    }
    
    return result;
}

MarginalAllelePosteriorMap
compute_germline_allele_posteriors(const GermlineGenotypeMarginalMap& germline_genotype_posteriors,
                                   const std::vector<Allele>& alleles)
{
    MarginalAllelePosteriorMap result {};
    result.reserve(alleles.size());
    
    for (const auto& allele : alleles) {
        result.emplace(allele, marginalise(allele, germline_genotype_posteriors));
    }
    
    return result;
}

static double marginalise(const Allele& allele, const CancerHaplotypeMarginalMap& cancer_haplotype_posteriors)
{
    double result {0.0};
    
    for (const auto& haplotype_posterior : cancer_haplotype_posteriors) {
        if (contains(haplotype_posterior.first, allele)) result += haplotype_posterior.second;
    }
    
    return result;
}

MarginalAllelePosteriorMap
compute_cancer_allele_posteriors(const CancerHaplotypeMarginalMap& cancer_haplotype_posteriors,
                                 const std::vector<Allele>& alleles)
{
    MarginalAllelePosteriorMap result {};
    result.reserve(alleles.size());
    
    for (const auto& allele : alleles) {
        result.emplace(allele, marginalise(allele, cancer_haplotype_posteriors));
    }
    
    return result;
}

AllelePosteriorMap
compute_germline_allele_posteriors(const GermlineGenotypeMarginalMap& germline_genotype_posteriors,
                                   const GM::GenotypeMixtures& genotype_mixtures,
                                   const std::vector<Allele>& alleles)
{
    AllelePosteriorMap result {};
    result.reserve(alleles.size());
    
    for (const auto& sample_mixtures : genotype_mixtures) {
        MarginalAllelePosteriorMap marginals {};
        marginals.reserve(alleles.size());
        
        const double germline_fraction = 1.0 - genotype_mixtures.at(sample_mixtures.first).back();
        
        for (const auto& allele : alleles) {
            marginals.emplace(allele, germline_fraction * marginalise(allele, germline_genotype_posteriors));
        }
        
        result.emplace(sample_mixtures.first, std::move(marginals));
    }
    
    return result;
}

AllelePosteriorMap
compute_cancer_allele_posteriors(const CancerHaplotypeMarginalMap& cancer_haplotype_posteriors,
                                 const GM::GenotypeMixtures& genotype_mixtures,
                                 const std::vector<Allele>& alleles)
{
    AllelePosteriorMap result {};
    result.reserve(genotype_mixtures.size());
    
    for (const auto& sample_mixtures : genotype_mixtures) {
        MarginalAllelePosteriorMap marginals {};
        marginals.reserve(alleles.size());
        
        const double cancer_fraction    = genotype_mixtures.at(sample_mixtures.first).back();
        const double cancer_probability = (2 * cancer_fraction) / ((1.0 - cancer_fraction) + 2 * cancer_fraction);
        
        for (const auto& allele : alleles) {
            marginals.emplace(allele, cancer_probability * marginalise(allele, cancer_haplotype_posteriors));
        }
        
        result.emplace(sample_mixtures.first, std::move(marginals));
    }
    
    return result;
}

GermlineGenotypeCalls
call_germline_genotypes(const Genotype<Haplotype>& map_germline_genotype,
                        const std::vector<GenomicRegion>& segments)
{
    GermlineGenotypeCalls result {};
    result.reserve(segments.size());
    
    for (const auto& region : segments) {
        result.emplace_back(splice<Allele>(map_germline_genotype, region), 1.0);
    }
    
    return result;
}

double max_posterior(const std::vector<Variant>& variants, const MarginalAllelePosteriorMap& allele_posteriors)
{
    double result {};
    
    for (const auto& variant : variants) {
        auto curr = allele_posteriors.at(variant.get_alt_allele());
        if (curr > result) result = curr;
    }
    
    return result;
}

std::vector<Variant> call_segment_variants(const std::vector<Variant>& variants,
                                           const MarginalAllelePosteriorMap& allele_posteriors,
                                           const double min_posterior)
{
    std::vector<Variant> result {};
    result.reserve(variants.size());
    
    std::copy_if(std::cbegin(variants), std::cend(variants), std::back_inserter(result),
                 [&allele_posteriors, min_posterior] (const auto& variant) {
                     return allele_posteriors.at(variant.get_alt_allele()) >= min_posterior;
                 });
    
    result.shrink_to_fit();
    
    return result;
}

VariantCalls call_germline_variants(const std::vector<std::vector<Variant>>& segments,
                                    const MarginalAllelePosteriorMap& normal_sample_germline_allele_posteriors,
                                    const double min_posterior)
{
    VariantCalls result {};
    
    for (const auto& segment : segments) {
        auto calls = call_segment_variants(segment, normal_sample_germline_allele_posteriors, min_posterior);
        if (!calls.empty()) {
            result.emplace_back(std::move(calls), max_posterior(calls, normal_sample_germline_allele_posteriors));
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

double max_posterior(const Allele& allele, const AllelePosteriorMap& allele_posteriors)
{
    double result {0.0};
    
    for (const auto& sample_allele_posteriors : allele_posteriors) {
        auto p = sample_allele_posteriors.second.at(allele);
        if (p > result) result = p;
    }
    
    return result;
}

std::unordered_set<Allele> get_called_germline_alleles(const GermlineGenotypeCalls& germline_genotype_calls)
{
    std::unordered_set<Allele> result {};
    
    if (germline_genotype_calls.empty()) return result;
    
    result.reserve(germline_genotype_calls.size() * germline_genotype_calls.front().genotype.ploidy());
    
    for (const auto& call : germline_genotype_calls) {
        auto alleles = call.genotype.copy_unique();
        result.insert(std::make_move_iterator(std::begin(alleles)), std::make_move_iterator(std::end(alleles)));
    }
    
    return result;
}

bool is_called_allele(const Allele& allele, const std::unordered_set<Allele>& called_alleles)
{
    return called_alleles.count(allele) == 1;
}

// TODO: this doesn't quite work
SomaticCalls call_somatic_mutations(const std::vector<Allele>& alleles,
                                    const AllelePosteriorMap& cancer_allele_posteriors,
                                    const GermlineGenotypeCalls& germline_genotype_calls,
                                    const double min_posterior)
{
    SomaticCalls result {};
    
    auto called_germline_alleles = get_called_germline_alleles(germline_genotype_calls);
    
    result.reserve(alleles.size() - called_germline_alleles.size());
    
    for (const auto& allele : alleles) {
        if (!is_called_allele(allele, called_germline_alleles)) {
            auto max_sample_allele_posterior = max_posterior(allele, cancer_allele_posteriors);
            
            if (max_sample_allele_posterior >= min_posterior) {
                result.emplace_back(allele, max_sample_allele_posterior);
            }
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

void parsimonise_germline_variant_calls(VariantCalls& germline_variant_calls,
                                        const ReferenceGenome& reference)
{
    for (auto& call : germline_variant_calls) {
        call.variants = parsimonise_together(call.variants, reference);
    }
}

void parsimonise_somatic_mutation_calls(SomaticCalls& somatic_mutation_calls,
                                        const ReferenceGenome& reference)
{
    for (auto& call : somatic_mutation_calls) {
        if (is_indel(call.allele) && !is_reference(call.allele, reference)) {
            Variant v {make_reference_allele(call.allele.get_region(), reference), call.allele};
            call.allele = make_parsimonious(v, reference).get_alt_allele();
        }
    }
}

std::vector<GenomicRegion> get_called_regions(const VariantCalls& germline_variant_calls,
                                              const SomaticCalls& somatic_mutation_calls)
{
    std::vector<GenomicRegion> result {};
    result.reserve(germline_variant_calls.size() + somatic_mutation_calls.size());
    
    for (const auto& call : germline_variant_calls) {
        result.emplace_back(encompassing_region(call.variants));
    }
    
    for (const auto& call : somatic_mutation_calls) {
        result.emplace_back(call.allele.get_region());
    }
    
    std::sort(std::begin(result), std::end(result));
    
    result.shrink_to_fit();
    
    return result;
}

auto get_regions(const VariantCalls& variant_calls)
{
    std::vector<GenomicRegion> result {};
    result.reserve(variant_calls.size());
    
    for (const auto& segment_calls : variant_calls) {
        result.emplace_back(encompassing_region(segment_calls.variants));
    }
    
    return result;
}

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

//std::vector<VcfRecord>
//CancerVariantCaller::call_variants(const GenomicRegion& region,
//                                   const std::vector<Variant>& candidates,
//                                   const ReadMap& reads) const
//{
//    std::vector<VcfRecord> result {};
//    
//    if (is_empty_region(region) || (candidates.empty() && refcall_type_ == RefCallType::None)) {
//        return result;
//    }
//    
//    HaplotypePhaser phaser {region.get_contig_name(), reference_, candidates, reads, 64, 2};
//    
//    while (!phaser.done()) {
//        auto haplotypes = phaser.get_haplotypes();
//        
//        make_unique(haplotypes, haplotype_prior_model_);
//        
//        phaser.set_haplotypes(haplotypes);
//        
//        std::cout << "there are " << haplotypes.size() << " unique haplotypes" << std::endl;
//        
//        auto haplotype_region = mapped_region(haplotypes.front());
//        
//        std::cout << "haplotype region is " << haplotype_region << std::endl;
//        
//        auto haplotype_region_reads = copy_overlapped(reads, haplotype_region);
//        
//        std::cout << "there are " << count_reads(haplotype_region_reads) << " reads in haplotype region" << std::endl;
//        
//        auto latents = genotype_model_.infer_latents(haplotypes, haplotype_region_reads, reference_);
//        
//        if (latents.genotype_posteriors.empty()) return result;
//        
//        //remove_low_posteriors(latents.genotype_posteriors, 1.0e-06);
//        
//        auto num_haplotypes = static_cast<unsigned>(haplotypes.size());
//        
//        auto map_cancer_genotype = find_map_genotype(latents.genotype_posteriors);
//        
//        std::cout << "map cancer genotype: " << std::endl;
//        print_variant_alleles(map_cancer_genotype.first);
//        std::cout << " " << map_cancer_genotype.second << std::endl;
//        
//        std::cout << "posterior genotype mixtures: " << std::endl;
//        for (const auto& sw : latents.genotype_mixtures) {
//            std::cout << sw.first << ": " << sw.second << std::endl;
//        }
//        
//        auto germline_genotype_posteriors = marginalise_germline_genotypes(latents.genotype_posteriors, num_haplotypes);
//        
//        auto phased_gps = phaser.phase(haplotypes, {{"germline", germline_genotype_posteriors}});
//        
//        auto map_germline_genotype = find_map_genotype(germline_genotype_posteriors);
//        
//        std::cout << "germline genotype posteriors: " << std::endl;
//        for (const auto& gp : germline_genotype_posteriors) {
//            print_variant_alleles(gp.first);
//            std::cout << " " << gp.second << std::endl;
//        }
//        
//        if (map_germline_genotype.second < 0.9) {
//            std::cout << "too much uncertainty in germline genotype to call variants" << std::endl;
//            continue;
//        }
//        
//        auto cancer_haplotype_posteriors  = marginalise_cancer_haplotypes(latents.genotype_posteriors, num_haplotypes);
//        
////            std::cout << "cancer haplotype posteriors: " << std::endl;
////            for (const auto& hp : cancer_haplotype_posteriors) {
////                print_variant_alleles(hp.first);
////                std::cout << " " << hp.second << std::endl;
////            }
//        
//        auto alleles = generate_callable_alleles(region, candidates, refcall_type_, reference_);
//        
//        auto germline_allele_posteriors = compute_germline_allele_posteriors(germline_genotype_posteriors,
//                                                                             latents.genotype_mixtures, alleles);
//        
//        for (const auto& sa : germline_allele_posteriors) {
//            std::cout << "germline allele posteriors for sample: " << sa.first << std::endl;
//            for (const auto& ap : sa.second) {
//                std::cout << ap.first << " " << ap.second << std::endl;
//            }
//        }
//        
//        auto cancer_allele_posteriors = compute_cancer_allele_posteriors(cancer_haplotype_posteriors,
//                                                                         latents.genotype_mixtures, alleles);
//        
//        for (const auto& sa : cancer_allele_posteriors) {
//            std::cout << "cancer allele posteriors for sample: " << sa.first << std::endl;
//            for (const auto& ap : sa.second) {
//                std::cout << ap.first << " " << ap.second << std::endl;
//            }
//        }
//        
//        auto segments = segment_overlapped(candidates);
//        
//        auto segment_regions = get_segment_regions(segments);
//        
//        auto germline_variant_calls = call_germline_variants(segments, germline_allele_posteriors.at(normal_sample_),
//                                                             min_variant_posterior_);
//        
//        auto germline_genotype_calls = call_germline_genotypes(map_germline_genotype.first, get_regions(germline_variant_calls));
//        
//        std::cout << "called germline variants" << std::endl;
//        for (const auto& call : germline_variant_calls) {
//            for (const auto& variant : call.variants) {
//                std::cout << variant << " " << call.posterior << std::endl;
//            }
//        }
//      
//        std::cout << "germline genotype (allele) calls" << std::endl;
//        for (const auto& call : germline_genotype_calls) {
//            std::cout << call.genotype << std::endl;
//        }
//        
//        // need germline genotype calls to catch reversion to reference somatics
//        auto somatic_mutation_calls = call_somatic_mutations(alleles, cancer_allele_posteriors,
//                                                             germline_genotype_calls,
//                                                             min_somatic_mutation_posterior_);
//        
//        std::cout << "called somatic mutations" << std::endl;
//        for (const auto& call : somatic_mutation_calls) {
//            std::cout << call.allele << " " << call.posterior << std::endl;
//        }
//        
//        parsimonise_germline_variant_calls(germline_variant_calls, reference_);
//        parsimonise_somatic_mutation_calls(somatic_mutation_calls, reference_);
//        
//        auto called_regions = get_called_regions(germline_variant_calls, somatic_mutation_calls);
//        
//        auto candidate_refcall_alleles = generate_candidate_reference_alleles(alleles, called_regions,
//                                                                              candidates, refcall_type_);
//        
////            std::cout << "candidate refcall alleles:" << std::endl;
////            for (const auto& allele : candidate_refcall_alleles) {
////                std::cout << allele << std::endl;
////            }
//        
//        const auto phase_region = (called_regions.empty()) ? head_region(region) : encompassing_region(called_regions.front(), called_regions.back());
//        
//        result.reserve(germline_variant_calls.size() + somatic_mutation_calls.size());
//        
//        if (call_somatics_only_) {
//            germline_variant_calls.clear();
//            germline_genotype_calls.clear();
//        } else {
//            germline_genotype_calls = call_germline_genotypes(map_germline_genotype.first,
//                                                              get_regions(germline_variant_calls));
//        }
//        
//        merge_transform(germline_variant_calls, germline_genotype_calls,
//                        somatic_mutation_calls, std::back_inserter(result),
//                        [this, &reads, &phase_region] (const auto& variant_call, const auto& genotype_call) {
//                            return output_germline_variant_call(variant_call.variants.front().get_ref_allele(),
//                                                                variant_call.variants, genotype_call,
//                                                                variant_call.posterior, reference_, reads, phase_region);
//                        },
//                        [this, &reads] (const auto& somatic_call) {
//                            return output_somatic_variant_call(somatic_call.allele, somatic_call.posterior, reference_, reads);
//                        },
//                        [] (const auto& lhs, const auto& rhs) {
//                            return is_before(lhs.allele, rhs.variants.front());
//                        });
//    }
//    
//    return result;
//}

std::unique_ptr<CancerVariantCaller::CallerLatents>
CancerVariantCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                                   const HaplotypePriorMap& haplotype_priors,
                                   HaplotypeLikelihoodCache& haplotype_likelihoods,
                                   const ReadMap& reads) const
{
    auto model_latents = genotype_model_.infer_latents(haplotypes, haplotype_priors,
                                                       haplotype_likelihoods, reads);
    
    auto num_haplotypes = static_cast<unsigned>(haplotypes.size());
    
    auto germline_genotype_posteriors = marginalise_germline_genotypes(model_latents.genotype_posteriors,
                                                                       num_haplotypes);
    
    auto cancer_haplotype_posteriors  = marginalise_cancer_haplotypes(model_latents.genotype_posteriors,
                                                                      num_haplotypes);
    
    return std::make_unique<Latents>(germline_genotype_posteriors, cancer_haplotype_posteriors,
                                     model_latents.genotype_mixtures);
}

std::vector<VcfRecord::Builder>
CancerVariantCaller::call_variants(const std::vector<Variant>& candidates,
                                   const std::vector<Allele>& callable_alleles,
                                   CallerLatents* latents,
                                   const HaplotypePhaser::PhaseSet& phase_set,
                                   const ReadMap& reads) const
{
    const auto dlatents = dynamic_cast<Latents*>(latents);
    
    auto germline_allele_posteriors = compute_germline_allele_posteriors(dlatents->germline_genotype_posteriors,
                                                                         dlatents->genotype_mixtures,
                                                                         callable_alleles);
    
//    for (const auto& sa : germline_allele_posteriors) {
//        std::cout << "germline allele posteriors for sample: " << sa.first << std::endl;
//        for (const auto& ap : sa.second) {
//            std::cout << ap.first << " " << ap.second << std::endl;
//        }
//    }
    
    auto cancer_allele_posteriors = compute_cancer_allele_posteriors(dlatents->cancer_haplotype_posteriors,
                                                                     dlatents->genotype_mixtures,
                                                                     callable_alleles);
    
//    for (const auto& sa : cancer_allele_posteriors) {
//        std::cout << "cancer allele posteriors for sample: " << sa.first << std::endl;
//        for (const auto& ap : sa.second) {
//            std::cout << ap.first << " " << ap.second << std::endl;
//        }
//    }
    
    auto segments = segment_overlapped(candidates);
    
    auto segment_regions = get_segment_regions(segments);
    
    auto germline_variant_calls = call_germline_variants(segments, germline_allele_posteriors.at(normal_sample_),
                                                         min_variant_posterior_);
    
    auto map_germline_genotype = find_map_genotype(dlatents->germline_genotype_posteriors);
    
    auto germline_genotype_calls = call_germline_genotypes(map_germline_genotype.first, get_regions(germline_variant_calls));
    
//    std::cout << "called germline variants" << std::endl;
//    for (const auto& call : germline_variant_calls) {
//        for (const auto& variant : call.variants) {
//            std::cout << variant << " " << call.posterior << std::endl;
//        }
//    }
//    
//    std::cout << "germline genotype (allele) calls" << std::endl;
//    for (const auto& call : germline_genotype_calls) {
//        std::cout << call.genotype << std::endl;
//    }
    
    // need germline genotype calls to catch reversion to reference somatics
    auto somatic_mutation_calls = call_somatic_mutations(callable_alleles, cancer_allele_posteriors,
                                                         germline_genotype_calls,
                                                         min_somatic_mutation_posterior_);
    
//    std::cout << "called somatic mutations" << std::endl;
//    for (const auto& call : somatic_mutation_calls) {
//        std::cout << call.allele << " " << call.posterior << std::endl;
//    }
    
    parsimonise_germline_variant_calls(germline_variant_calls, reference_);
    parsimonise_somatic_mutation_calls(somatic_mutation_calls, reference_);
    
    auto called_regions = get_called_regions(germline_variant_calls, somatic_mutation_calls);
    
    auto candidate_refcall_alleles = generate_candidate_reference_alleles(callable_alleles, called_regions,
                                                                          candidates, refcall_type_);
    
    //            std::cout << "candidate refcall alleles:" << std::endl;
    //            for (const auto& allele : candidate_refcall_alleles) {
    //                std::cout << allele << std::endl;
    //            }
    
    std::vector<VcfRecord::Builder> result {};
    
    result.reserve(germline_variant_calls.size() + somatic_mutation_calls.size());
    
    if (call_somatics_only_) {
        germline_variant_calls.clear();
        germline_genotype_calls.clear();
    } else {
        germline_genotype_calls = call_germline_genotypes(map_germline_genotype.first,
                                                          get_regions(germline_variant_calls));
    }
    
    auto phase_region = phase_set.region;
    
    merge_transform(germline_variant_calls, germline_genotype_calls,
                    somatic_mutation_calls, std::back_inserter(result),
                    [this, &reads, &phase_region] (const auto& variant_call, const auto& genotype_call) {
                        return output_germline_variant_call(variant_call.variants.front().get_ref_allele(),
                                                            variant_call.variants, genotype_call,
                                                            variant_call.posterior, reference_, reads, phase_region);
                    },
                    [this, &reads] (const auto& somatic_call) {
                        return output_somatic_variant_call(somatic_call.allele, somatic_call.posterior, reference_, reads);
                    },
                    [] (const auto& lhs, const auto& rhs) {
                        return is_before(lhs.allele, rhs.variants.front());
                    });
    
    return result;
}
} // namespace Octopus
