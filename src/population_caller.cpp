//
//  basic_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_caller.hpp"

#include <unordered_map>
#include <unordered_set>
#include <unordered_map>
#include <utility>   // std::pair, std::make_pair, std::move, std::forward
#include <algorithm> // std::max_element, std::copy_if, std::copy, std::transform
#include <iterator>  // std::begin, std::end, std::cbegin, std::cend, std::back_inserter, std::advance, std::next

#include "common.hpp"
#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "vcf_record.hpp"
#include "maths.hpp"
#include "mappable_algorithms.hpp"
#include "variant_utils.hpp"
#include "genotype_model.hpp"
#include "population_genotype_model.hpp"
#include "read_utils.hpp"
#include "string_utils.hpp"
#include "sequence_utils.hpp"
#include "random_candidate_variant_generator.hpp"

#include <iostream> // TEST

namespace Octopus
{

// public methods

PopulationVariantCaller::PopulationVariantCaller(ReferenceGenome& reference, CandidateVariantGenerator& candidate_generator,
                                                 RefCallType refcall_type, double min_variant_posterior,
                                                 double min_refcall_posterior, unsigned ploidy)
:
VariantCaller {reference, candidate_generator, refcall_type},
phaser_ {reference, 1000, 0},
ploidy_ {ploidy},
min_variant_posterior_ {min_variant_posterior},
min_refcall_posterior_ {min_refcall_posterior}
{}

std::string PopulationVariantCaller::do_get_details() const
{
    return "population caller. ploidy = " + std::to_string(ploidy_);
}

// non member methods

using GenotypeCalls = std::unordered_map<Octopus::SampleIdType, std::pair<Genotype<Allele>, double>>;
using VariantCalls = std::vector<std::pair<std::vector<Variant>, double>>;
using AllelePosteriors = std::unordered_map<Octopus::SampleIdType, std::unordered_map<Allele, double>>;

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

    namespace debug {
        void print_top_genotype_posteriors(const GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors, const size_t n = 5);
        void print_top_allele_posteriors(const AllelePosteriors& allele_posteriors, const size_t n = 10);
    } // namespace debug
    
//std::unordered_map<Haplotype, double>
//compute_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
//                             const GenotypeModel::Population::SampleGenotypeProbabilities& genotype_posteriors)
//{
//    std::unordered_map<Haplotype, double> result {};
//    result.reserve(haplotypes.size());
//    
//    for (const auto& haplotype : haplotypes) {
//        for (const auto& genotype_posterior: genotype_posteriors) {
//            if (genotype_posterior.first.contains(haplotype)) {
//                result[haplotype] += genotype_posterior.second;
//            }
//        }
//    }
//    
//    return result;
//}

void remove_low_posterior_genotypes(GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors,
                                    double min_posterior)
{
    for (auto& sample_genotype_posteriors : genotype_posteriors) {
        for (auto itr = std::begin(sample_genotype_posteriors.second); itr != std::end(sample_genotype_posteriors.second); ) {
            if (itr->second < min_posterior) {
                itr = sample_genotype_posteriors.second.erase(itr);
            } else {
                ++itr;
            }
        }
    }
}
    
double marginalise(const Allele& allele, const GenotypeModel::Population::SampleGenotypeProbabilities& genotype_posteriors)
{
    double result {0.0};
    
    for (const auto& genotype_posterior : genotype_posteriors) {
        if (contains(genotype_posterior.first, allele)) result += genotype_posterior.second;
    }
    
    return result;
}

std::unordered_map<Allele, double>
compute_sample_allele_posteriors(const GenotypeModel::Population::SampleGenotypeProbabilities& genotype_posteriors,
                                 const std::vector<Allele>& alleles)
{
    std::unordered_map<Allele, double> result {};
    result.reserve(alleles.size());
    
    for (const auto& allele : alleles) {
        result.emplace(allele, marginalise(allele, genotype_posteriors));
    }
    
    return result;
}

AllelePosteriors
compute_allele_posteriors(const GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors,
                          const std::vector<Allele>& alleles)
{
    AllelePosteriors result {};
    result.reserve(genotype_posteriors.size());
    
    for (const auto sample_genotype_posteriors : genotype_posteriors) {
        result.emplace(sample_genotype_posteriors.first,
                       compute_sample_allele_posteriors(sample_genotype_posteriors.second, alleles));
    }
    
    return result;
}

std::unordered_map<Genotype<Allele>, double>
marginalise(const GenomicRegion& region, const GenotypeModel::Population::SampleGenotypeProbabilities& genotype_posteriors)
{
    std::unordered_map<Genotype<Allele>, double> result {};
    result.reserve(genotype_posteriors.size());
    
    for (const auto& genotype_posterior : genotype_posteriors) {
        result[splice<Allele>(genotype_posterior.first, region)] += genotype_posterior.second;
    }
    
    return result;
}

template <typename Map>
auto call_genotype(const Map& genotype_posteriors)
{
    return *std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                             [] (const auto& lhs, const auto& rhs) {
                                 return lhs.second < rhs.second;
                             });
}

std::vector<GenotypeCalls>
call_genotypes(const GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors,
               const std::vector<GenomicRegion>& regions)
{
    std::vector<GenotypeCalls> result(regions.size());
    
    for (const auto sample_genotype_posteriors : genotype_posteriors) {
        // A better way would be to splice out the region and marginalise over this spliced region
        const auto& sample_genotype_call = call_genotype(sample_genotype_posteriors.second);
        for (size_t i {}; i < regions.size(); ++i) {
            result[i].emplace(sample_genotype_posteriors.first,
                              std::make_pair(splice<Allele>(sample_genotype_call.first, regions[i]),
                                             sample_genotype_call.second));
        }
    }
    
    return result;
}
    
//    std::vector<GenotypeCalls>
//    call_genotypes(const GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors,
//                   const std::vector<GenomicRegion>& regions)
//    {
//        std::vector<GenotypeCalls> result(regions.size());
//        
//        for (const auto sample_genotype_posteriors : genotype_posteriors) {
//            for (size_t i {}; i < regions.size(); ++i) {
//                result[i].emplace(sample_genotype_posteriors.first,
//                                  call_genotype(marginalise(regions[i], sample_genotype_posteriors.second)));
//            }
//        }
//        
//        return result;
//    }

double max_posterior(const Allele& allele, const AllelePosteriors& allele_posteriors)
{
    double result {};
    
    for (const auto& sample_allele_posteriors : allele_posteriors) {
        auto curr = sample_allele_posteriors.second.at(allele);
        if (curr > result) result = curr;
    }
    
    return result;
}

double max_posterior(const std::vector<Variant>& variants, const AllelePosteriors& allele_posteriors)
{
    double result {};
    
    for (const auto& sample_allele_posteriors : allele_posteriors) {
        for (const auto& variant : variants) {
            auto curr = sample_allele_posteriors.second.at(variant.get_alternative_allele());
            if (curr > result) result = curr;
        }
    }
    
    return result;
}

std::vector<Variant> call_segment_variants(const std::vector<Variant>& variants,
                                           const AllelePosteriors& allele_posteriors,
                                           double min_posterior)
{
    std::vector<Variant> result {};
    result.reserve(variants.size());
    
    std::copy_if(std::cbegin(variants), std::cend(variants), std::back_inserter(result),
                 [&allele_posteriors, min_posterior] (const auto& variant) {
                     return max_posterior(variant.get_alternative_allele(), allele_posteriors) >= min_posterior;
                 });
    
    result.shrink_to_fit();
    
    return result;
}

VariantCalls call_segment_variants(const std::vector<std::vector<Variant>>& segments,
                                   const AllelePosteriors& allele_posteriors, double min_posterior)
{
    VariantCalls result {};
    result.reserve(segments.size());
    
    for (const auto& segment : segments) {
        auto calls = call_segment_variants(segment, allele_posteriors, min_posterior);
        if (!calls.empty()) {
            result.emplace_back(std::move(calls), max_posterior(calls, allele_posteriors));
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

std::vector<Variant> generate_random_variants(const GenomicRegion& region, ReferenceGenome& reference)
{
    RandomCandidateVariantGenerator generator {reference};
    return generator.get_candidates(region);
}

auto get_regions(const VariantCalls& variant_calls)
{
    std::vector<GenomicRegion> result {};
    result.reserve(variant_calls.size());
    
    for (const auto& segment_calls : variant_calls) {
        result.emplace_back(get_encompassing_region(segment_calls.first));
    }
    
    return result;
}

void parsimonise_variant_calls(VariantCalls& variant_calls, ReferenceGenome& reference)
{
    for (auto& segment_calls : variant_calls) {
        segment_calls.first = parsimonise_together(segment_calls.first, reference);
    }
}

double compute_homozygous_reference_posterior(const Allele& reference_allele, const GenotypeCalls& genotype_calls)
{
    if (genotype_calls.empty()) return 0.0;
    
    double result {1.0};
    
    for (const auto& sample_genotype_call : genotype_calls) {
        if (is_homozygous(sample_genotype_call.second.first, reference_allele)) {
            auto curr = sample_genotype_call.second.second;
            if (curr < result) {
                result = curr;
            }
        } else {
            return 0.0;
        }
    }
    
    return result;
}

std::vector<RefCall>
call_reference(const GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors,
               const std::vector<Allele>& reference_alleles, double min_posterior)
{
    std::vector<RefCall> result {};
    
    if (reference_alleles.empty()) return result;
    
    result.reserve(reference_alleles.size());
    
    auto genotype_calls = call_genotypes(genotype_posteriors, get_regions(reference_alleles));
    
    auto ref_allele_itr = std::cbegin(reference_alleles);
    
    for (const auto& genotype_call : genotype_calls) {
        auto refcall_posterior = compute_homozygous_reference_posterior(*ref_allele_itr, genotype_call);
        
        //std::cout << refcall_posterior << std::endl;
        
        if (refcall_posterior >= min_posterior) {
            std::vector<std::pair<SampleIdType, double>> sample_posteriors {};
            sample_posteriors.reserve(genotype_call.size());
            
            for (const auto& sample_posterior : genotype_call) {
                sample_posteriors.emplace_back(sample_posterior.first, sample_posterior.second.second);
            }
            
            result.emplace_back(*ref_allele_itr, refcall_posterior, std::move(sample_posteriors));
        }
        
        std::advance(ref_allele_itr, 1);
    }
    
    result.shrink_to_fit();
    
    return result;
}

double marginalise_reference_genotype(const Allele& reference_allele,
                                      const GenotypeModel::Population::SampleGenotypeProbabilities& sample_genotype_posteriors)
{
    double result {};
    
    for (const auto& genotype_posterior : sample_genotype_posteriors) {
        if (is_homozygous(genotype_posterior.first, reference_allele)) {
            result += genotype_posterior.second;
        }
    }
    
    return result;
}

std::vector<RefCall>
call_reference_better(const GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors,
                      const std::vector<Allele>& reference_alleles, double min_posterior)
{
    std::vector<RefCall> result {};
    
    if (reference_alleles.empty()) return result;
    
    result.reserve(reference_alleles.size());
    
    for (const auto& reference_allele : reference_alleles) {
        double min_sample_posteior {1.0};
        
        std::vector<std::pair<SampleIdType, double>> sample_posteriors {};
        sample_posteriors.reserve(genotype_posteriors.size());
        
        for (const auto& sample_genotype_posteriors : genotype_posteriors) {
            auto sample_posterior = marginalise_reference_genotype(reference_allele, sample_genotype_posteriors.second);
            
            if (sample_posterior < min_posterior) {
                min_sample_posteior = sample_posterior;
                break; // to avoid computing the rest
            } else if (sample_posterior < min_sample_posteior) {
                min_sample_posteior = sample_posterior;
            }
            
            sample_posteriors.emplace_back(sample_genotype_posteriors.first, sample_posterior);
        }
        
        if (min_sample_posteior >= min_posterior) {
            result.emplace_back(reference_allele, min_sample_posteior, std::move(sample_posteriors));
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

unsigned count_alleles(const GenotypeCalls& genotype_calls)
{
    if (genotype_calls.empty()) return 0;
    
    std::unordered_set<Allele> unique_alleles {};
    
    for (const auto& sample_call : genotype_calls) {
        for (const auto& allele : sample_call.second.first) {
            unique_alleles.emplace(allele);
        }
    }
    
    return static_cast<unsigned>(unique_alleles.size());
}

std::vector<unsigned> count_alleles(const std::vector<Variant>& variants, const GenotypeCalls& genotype_calls)
{
    std::unordered_map<Allele, unsigned> allele_counts {};
    allele_counts.reserve(variants.size());
    
    for (const auto& sample_call : genotype_calls) {
        for (const auto& allele : sample_call.second.first) {
            ++allele_counts[allele];
        }
    }
    
    std::vector<unsigned> result {};
    result.reserve(variants.size());
    
    for (const auto& variant : variants) {
        const auto& alt_allele = variant.get_alternative_allele();
        result.push_back(allele_counts[alt_allele]);
    }
    
    return result;
}

std::vector<VcfRecord::SequenceType> to_vcf_genotype(const Genotype<Allele>& genotype)
{
    std::vector<VcfRecord::SequenceType> result {};
    result.reserve(genotype.ploidy());
    for (const auto& allele : genotype) result.push_back(allele.get_sequence());
    return result;
}

VcfRecord output_variant_call(const Allele& reference_allele, const std::vector<Variant>& variants,
                              const GenotypeCalls& genotype_calls, double posterior,
                              ReferenceGenome& reference, const ReadMap& reads,
                              const GenomicRegion& phase_region)
{
    auto result = VcfRecord::Builder();
    
    auto phred_quality = Maths::probability_to_phred(posterior);
    
    const auto& region = get_region(reference_allele);
    
    result.set_chromosome(get_contig_name(region));
    result.set_position(get_begin(region));
    result.set_ref_allele(reference_allele.get_sequence());
    result.set_alt_alleles(get_alt_allele_sequences(variants));
    result.set_quality(phred_quality);
    
    //result.set_filters({"PASS"}); // TODO
    
    result.add_info("AC", to_strings(count_alleles(variants, genotype_calls)));
    result.add_info("AN", std::to_string(count_alleles(genotype_calls)));
    result.add_info("NS", std::to_string(count_samples_with_coverage(reads, region)));
    result.add_info("DP", std::to_string(sum_max_coverages(reads, region)));
    result.add_info("SB", to_string(strand_bias(reads, region), 2));
    result.add_info("BQ", std::to_string(static_cast<unsigned>(rmq_base_quality(reads, region))));
    result.add_info("MQ", std::to_string(static_cast<unsigned>(rmq_mapping_quality(reads, region))));
    result.add_info("MQ0", std::to_string(count_mapq_zero(reads, region)));
    
    result.set_format({"GT", "FT", "GP", "PS", "PQ", "DP", "BQ", "MQ"});
    
    for (const auto& sample_call : genotype_calls) {
        const auto& sample = sample_call.first;
        result.add_genotype(sample, to_vcf_genotype(sample_call.second.first), VcfRecord::Builder::Phasing::Phased);
        
        result.add_genotype_field(sample, "FT", "."); // TODO
        result.add_genotype_field(sample, "GP", std::to_string(Maths::probability_to_phred(sample_call.second.second)));
        result.add_genotype_field(sample, "PS", std::to_string(get_begin(phase_region) + 1));
        result.add_genotype_field(sample, "PQ", "60"); // TODO
        result.add_genotype_field(sample, "DP", std::to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", std::to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", std::to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result.build_once();
}

VcfRecord output_reference_call(RefCall call, ReferenceGenome& reference, const ReadMap& reads, unsigned ploidy)
{
    auto result = VcfRecord::Builder();
    
    auto phred_quality = Maths::probability_to_phred(call.posterior);
    
    const auto& region = get_region(call.reference_allele);
    
    result.set_chromosome(get_contig_name(region));
    result.set_position(get_begin(region));
    result.set_ref_allele(call.reference_allele.get_sequence().front());
    
    result.set_quality(phred_quality);
    
    result.set_filters({"REFCALL"});
    if (size(region) > 1) {
        result.add_info("END", std::to_string(get_end(region))); // - 1 as VCF uses closed intervals
    }
    
    result.add_info("NS", std::to_string(count_samples_with_coverage(reads, region)));
    result.add_info("DP", std::to_string(sum_max_coverages(reads, region)));
    result.add_info("SB", to_string(strand_bias(reads, region), 2));
    result.add_info("BQ", std::to_string(static_cast<unsigned>(rmq_base_quality(reads, region))));
    result.add_info("MQ", std::to_string(static_cast<unsigned>(rmq_mapping_quality(reads, region))));
    result.add_info("MQ0", std::to_string(count_mapq_zero(reads, region)));
    
    result.set_format({"GT", "GP", "DP", "BQ", "MQ"});
    
    for (const auto& sample_posteior : call.sample_posteriors) {
        const auto& sample = sample_posteior.first;
        result.add_homozygous_ref_genotype(sample, ploidy);
        result.add_genotype_field(sample, "GP", std::to_string(Maths::probability_to_phred(sample_posteior.second)));
        result.add_genotype_field(sample, "DP", std::to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", std::to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", std::to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result.build_once();
}

template <typename InputIt1, typename InputIt2, typename InputIt3, typename OutputIt,
typename BinaryOperation, typename UnaryOperation, typename Compare>
OutputIt merge_transform(InputIt1 first1, InputIt1 last1, InputIt2 first2,
                         InputIt3 first3, InputIt3 last3,
                         OutputIt d_first,
                         BinaryOperation binary_op, UnaryOperation unary_op,
                         Compare cmp)
{
    for (; first1 != last1; ++d_first) {
        if (first3 == last3) {
            return std::transform(first1, last1, first2, d_first, binary_op);
        }
        if (cmp(*first3, *first1)) {
            *d_first = unary_op(*first3);
            ++first3;
        } else {
            *d_first = binary_op(*first1, *first2);
            ++first1;
            ++first2;
        }
    }
    return std::transform(first3, last3, d_first, unary_op);
}

std::vector<VcfRecord>
PopulationVariantCaller::call_variants(const GenomicRegion& region, const std::vector<Variant>& candidates,
                                       const ReadMap& reads)
{
    std::vector<VcfRecord> result {};
    
    if (empty(region)) return result;
    
    phaser_.setup(candidates, reads);
    
    while (!phaser_.expended_candidates()) {
        auto haplotypes = phaser_.get_haplotypes();
        
        //    if (candidates.empty()) {
        //        if (refcalls_requested()) {
        //            extend_tree(generate_random_variants(region, reference_), tree);
        //        } else {
        //            return result;
        //        }
        //    } else {
        //        extend_tree(candidates, tree);
        //    }
        
        std::cout << "there are " << haplotypes.size() << " unique haplotypes" << std::endl;
        
        auto haplotype_region = get_region(haplotypes.front());
        
        std::cout << "haplotype region is " << haplotype_region << std::endl;
        
        auto haplotype_region_reads = copy_overlapped(reads, haplotype_region);
        
        std::cout << "there are " << count_reads(haplotype_region_reads) << " reads in haplotype region" << std::endl;
        
        GenotypeModel::Population genotype_model {ploidy_};
        
        auto genotype_posteriors = genotype_model.evaluate(haplotypes, haplotype_region_reads, reference_).genotype_posteriors;
        
        auto phased_gps = phaser_.phase(haplotypes, genotype_posteriors, reads);
        
        //remove_low_posterior_genotypes(genotype_posteriors, 0.0000000001);
        
        auto alleles = generate_callable_alleles(region, candidates, refcall_type_, reference_);
        
        auto allele_posteriors = compute_allele_posteriors(genotype_posteriors, alleles);
        
        auto variant_calls = call_segment_variants(segment_overlapped(candidates), allele_posteriors, min_variant_posterior_);
        
        parsimonise_variant_calls(variant_calls, reference_);
        
        auto called_regions = get_regions(variant_calls);
        
        auto variant_genotype_calls = call_genotypes(genotype_posteriors, called_regions);
        
        //debug::print_top_genotype_posteriors(genotype_posteriors);
        
        //debug::print_top_allele_posteriors(allele_posteriors);
        
//        std::cout << "variants" << std::endl;
//        for (auto v : variant_calls) {
//            std::cout << v.first.front() << std::endl;
//        }
//        
//        std::cout << "variant regions" << std::endl;
//        for (auto r : called_regions) {
//            std::cout << r << std::endl;
//        }
//        
//        std::cout << "variant genotypes" << std::endl;
//        for (auto g : variant_genotype_calls) {
//            for (auto sc : g) {
//                std::cout << sc.first << std::endl;
//                std::cout << sc.second.first << std::endl;
//            }
//        }
        
        auto candidate_ref_alleles = generate_candidate_reference_alleles(alleles, called_regions, candidates, refcall_type_);
        
        //    std::cout << "candidate refcall alleles are" << std::endl;
        //    for (const auto& allele : candidate_ref_alleles) {
        //        std::cout << allele << std::endl;
        //    }
        
        //auto refcalls = call_reference(genotype_posteriors, candidate_ref_alleles, min_refcall_posterior_);
        auto refcalls = call_reference_better(genotype_posteriors, candidate_ref_alleles, min_refcall_posterior_);
        
        auto phase_region = (called_regions.empty()) ? get_head(region) : get_encompassing(called_regions.front(), called_regions.back());
        
        result.reserve(variant_calls.size() + refcalls.size());
        
        merge_transform(std::cbegin(variant_calls), std::cend(variant_calls), std::cbegin(variant_genotype_calls),
                        std::cbegin(refcalls), std::cend(refcalls), std::back_inserter(result),
                        [this, &reads, &phase_region] (const auto& variant_call, const auto& genotype_call) {
                            return output_variant_call(variant_call.first.front().get_reference_allele(),
                                                       variant_call.first, genotype_call,
                                                       variant_call.second, reference_, reads, phase_region);
                        },
                        [this, &reads] (const auto& refcall) {
                            return output_reference_call(refcall, reference_, reads, ploidy_);
                        },
                        [] (const auto& lhs, const auto& rhs) {
                            return is_before(lhs.reference_allele, rhs.first.front());
                        });
    }
    
    return result;
}
    
    namespace debug {
        
        void print_top_genotype_posteriors(const GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors, const size_t n)
        {
            for (const auto& sample_posteriors : genotype_posteriors) {
                auto m = std::min(n, sample_posteriors.second.size());
                std::cout << "printing top " << m << " genotype posteriors for sample " << sample_posteriors.first << std::endl;
                
                std::vector<std::pair<Genotype<Haplotype>, double>> v {};
                v.reserve(sample_posteriors.second.size());
                
                std::copy(std::cbegin(sample_posteriors.second), std::cend(sample_posteriors.second), std::back_inserter(v));
                
                std::sort(std::begin(v), std::end(v), [] (const auto& lhs, const auto& rhs) {
                    return lhs.second > rhs.second;
                });
                
                for (size_t i {}; i < m; ++i) {
                    std::cout << "\t* ";
                    print_variant_alleles(v[i].first);
                    std::cout << " " << v[i].second << std::endl;
                }
            }
        }
        
        void print_top_allele_posteriors(const AllelePosteriors& allele_posteriors, const size_t n)
        {
            for (const auto& sample_posteriors : allele_posteriors) {
                auto m = std::min(n, sample_posteriors.second.size());
                std::cout << "printing top " << m << " allele posteriors for sample " << sample_posteriors.first << std::endl;
                
                std::vector<std::pair<Allele, double>> v {};
                v.reserve(sample_posteriors.second.size());
                
                std::copy(std::cbegin(sample_posteriors.second), std::cend(sample_posteriors.second), std::back_inserter(v));
                
                std::sort(std::begin(v), std::end(v), [] (const auto& lhs, const auto& rhs) {
                    return lhs.second > rhs.second;
                });
                
                for (size_t i {}; i < m; ++i) {
                    std::cout << "\t* ";
                    std::cout << v[i].first << " " << v[i].second << std::endl;
                }
            }
        }
    } // namespace debug
    
} // namespace Octopus
