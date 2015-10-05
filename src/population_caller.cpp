//
//  basic_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_caller.hpp"

#include <unordered_map>
#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "common.hpp"
#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "haplotype_tree.hpp"
#include "search_regions.hpp"
#include "vcf_record.hpp"
#include "maths.hpp"
#include "mappable_algorithms.hpp"
#include "variant_utils.hpp"
#include "genotype_model.hpp"
#include "population_genotype_model.hpp"
#include "read_utils.hpp"
#include "string_utils.hpp"
#include "sequence_utils.hpp"
#include "search_regions.hpp"
#include "random_candidate_variant_generator.hpp"

#include <iostream> // TEST

namespace Octopus
{

PopulationVariantCaller::PopulationVariantCaller(ReferenceGenome& reference, CandidateVariantGenerator& candidate_generator,
                                                 RefCall refcalls, double min_posterior, unsigned ploidy)
:
VariantCaller {reference, candidate_generator, refcalls},
ploidy_ {ploidy},
min_posterior_ {min_posterior}
{}

std::string PopulationVariantCaller::do_get_details() const
{
    return "population caller. ploidy = " + std::to_string(ploidy_);
}

GenomicRegion PopulationVariantCaller::get_init_region(const GenomicRegion& region, const ReadMap& reads,
                                                       const std::vector<Variant>& candidates)
{
    return advance_region(get_head(region, 0), reads, candidates, 6, 0);
}

GenomicRegion PopulationVariantCaller::get_next_region(const GenomicRegion& current_region, const ReadMap& reads,
                                                       const std::vector<Variant>& candidates)
{
    return advance_region(current_region, reads, candidates, 5, 0);
}

std::unordered_map<Haplotype, double>
compute_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                             const GenotypeModel::Population::SampleGenotypeProbabilities& genotype_posteriors)
{
    std::unordered_map<Haplotype, double> result {};
    result.reserve(haplotypes.size());
    
    for (const auto& haplotype : haplotypes) {
        for (const auto& genotype_posterior: genotype_posteriors) {
            if (genotype_posterior.first.contains(haplotype)) {
                result[haplotype] += genotype_posterior.second;
            }
        }
    }
    
    return result;
}

double marginalise(const Allele& allele, const GenotypeModel::Population::SampleGenotypeProbabilities& genotype_posteriors)
{
    double result {};
    
    for (const auto& genotype_posterior : genotype_posteriors) {
        if (std::any_of(std::cbegin(genotype_posterior.first), std::cend(genotype_posterior.first),
                        [&allele] (const auto& haplotype) {
                            return haplotype.contains(allele);
                        })) {
                            result += genotype_posterior.second;
                        }
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

auto call_genotype(const GenotypeModel::Population::SampleGenotypeProbabilities& genotype_posteriors)
{
    return *std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                             [] (const auto& lhs, const auto& rhs) {
                                 return lhs.second < rhs.second;
                             });
}

using GenotypeCalls = std::unordered_map<Octopus::SampleIdType, std::pair<Genotype<Allele>, double>>;

std::vector<GenotypeCalls>
call_genotypes(const GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors,
               const std::vector<GenomicRegion>& regions)
{
    std::vector<GenotypeCalls> result(regions.size());
    
    for (const auto sample_genotype_posteriors : genotype_posteriors) {
        const auto& sample_genotype_call = call_genotype(sample_genotype_posteriors.second);
        for (size_t i {}; i < regions.size(); ++i) {
            result[i].emplace(sample_genotype_posteriors.first,
                              std::make_pair(splice<Allele>(sample_genotype_call.first, regions[i]),
                                             sample_genotype_call.second));
        }
    }
    
    return result;
}

using AllelePosteriors = std::unordered_map<Octopus::SampleIdType, std::unordered_map<Allele, double>>;

AllelePosteriors compute_allele_posteriors(const GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors,
                                           const std::vector<Haplotype>& haplotypes,
                                           const std::vector<Allele>& alleles)
{
    AllelePosteriors result {};
    result.reserve(genotype_posteriors.size());
    
    for (const auto sample_genotype_posteriors : genotype_posteriors) {
        auto haplotype_posteriors = compute_haplotype_posteriors(haplotypes, sample_genotype_posteriors.second);
        result.emplace(sample_genotype_posteriors.first,
                       compute_sample_allele_posteriors(sample_genotype_posteriors.second, alleles));
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

double max_posterior(const Allele& allele, const AllelePosteriors& allele_posteriors)
{
    double result {};
    
    for (const auto& sample_allele_posteriors : allele_posteriors) {
        auto p = sample_allele_posteriors.second.at(allele);
        if (p > result) result = p;
    }
    
    return result;
}

double max_posterior(const std::vector<Variant>& variants, const AllelePosteriors& allele_posteriors)
{
    double result {};
    
    for (const auto& sample_allele_posteriors : allele_posteriors) {
        for (const auto& variant : variants) {
            auto p = sample_allele_posteriors.second.at(variant.get_alternative_allele());
            if (p > result) result = p;
        }
    }
    
    return result;
}

template <typename T>
static std::vector<GenomicRegion> get_segment_regions(const std::vector<std::vector<T>>& segments)
{
    std::vector<GenomicRegion> result {};
    result.reserve(segments.size());
    for (const auto& segment : segments) {
        result.push_back(get_encompassing_region(segment));
    }
    return result;
}

std::vector<Variant>
call_segment_variants(const std::vector<Variant>& variants, const AllelePosteriors& allele_posteriors,
                      double min_posterior)
{
    std::vector<Variant> result {};
    result.reserve(variants.size());
    
    std::copy_if(std::cbegin(variants), std::cend(variants), std::back_inserter(result),
                 [&allele_posteriors, min_posterior] (const auto& variant) {
                     return max_posterior(variant.get_alternative_allele(), allele_posteriors) >= min_posterior;
                 });
    
    return result;
}

double get_homozygous_reference_posterior(const Allele& reference_allele, const GenotypeCalls& genotype_calls)
{
    if (genotype_calls.empty()) return 0.0;
    
    double min_posterior {1.0};
    
    for (const auto& sample_genotype_call : genotype_calls) {
        if (is_homozygous_reference(sample_genotype_call.second.first, reference_allele)) {
            if (sample_genotype_call.second.second < min_posterior) {
                min_posterior = sample_genotype_call.second.second;
            }
        } else {
            return 0.0;
        }
    }
    
    return min_posterior;
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

VcfRecord call_variant(const Allele& reference_allele, const std::vector<Variant>& variants,
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
        
        result.add_genotype_field(sample, "FT", "PASS"); // TODO
        result.add_genotype_field(sample, "GP", std::to_string(Maths::probability_to_phred(sample_call.second.second)));
        result.add_genotype_field(sample, "PS", std::to_string(get_begin(phase_region)));
        result.add_genotype_field(sample, "PQ", "60"); // TODO
        result.add_genotype_field(sample, "DP", std::to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", std::to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", std::to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result.build_once();
}
    
VcfRecord call_reference(const Allele& reference_allele, double posterior,
                         ReferenceGenome& reference, const ReadMap& reads)
{
    auto result = VcfRecord::Builder();
    
    auto phred_quality = Maths::probability_to_phred(posterior);
    
    const auto& region = get_region(reference_allele);
    
    result.set_chromosome(get_contig_name(region));
    result.set_position(get_begin(region));
    result.set_ref_allele(reference_allele.get_sequence().front());
    
    result.set_quality(phred_quality);
    
    result.set_filters({"REFCALL"});
    if (size(region) > 1) {
        result.add_info("END", std::to_string(get_end(region)));
    }
    
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
        
        result.add_genotype_field(sample, "FT", "PASS"); // TODO
        result.add_genotype_field(sample, "GP", std::to_string(Maths::probability_to_phred(sample_call.second.second)));
        result.add_genotype_field(sample, "PS", std::to_string(get_begin(phase_region)));
        result.add_genotype_field(sample, "PQ", "60"); // TODO
        result.add_genotype_field(sample, "DP", std::to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", std::to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", std::to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result.build_once();
}

std::vector<Variant> generate_random_variants(const GenomicRegion& region, ReferenceGenome& reference)
{
    RandomCandidateVariantGenerator generator {reference};
    return generator.get_candidates(region);
}

using VariantCalls = std::vector<std::pair<std::vector<Variant>, double>>;
    
VariantCalls
call_segment_variants(const std::vector<std::vector<Variant>>& segments,
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
    
    auto get_variant_call_regions(const VariantCalls& variant_calls)
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
    
    void append_allele(std::vector<Allele>& alleles, const Allele& allele, VariantCaller::RefCall refcalls)
    {
        if (refcalls == VariantCaller::RefCall::Blocked && !alleles.empty() && are_adjacent(alleles.back(), allele)) {
            alleles.back() = Allele {
                get_encompassing(alleles.back(), allele),
                alleles.back().get_sequence() + allele.get_sequence()
            };
        } else {
            alleles.push_back(allele);
        }
    }
    
    std::vector<Allele>
    generate_reference_alleles(const std::vector<Allele>& callable_alleles,
                               const VariantCalls& variant_calls,
                               const std::vector<Variant>& candidates,
                               VariantCaller::RefCall refcalls)
    {
        if (callable_alleles.empty() || refcalls == VariantCaller::RefCall::None) return {};
        
        if (candidates.empty()) return callable_alleles;
        
        auto called_regions = get_variant_call_regions(variant_calls);
        
        auto allele_itr        = std::cbegin(callable_alleles);
        auto allele_end_itr    = std::cend(callable_alleles);
        auto called_itr        = std::cbegin(called_regions);
        auto called_end_itr    = std::cend(called_regions);
        auto candidate_itr     = std::cbegin(candidates);
        auto candidate_end_itr = std::cend(candidates);
        
        std::vector<Allele> result {};
        result.reserve(callable_alleles.size());
        
        while (allele_itr != allele_end_itr) {
            if (candidate_itr == candidate_end_itr) {
                append_allele(result, *allele_itr, refcalls);
                std::copy(std::next(allele_itr), allele_end_itr, std::back_inserter(result));
                break;
            }
            
            if (called_itr == called_end_itr) {
                append_allele(result, *allele_itr, refcalls);
                if (begins_before(*allele_itr, *candidate_itr)) {
                    std::advance(allele_itr, 1);
                } else {
                    std::advance(allele_itr, 2); // as next allele is alt from candidate
                    std::advance(candidate_itr, 1);
                }
            } else {
                if (begins_before(*called_itr, *candidate_itr)) {
                    if (!overlaps(*allele_itr, *called_itr)) {
                        append_allele(result, *allele_itr, refcalls);
                        std::advance(allele_itr, 1);
                    } else {
                        if (begins_before(*allele_itr, *called_itr)) { // when variant has been padded
                            append_allele(result, splice(*allele_itr, get_left_overhang(*allele_itr, *called_itr)), refcalls);
                        }
                        
                        // skip contained alleles and candidates as they include called variants
                        allele_itr    = contained_range(allele_itr, allele_end_itr, *called_itr).end().base();
                        candidate_itr = contained_range(candidate_itr, candidate_end_itr, *called_itr).end().base();
                        
                        std::advance(called_itr, 1);
                    }
                } else {
                    append_allele(result, *allele_itr, refcalls);
                    if (begins_before(*allele_itr, *candidate_itr)) {
                        std::advance(allele_itr, 1);
                    } else {
                        std::advance(allele_itr, 2); // as next allele is alt from candidate
                        std::advance(candidate_itr, 1);
                    }
                }
            }
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    std::vector<std::pair<Allele, double>>
    make_refcalls(const GenotypeModel::Population::GenotypeProbabilities& genotype_posteriors,
                  const std::vector<Allele>& reference_alleles, double min_posterior)
    {
        std::vector<std::pair<Allele, double>> result {};
        result.reserve(reference_alleles.size());
        
        // TODO: do we need to call here? could just check the posteriors for homo-ref
        auto genotype_calls = call_genotypes(genotype_posteriors, get_regions(reference_alleles));
        
        auto ref_allele_itr = std::cbegin(reference_alleles);
        
        for (const auto& genotype_call : genotype_calls) {
            auto posterior = get_homozygous_reference_posterior(*ref_allele_itr, genotype_call);
            if (posterior >= min_posterior) {
                result.emplace_back(*ref_allele_itr, posterior);
            }
            std::advance(ref_allele_itr, 1);
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    template <typename InputIt1, typename InputIt2, typename InputIt3, typename OutputIt,
              typename BinaryOperation, typename UnaryOperation, typename Compare>
    OutputIt merge_transform(InputIt1 first1, InputIt1 last1,
                             InputIt2 first2, InputIt2 last2,
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
        return std::transform(first3, first3, d_first, unary_op);
    }
    
std::vector<VcfRecord>
PopulationVariantCaller::call_variants(const GenomicRegion& region, const std::vector<Variant>& candidates,
                                       const ReadMap& reads)
{
    std::vector<VcfRecord> result {};
    
    if (empty(region)) return result;
    
    HaplotypeTree tree {reference_};
    
    if (candidates.empty()) {
        if (refcalls_requested()) {
            extend_tree(generate_random_variants(region, reference_), tree);
        } else {
            return result;
        }
    } else {
        extend_tree(candidates, tree);
    }
    
    auto haplotypes = tree.get_haplotypes(region);
    
    std::cout << "there are " << haplotypes.size() << " haplotypes" << std::endl;
    
    GenotypeModel::Population genotype_model {ploidy_};
    
    auto genotype_posteriors = genotype_model.evaluate(haplotypes, reads, reference_).genotype_posteriors;
    
    auto alleles = generate_callable_alleles(region, candidates, refcalls_, reference_);
    
    auto allele_posteriors = compute_allele_posteriors(genotype_posteriors, haplotypes, alleles);
    
    auto variant_calls = call_segment_variants(segment_overlapped(candidates), allele_posteriors, min_posterior_);
    
    parsimonise_variant_calls(variant_calls, reference_);
    
    auto called_regions = get_variant_call_regions(variant_calls);
    
    auto candidate_refcall_alleles = generate_reference_alleles(alleles, variant_calls, candidates, refcalls_);
    
    auto variant_genotype_calls = call_genotypes(genotype_posteriors, called_regions);
    
    auto refcalls = make_refcalls(genotype_posteriors, candidate_refcall_alleles, min_posterior_);
    
    auto phase_region = get_encompassing(called_regions.front(), called_regions.back());
    
    result.reserve(variant_calls.size() + refcalls.size());
    
    merge_transform(std::cbegin(variant_calls), std::cend(variant_calls),
                    std::cbegin(variant_genotype_calls), std::cend(variant_genotype_calls),
                    std::cbegin(refcalls), std::cend(refcalls), std::back_inserter(result),
                    [this, &reads, &phase_region] (const auto& variant_call, const auto& genotype_call) {
                        return call_variant(variant_call.first.front().get_reference_allele(),
                                            variant_call.first, genotype_call,
                                            variant_call.second, reference_, reads, phase_region);
                    },
                    [this, &reads] (const auto& refcall) {
                        return call_reference(refcall.first, refcall.second, reference_, reads);
                    },
                    [] (const auto& lhs, const auto& rhs) {
                        return get_region(lhs.first) < get_region(rhs.first.front());
                    });
    
    return result;
}
    
} // namespace Octopus
