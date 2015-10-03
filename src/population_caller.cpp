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

#include "mappable_algorithms.hpp"
#include "variant_utils.hpp"
#include "genotype_model.hpp"
#include "population_genotype_model.hpp"
#include "read_utils.hpp"

#include "sequence_utils.hpp"
#include "search_regions.hpp"

#include <iostream> // TEST

namespace Octopus
{

PopulationVariantCaller::PopulationVariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                                                 ReadFilter read_filter, ReadTransform read_transform,
                                                 CandidateVariantGenerator& candidate_generator,
                                                 unsigned ploidy)
:
VariantCaller {reference, read_manager, read_filter, read_transform, candidate_generator},
ploidy_ {ploidy}
{}

GenomicRegion PopulationVariantCaller::get_init_region(const GenomicRegion& region)
{
    return region;
}

GenomicRegion PopulationVariantCaller::get_next_region(const GenomicRegion& current_region)
{
    return GenomicRegion {"TEST", 0, 0};
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
               const std::vector<GenomicRegion>& segments)
{
    std::vector<GenotypeCalls> result(segments.size());
    
    for (const auto sample_genotype_posteriors : genotype_posteriors) {
        const auto& sample_genotype_call = call_genotype(sample_genotype_posteriors.second);
        for (size_t i {}; i < segments.size(); ++i) {
            result[i].emplace(sample_genotype_posteriors.first,
                              std::make_pair(splice<Allele>(sample_genotype_call.first, segments[i]),
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

std::vector<VcfRecord::SequenceType> get_alt_allele_sequences(const std::vector<Variant>& segment)
{
    std::vector<VcfRecord::SequenceType> result {};
    result.reserve(segment.size());
    std::transform(std::cbegin(segment), std::cend(segment), std::back_inserter(result),
                   [] (const auto& variant) { return variant.get_alternative_allele_sequence(); });
    return result;
}

static std::vector<GenomicRegion> get_segment_regions(const std::vector<std::vector<Allele>>& segments)
{
    std::vector<GenomicRegion> result {};
    result.reserve(segments.size());
    for (const auto& segment : segments) result.push_back(segment.front().get_region());
    return result;
}

std::vector<Variant> call_segment_variants(const std::vector<Allele>& segment,
                                           const AllelePosteriors& allele_posteriors, double min_posterior)
{
    std::vector<Variant> result {};
    
    const auto& ref_allele = segment.front();
    
    std::for_each(std::next(std::cbegin(segment)), std::cend(segment),
                  [&ref_allele, &allele_posteriors, min_posterior, &result] (const auto& allele) {
                      if (max_posterior(allele, allele_posteriors) >= min_posterior) {
                          result.emplace_back(ref_allele, allele);
                      }
                  });
    
    return result;
}

double get_homozygous_reference_posterior(const Allele& reference, const GenotypeCalls& genotype_calls)
{
    if (genotype_calls.empty()) return 0.0;
    
    double min_posterior {1.0};
    
    for (const auto& sample_genotype_call : genotype_calls) {
        if (is_homozygous_reference(sample_genotype_call.second.first, reference)) {
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

template <typename T>
std::string to_string(const T val, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << val;
    return out.str();
}

template <typename T>
std::vector<std::string> to_string(const std::vector<T>& values)
{
    std::vector<std::string> result {};
    result.reserve(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::back_inserter(result),
                   [] (auto value) { return std::to_string(value); });
    return result;
}

VcfRecord call_segment(const Allele& reference_allele, const std::vector<Variant>& variants,
                       const GenotypeCalls& genotype_calls, unsigned phred_quality,
                       ReferenceGenome& reference, const ReadMap& reads,
                       const GenomicRegion& phase_region)
{
    auto result = VcfRecord::Builder();
    
    const auto& region = get_region(reference_allele);
    
    result.set_chromosome(get_contig_name(region));
    
    //auto parsimonious_variants = make_parsimonious(variants, reference);
    
    result.set_position(get_begin(region));
    result.set_ref_allele(reference_allele.get_sequence());
    result.set_alt_alleles(get_alt_allele_sequences(variants));
    
    result.set_quality(phred_quality);
    
    result.add_info("AC", to_string(count_alleles(variants, genotype_calls)));
    result.add_info("AN", std::to_string(count_alleles(genotype_calls)));
    result.add_info("NS", std::to_string(count_samples_with_coverage(reads, region)));
    result.add_info("DP", std::to_string(sum_max_coverages(reads, region)));
    result.add_info("SB", to_string(strand_bias(reads, region), 2));
    result.add_info("BQ", std::to_string(static_cast<unsigned>(rmq_base_quality(reads, region))));
    result.add_info("MQ", std::to_string(static_cast<unsigned>(rmq_mapping_quality(reads, region))));
    result.add_info("MQ0", std::to_string(count_mapq_zero(reads, region)));
    
    result.set_format({"GT", "FT", "GP", "PS", "PQ", "DP", "BQ", "MQ"});
    
    if (variants.empty()) {
        result.set_filters({"REFCALL"});
    }
    
    for (const auto& sample_call : genotype_calls) {
        const auto& sample = sample_call.first;
        result.add_genotype(sample, to_vcf_genotype(sample_call.second.first), VcfRecord::Builder::Phasing::Phased);
        
        result.add_genotype_field(sample, "FT", "PASS"); // TODO
        result.add_genotype_field(sample, "GP", std::to_string(to_phred_quality(sample_call.second.second)));
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
    auto positions = decompose(region);
    
    std::vector<Variant> result {};
    
    auto position = positions[positions.size() / 2];
    
    auto reference_allele = get_reference_allele(position, reference);
    
    Allele mutation {position, reverse_complement(reference_allele.get_sequence())};
    
    result.emplace_back(reference_allele, mutation);
    
    return result;
}

std::vector<VcfRecord>
PopulationVariantCaller::call_variants(const GenomicRegion& region, const std::vector<Variant>& candidates,
                                       const ReadMap& reads)
{
    std::vector<VcfRecord> result {};
    
    if (empty(region)) return result;
    
    HaplotypeTree tree {reference_};
    
    if (candidates.empty()) {
        if (make_ref_calls_) {
            extend_tree(generate_random_variants(region, reference_), tree);
        } else {
            return result;
        }
    } else {
        extend_tree(candidates, tree);
    }
    
    auto haplotypes = tree.get_haplotypes(region);
    
    GenotypeModel::Population genotype_model {ploidy_};
    
    auto genotype_posteriors = genotype_model.evaluate(haplotypes, reads, reference_).genotype_posteriors;
    
    auto alleles = decompose(candidates);
    
    auto allele_posteriors = compute_allele_posteriors(genotype_posteriors, haplotypes, alleles);
    
    auto segments = segment(alleles);
    auto regions  = get_segment_regions(segments);
    
    auto genotype_calls     = call_genotypes(genotype_posteriors, regions);
    auto genotype_calls_itr = std::cbegin(genotype_calls);
    
    auto phase_region = get_encompassing(regions.front(), regions.back());
    
    for (const auto& segment : segments) {
        auto variants = call_segment_variants(segment, allele_posteriors, min_posterior_);
        
        const auto& reference_allele = segment.front();
        
        if (variants.empty() && make_ref_calls_) {
            auto homozygous_reference_posterior = get_homozygous_reference_posterior(reference_allele, *genotype_calls_itr);
            
            if (homozygous_reference_posterior >= min_posterior_) {
                auto segment_quality = to_phred_quality(homozygous_reference_posterior);
                result.emplace_back(call_segment(reference_allele, variants, *genotype_calls_itr,
                                                 segment_quality, reference_, reads, phase_region));
            }
        } else {
            auto segment_quality = to_phred_quality(max_posterior(variants, allele_posteriors));
            result.emplace_back(call_segment(reference_allele, variants, *genotype_calls_itr,
                                             segment_quality, reference_, reads, phase_region));
        }
        
        ++genotype_calls_itr;
    }
    
    return result;
}
    
} // namespace Octopus
