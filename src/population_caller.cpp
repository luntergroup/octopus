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

VcfRecord call_segment(const Allele& reference_allele, const std::vector<Variant>& variants,
                       const GenotypeCalls& genotype_calls, unsigned phred_quality,
                       ReferenceGenome& reference)
{
    auto result = VcfRecord::Builder();
    
    result.set_chromosome(get_contig_name(reference_allele));
    
    //auto parsimonious_variants = make_parsimonious(variants, reference);
    
    result.set_position(get_begin(reference_allele));
    result.set_ref_allele(reference_allele.get_sequence());
    result.set_alt_alleles(get_alt_allele_sequences(variants));
    
    result.set_quality(phred_quality);
    result.add_info("NS", std::to_string(genotype_calls.size()));
    result.set_format({"GT", "GP"});
    
    if (variants.empty()) {
        result.set_filters({"REFCALL"});
    }
    
    for (const auto& sample_call : genotype_calls) {
        const auto& sample = sample_call.first;
        result.add_genotype(sample, to_vcf_genotype(sample_call.second.first), VcfRecord::Builder::Phasing::Phased);
        result.add_genotype_field(sample, "GP", std::to_string(to_phred_quality(sample_call.second.second)));
    }
    
    return result.build_once();
}

std::vector<VcfRecord> PopulationVariantCaller::call_variants(const GenomicRegion& region,
                                                              const std::vector<Variant>& candidates,
                                                              const ReadMap& reads)
{
    std::vector<VcfRecord> result {};
    
    HaplotypeTree tree {reference_};
    extend_tree(candidates, tree);
    
    auto haplotypes = tree.get_haplotypes(region);
    
    //std::cout << "there are " << haplotypes.size() << " haplotypes" << std::endl;
    
    GenotypeModel::Population genotype_model {ploidy_};
    
    auto genotype_posteriors = genotype_model.evaluate(haplotypes, reads, reference_).genotype_posteriors;
    
//    for (auto& g : genotype_posteriors.at("HG00101")) {
//        print_alleles(g.first);
//        std::cout << g.second << std::endl;
//    }
    
    auto alleles = decompose(candidates);
    
    auto allele_posteriors = compute_allele_posteriors(genotype_posteriors, haplotypes, alleles);
    
//    for (const auto& ap : allele_posteriors.at("HG00102")) {
//        std::cout << ap.first << " " << ap.second << std::endl;
//    }
    
    auto segments = segment(alleles);
    auto regions  = get_segment_regions(segments);
    
    auto genotype_calls = call_genotypes(genotype_posteriors, regions);
    auto genotype_calls_itr = std::cbegin(genotype_calls);
    
    for (const auto& segment : segments) {
        auto variants = call_segment_variants(segment, allele_posteriors, min_posterior_);
        
        const auto& reference_allele = segment.front();
        
        if (variants.empty() && make_ref_calls_) {
            auto homozygous_reference_posterior = get_homozygous_reference_posterior(reference_allele, *genotype_calls_itr);
            
            if (homozygous_reference_posterior >= min_posterior_) {
                auto segment_quality = to_phred_quality(homozygous_reference_posterior);
                result.emplace_back(call_segment(reference_allele, variants, *genotype_calls_itr, segment_quality, reference_));
            }
        } else {
            auto segment_quality = to_phred_quality(max_posterior(variants, allele_posteriors));
            result.emplace_back(call_segment(reference_allele, variants, *genotype_calls_itr, segment_quality, reference_));
        }
        
        ++genotype_calls_itr;
    }
    
    return result;
}
    
} // namespace Octopus
