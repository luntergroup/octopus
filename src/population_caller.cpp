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
#include <utility>
#include <algorithm>
#include <iterator>

#include <boost/optional.hpp>

#include "common.hpp"
#include "genomic_region.hpp"
#include "read_pipe.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "merge_transform.hpp"
#include "vcf_record.hpp"
#include "maths.hpp"
#include "mappable_algorithms.hpp"
#include "variant_utils.hpp"
#include "read_utils.hpp"
#include "string_utils.hpp"
#include "sequence_utils.hpp"
#include "random_candidate_variant_generator.hpp"
#include "haplotype_prior_model.hpp"
#include "haplotype_phaser.hpp"

#include <iostream> // TEST

namespace Octopus
{

// public methods

PopulationVariantCaller::PopulationVariantCaller(const ReferenceGenome& reference,
                                                 ReadPipe& read_pipe,
                                                 CandidateVariantGenerator&& candidate_generator,
                                                 RefCallType refcall_type, double min_variant_posterior,
                                                 double min_refcall_posterior, unsigned ploidy)
:
VariantCaller {reference, read_pipe, std::move(candidate_generator), refcall_type},
genotype_model_ {ploidy},
ploidy_ {ploidy},
min_variant_posterior_ {min_variant_posterior},
min_refcall_posterior_ {min_refcall_posterior}
{}

std::string PopulationVariantCaller::do_get_details() const
{
    return "population caller. ploidy = " + std::to_string(ploidy_);
}

// non member methods

namespace {

using GM = GenotypeModel::Population;

struct GenotypeCall
{
    GenotypeCall() = default;
    template <typename T>
    GenotypeCall(T&& genotype, double posterior) : genotype {std::forward<T>(genotype)}, posterior {posterior} {}
    
    //const GenomicRegion& get_region() const { return genotype.get_region(); }
    
    Genotype<Allele> genotype;
    double posterior;
    
    GenomicRegion phase_region;
    double phase_score;
};

using GenotypeCalls = std::unordered_map<SampleIdType, GenotypeCall>;

using AllelePosteriors = std::unordered_map<SampleIdType, std::unordered_map<Allele, double>>;

struct VariantCall : public Mappable<VariantCall>
{
    VariantCall() = default;
    template <typename T>
    VariantCall(T&& variant, double posterior) : variant {std::forward<T>(variant)}, posterior {posterior} {}
    
    GenomicRegion get_region() const { return variant.get_region(); }
    
    Variant variant;
    double posterior;
};

using VariantCalls = std::vector<VariantCall>;

struct VariantCallBlock : public Mappable<VariantCallBlock>
{
    VariantCallBlock() = default;
    template <typename T>
    VariantCallBlock(T&& variants, double posterior) : variants {std::forward<T>(variants)}, posterior {posterior} {}
    
    GenomicRegion get_region() const { return variants.front().get_region(); }
    
    std::vector<Variant> variants;
    double posterior;
};

using VariantCallBlocks = std::vector<VariantCallBlock>;

struct RefCall : public Mappable<RefCall>
{
    RefCall() = default;
    template <typename A, typename T>
    RefCall(A&& reference_allele, double posterior, T&& sample_posteriors)
    :
    reference_allele {std::forward<A>(reference_allele)},
    posterior {posterior},
    sample_posteriors {std::forward<T>(sample_posteriors)}
    {}
    
    const GenomicRegion& get_region() const { return reference_allele.get_region(); }
    
    Allele reference_allele;
    double posterior;
    std::vector<std::pair<SampleIdType, double>> sample_posteriors;
};

using RefCalls = std::vector<RefCall>;

} // namespace

namespace debug {
    void print_genotype_posteriors(const GM::GenotypeProbabilities& genotype_posteriors, size_t n = 5);
    void print_allele_posteriors(const AllelePosteriors& allele_posteriors, size_t n = 10);
    void print_variant_calls(const VariantCallBlocks& calls);
    void print_genotype_calls(const GenotypeCalls& calls);
    void print_genotype_calls(const std::vector<GenotypeCalls>& calls);
} // namespace debug

namespace {
    
void remove_low_posterior_genotypes(GM::GenotypeProbabilities& genotype_posteriors,
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

double marginalise(const Allele& allele,
                   const GM::SampleGenotypeProbabilities& genotype_posteriors)
{
    double result {0.0};
    
    for (const auto& genotype_posterior : genotype_posteriors) {
        if (contains(genotype_posterior.first, allele)) {
            result += genotype_posterior.second;
        }
    }
    
    return result;
}

std::unordered_map<Allele, double>
compute_sample_allele_posteriors(const GM::SampleGenotypeProbabilities& genotype_posteriors,
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
compute_allele_posteriors(const GM::GenotypeProbabilities& genotype_posteriors,
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
marginalise(const GenomicRegion& region,
            const GM::SampleGenotypeProbabilities& genotype_posteriors)
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

double marginalise(const Genotype<Allele>& genotype,
                          const GM::SampleGenotypeProbabilities& genotype_posteriors)
{
    return std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), 0.0,
                           [&genotype] (double curr, const auto& p) {
                               return curr + ((contains(p.first, genotype)) ? p.second : 0.0);
                           });
}

std::vector<GenotypeCalls>
call_genotypes(const GM::GenotypeProbabilities& genotype_posteriors,
               const std::vector<GenomicRegion>& regions)
{
    std::vector<GenotypeCalls> result(regions.size());
    
    for (const auto sample_genotype_posteriors : genotype_posteriors) {
        const auto& sample_genotype_call = call_genotype(sample_genotype_posteriors.second);
        
        for (size_t i {}; i < regions.size(); ++i) {
            const auto spliced_genotype = splice<Allele>(sample_genotype_call.first, regions[i]);
            
            const auto posterior = marginalise(spliced_genotype, sample_genotype_posteriors.second);
            
            result[i].emplace(sample_genotype_posteriors.first, GenotypeCall(std::move(spliced_genotype), posterior));
        }
    }
    
    return result;
}

double max_posterior(const Allele& allele, const AllelePosteriors& allele_posteriors)
{
    double result {};
    
    for (const auto& sample_allele_posteriors : allele_posteriors) {
        auto curr = sample_allele_posteriors.second.at(allele);
        if (curr > result) result = curr;
    }
    
    return result;
}

VariantCalls call_variants(const std::vector<Variant>& candidates,
                           const AllelePosteriors& allele_posteriors,
                           const double min_posterior)
{
    VariantCalls result {};
    result.reserve(candidates.size());
    
    for (const Variant& candidate : candidates) {
        auto posterior = max_posterior(candidate.get_alt_allele(), allele_posteriors);
        
        if (posterior >= min_posterior) {
            result.emplace_back(candidate, posterior);
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

double max_posterior(const VariantCalls& calls)
{
    return std::max_element(std::cbegin(calls), std::cend(calls),
                            [] (const auto& lhs, const auto& rhs) {
                                return lhs.posterior < rhs.posterior;
                            })->posterior;
}

VariantCallBlocks block_variant_calls(const std::vector<VariantCalls>& segmented_calls)
{
    VariantCallBlocks result {};
    result.reserve(segmented_calls.size());
    
    for (const auto& calls : segmented_calls) {
        std::vector<Variant> variants {};
        variants.reserve(calls.size());
        std::transform(std::cbegin(calls), std::cend(calls), std::back_inserter(variants),
                       [] (const auto& call) { return call.variant; });
        result.emplace_back(std::move(variants), max_posterior(calls));
    }
    
    return result;
}

VariantCallBlocks call_blocked_variants(const std::vector<Variant>& candidates,
                                        const AllelePosteriors& allele_posteriors,
                                        const double min_posterior)
{
    return block_variant_calls(segment_overlapped(call_variants(candidates, allele_posteriors, min_posterior)));
}

auto get_regions(const VariantCallBlocks& variant_calls)
{
    std::vector<GenomicRegion> result {};
    result.reserve(variant_calls.size());
    
    for (const auto& segment_calls : variant_calls) {
        result.emplace_back(get_encompassing_region(segment_calls.variants));
    }
    
    return result;
}

void parsimonise_variant_calls(VariantCallBlocks& variant_calls, const ReferenceGenome& reference)
{
    for (auto& segment_calls : variant_calls) {
        segment_calls.variants = parsimonise_together(segment_calls.variants, reference);
    }
}

double compute_homozygous_reference_posterior(const Allele& reference_allele,
                                              const GenotypeCalls& genotype_calls)
{
    if (genotype_calls.empty()) return 0.0;
    
    double result {1.0};
    
    for (const auto& sample_genotype_call : genotype_calls) {
        if (is_homozygous(sample_genotype_call.second.genotype, reference_allele)) {
            auto curr = sample_genotype_call.second.posterior;
            if (curr < result) {
                result = curr;
            }
        } else {
            return 0.0;
        }
    }
    
    return result;
}

double marginalise_reference_genotype(const Allele& reference_allele,
                                      const GM::SampleGenotypeProbabilities& sample_genotype_posteriors)
{
    double result {};
    
    for (const auto& genotype_posterior : sample_genotype_posteriors) {
        if (is_homozygous(genotype_posterior.first, reference_allele)) {
            result += genotype_posterior.second;
        }
    }
    
    return result;
}

RefCalls call_reference(const GM::GenotypeProbabilities& genotype_posteriors,
                        const std::vector<Allele>& reference_alleles,
                        const ReadMap& reads, double min_posterior)
{
    RefCalls result {};
    
    if (reference_alleles.empty()) return result;
    
    result.reserve(reference_alleles.size());
    
    for (const auto& reference_allele : reference_alleles) {
        double min_sample_posteior {1.0};
        
        std::vector<std::pair<SampleIdType, double>> sample_posteriors {};
        sample_posteriors.reserve(genotype_posteriors.size());
        
        for (const auto& sample_genotype_posteriors : genotype_posteriors) {
            double sample_posterior {};
            
            if (min_coverage(reads.at(sample_genotype_posteriors.first), get_region(reference_allele)) > 0) {
                sample_posterior = marginalise_reference_genotype(reference_allele, sample_genotype_posteriors.second);
                
                if (sample_posterior < min_posterior) {
                    min_sample_posteior = sample_posterior;
                    break; // to avoid computing the rest
                } else if (sample_posterior < min_sample_posteior) {
                    min_sample_posteior = sample_posterior;
                }
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
        for (const auto& allele : sample_call.second.genotype) {
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
        for (const auto& allele : sample_call.second.genotype) {
            ++allele_counts[allele];
        }
    }
    
    std::vector<unsigned> result {};
    result.reserve(variants.size());
    
    for (const auto& variant : variants) {
        result.push_back(allele_counts[variant.get_alt_allele()]);
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

VcfRecord output_variant_call(const VariantCallBlock& block, const GenotypeCalls& genotype_calls,
                              const ReferenceGenome& reference, const ReadMap& reads)
{
    using std::to_string;
    
    auto result = VcfRecord::Builder();
    
    const auto phred_quality = Maths::probability_to_phred<float>(block.posterior);
    
    const auto reference_allele = block.variants.front().get_ref_allele();
    
    const auto& region = get_region(reference_allele);
    
    result.set_chromosome(get_contig_name(region));
    result.set_position(get_begin(region));
    result.set_ref_allele(reference_allele.get_sequence());
    result.set_alt_alleles(get_alt_allele_sequences(block.variants));
    result.set_quality(phred_quality);
    
    //result.set_filters({"PASS"}); // TODO
    
    result.add_info("AC",  to_strings(count_alleles(block.variants, genotype_calls)));
    result.add_info("AN",  to_string(count_alleles(genotype_calls)));
    result.add_info("NS",  to_string(count_samples_with_coverage(reads, region)));
    result.add_info("DP",  to_string(sum_max_coverages(reads, region)));
    result.add_info("SB",  Octopus::to_string(strand_bias(reads, region), 2));
    result.add_info("BQ",  to_string(static_cast<unsigned>(rmq_base_quality(reads, region))));
    result.add_info("MQ",  to_string(static_cast<unsigned>(rmq_mapping_quality(reads, region))));
    result.add_info("MQ0", to_string(count_mapq_zero(reads, region)));
    
    result.set_format({"GT", "FT", "GQ", "PS", "PQ", "DP", "BQ", "MQ"});
    
    for (const auto& sample_call : genotype_calls) {
        const auto& sample = sample_call.first;
        result.add_genotype(sample, to_vcf_genotype(sample_call.second.genotype), VcfRecord::Builder::Phasing::Phased);
        
        result.add_genotype_field(sample, "FT", "."); // TODO
        result.add_genotype_field(sample, "GQ", Octopus::to_string(Maths::probability_to_phred<float>(sample_call.second.posterior), 2));
        result.add_genotype_field(sample, "PS", to_string(get_begin(sample_call.second.phase_region) + 1));
        result.add_genotype_field(sample, "PQ", Octopus::to_string(Maths::probability_to_phred<float>(sample_call.second.phase_score), 2)); // TODO
        result.add_genotype_field(sample, "DP", to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result.build_once();
}

VcfRecord output_reference_call(RefCall call, const ReferenceGenome& reference,
                                const ReadMap& reads, unsigned ploidy)
{
    using std::to_string;
    
    auto result = VcfRecord::Builder();
    
    const auto phred_quality = Maths::probability_to_phred(call.posterior);
    
    const auto& region = get_region(call.reference_allele);
    
    result.set_chromosome(get_contig_name(region));
    result.set_position(get_begin(region));
    result.set_ref_allele(call.reference_allele.get_sequence().front());
    
    result.set_alt_allele("<NON_REF>");
    
    result.set_quality(phred_quality);
    
    //result.set_filters({"REFCALL"});
    
    if (size(region) > 1) {
        result.add_info("END", to_string(get_end(region))); // - 1 as VCF uses closed intervals
    }
    
    result.add_info("NS",  to_string(count_samples_with_coverage(reads, region)));
    result.add_info("DP",  to_string(sum_max_coverages(reads, region)));
    result.add_info("SB",  Octopus::to_string(strand_bias(reads, region), 2));
    result.add_info("BQ",  to_string(static_cast<unsigned>(rmq_base_quality(reads, region))));
    result.add_info("MQ",  to_string(static_cast<unsigned>(rmq_mapping_quality(reads, region))));
    result.add_info("MQ0", to_string(count_mapq_zero(reads, region)));
    
    result.set_format({"GT", "GQ", "DP", "BQ", "MQ"});
    
    for (const auto& sample_posteior : call.sample_posteriors) {
        const auto& sample = sample_posteior.first;
        result.add_homozygous_ref_genotype(sample, ploidy);
        result.add_genotype_field(sample, "GQ", to_string(Maths::probability_to_phred(sample_posteior.second)));
        result.add_genotype_field(sample, "DP", to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result.build_once();
}
    
} // namespace

std::vector<VcfRecord>
PopulationVariantCaller::call_variants(const GenomicRegion& region,
                                       const std::vector<Variant>& candidates,
                                       const ReadMap& reads) const
{
    std::vector<VcfRecord> result {};
    
    if (empty(region) || (candidates.empty() && refcall_type_ == RefCallType::None)) {
        return result;
    }
    
    HaplotypePhaser phaser {region.get_contig_name(), reference_, candidates, reads, 2000, 3};
    
    while (!phaser.done()) {
        auto haplotypes = phaser.get_haplotypes();
        
        make_unique(haplotypes, haplotype_prior_model_);
        
        phaser.set_haplotypes(haplotypes);
        
        std::cout << "there are " << haplotypes.size() << " unique haplotypes" << std::endl;
        
        const auto haplotype_region = get_region(haplotypes.front());
        
        std::cout << "haplotype region is " << haplotype_region << std::endl;
        
        auto haplotype_region_reads = copy_overlapped(reads, haplotype_region);
        
        std::cout << "there are " << count_reads(haplotype_region_reads) << " reads in haplotype region" << std::endl;
        std::cout << "computing latents" << std::endl;
        
        auto latents = genotype_model_.evaluate(haplotypes, haplotype_region_reads, reference_);
        
        auto& genotype_posteriors = latents.genotype_posteriors;
        
        std::cout << "finished computing latents" << std::endl;
        
        const auto phase_set = phaser.phase(haplotypes, genotype_posteriors);
        
        if (phase_set) {
            std::cout << "phased region is " << phase_set->region << std::endl;
            
            auto overlapped_candidates = copy_overlapped(candidates, phase_set->region);
            
            debug::print_genotype_posteriors(genotype_posteriors);
            
            //remove_low_posterior_genotypes(genotype_posteriors, 0.0000000001);
            
            auto alleles = generate_callable_alleles(phase_set->region, overlapped_candidates,
                                                     refcall_type_, reference_);
            
            auto allele_posteriors = compute_allele_posteriors(genotype_posteriors, alleles);
            
            debug::print_allele_posteriors(allele_posteriors);
            
            auto variant_calls = call_blocked_variants(overlapped_candidates, allele_posteriors,
                                                       min_variant_posterior_);
            
            //debug::print_variant_calls(variant_calls);
            
            parsimonise_variant_calls(variant_calls, reference_);
            
            auto called_regions = get_regions(variant_calls);
            
            auto variant_genotype_calls = call_genotypes(genotype_posteriors, called_regions);
            
            // TODO
            for (auto& g : variant_genotype_calls) {
                for (auto& p : g) {
                    const auto& phase_regions = phase_set->phase_regions.at(p.first);
                    const auto& call_region = p.second.genotype[0].get_region();
                    auto it = std::find_if(std::cbegin(phase_regions), std::cend(phase_regions),
                                           [&call_region] (const auto& phase_region) {
                                               return contains(phase_region.region, call_region);
                                           });
                    
                    auto it2 = std::find_if(std::cbegin(called_regions), std::cend(called_regions),
                                            [it] (const auto& region) {
                                                return contains(it->region, region);
                                            });
                    
                    p.second.phase_region = *it2;
                    p.second.phase_score  = it->score;
                }
            }
            
            for (const auto& p : phase_set->phase_regions) {
                std::cout << "phase regions for sample " << p.first << std::endl;
                for (const auto& r : p.second) {
                    std::cout << r.region << " " << r.score << std::endl;
                }
            }
            
            debug::print_genotype_calls(variant_genotype_calls);
            
            auto candidate_ref_alleles = generate_candidate_reference_alleles(alleles, called_regions,
                                                                              overlapped_candidates,
                                                                              refcall_type_);
            
            //    std::cout << "candidate refcall alleles are" << std::endl;
            //    for (const auto& allele : candidate_ref_alleles) {
            //        std::cout << allele << std::endl;
            //    }
            
            auto refcalls = call_reference(genotype_posteriors, candidate_ref_alleles, reads,
                                           min_refcall_posterior_);
            
            //std::cout << called_regions.front() << std::endl;
            
            result.reserve(variant_calls.size() + refcalls.size());
            
            merge_transform(variant_calls, variant_genotype_calls, refcalls, std::back_inserter(result),
                            [this, &haplotype_region_reads] (const auto& variant_call,
                                                             const auto& genotype_call) {
                                return output_variant_call(variant_call, genotype_call, reference_,
                                                           haplotype_region_reads);
                            },
                            [this, &haplotype_region_reads] (const auto& refcall) {
                                return output_reference_call(refcall, reference_,
                                                             haplotype_region_reads, ploidy_);
                            });
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}
    
    namespace debug {
        
        void print_genotype_posteriors(const GM::GenotypeProbabilities& genotype_posteriors, const size_t n)
        {
            for (const auto& sample_posteriors : genotype_posteriors) {
                auto m = std::min(n, sample_posteriors.second.size());
                std::cout << "printing top " << m << " genotype posteriors for sample " << sample_posteriors.first << std::endl;
                
                std::vector<std::pair<Genotype<Haplotype>, double>> v {};
                v.reserve(sample_posteriors.second.size());
                
                std::copy(std::cbegin(sample_posteriors.second), std::cend(sample_posteriors.second),
                          std::back_inserter(v));
                
                std::sort(std::begin(v), std::end(v), [] (const auto& lhs, const auto& rhs) {
                    return lhs.second > rhs.second;
                });
                
                for (size_t i {}; i < m; ++i) {
                    std::cout << "\t* ";
                    print_variant_alleles(v[i].first);
                    std::cout << " " << std::setprecision(20) << v[i].second << std::endl;
                }
            }
        }
        
        void print_allele_posteriors(const AllelePosteriors& allele_posteriors, const size_t n)
        {
            for (const auto& sample_posteriors : allele_posteriors) {
                auto m = std::min(n, sample_posteriors.second.size());
                std::cout << "printing top " << m << " allele posteriors for sample " << sample_posteriors.first << std::endl;
                
                std::vector<std::pair<Allele, double>> v {};
                v.reserve(sample_posteriors.second.size());
                
                std::copy(std::cbegin(sample_posteriors.second), std::cend(sample_posteriors.second),
                          std::back_inserter(v));
                
                std::sort(std::begin(v), std::end(v), [] (const auto& lhs, const auto& rhs) {
                    return lhs.second > rhs.second;
                });
                
                for (size_t i {}; i < m; ++i) {
                    std::cout << "\t* ";
                    std::cout << v[i].first << " " << std::setprecision(20) << v[i].second << std::endl;
                }
            }
        }
        
        void print_variant_calls(const VariantCallBlocks& calls)
        {
            std::cout << "printing all variant calls" << std::endl;
            for (const auto& segment_calls : calls) {
                std::cout << "printing calls in segment " << get_encompassing_region(segment_calls.variants) << std::endl;
                std::copy(std::cbegin(segment_calls.variants), std::cend(segment_calls.variants),
                          std::ostream_iterator<Variant>(std::cout, "\n"));
            }
        }
        
        void print_genotype_calls(const GenotypeCalls& calls)
        {
            std::cout << "printing genotype calls for each sample" << std::endl;
            for (const auto& call : calls) {
                std::cout << call.first << " " << call.second.genotype << " " << call.second.posterior << std::endl;
            }
        }
        
        void print_genotype_calls(const std::vector<GenotypeCalls>& calls)
        {
            for (const auto& c : calls) {
                print_genotype_calls(c);
            }
        }
    } // namespace debug
    
} // namespace Octopus
