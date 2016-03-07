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
#include <deque>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <utility>

#include "genomic_region.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "vcf_record.hpp"
#include "maths.hpp"
#include "mappable_algorithms.hpp"
#include "read_utils.hpp"
#include "string_utils.hpp"
#include "sequence_utils.hpp"
#include "haplotype_prior_model.hpp"
#include "probability_matrix.hpp"
#include "merge_transform.hpp"

#include <iostream> // DEBUG

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

PopulationVariantCaller::Latents::Latents(GenotypeModel::Population::Latents&& model_latents)
:
haplotype_frequencies_ {std::make_shared<HaplotypePosteiorMap>(std::move(model_latents.haplotype_frequencies))},
genotype_posteriors_ {std::make_shared<GenotypePosteriorMap>(std::move(model_latents.genotype_posteriors))}
{}
    
PopulationVariantCaller::Latents::Latents(HaplotypePosteiorMap&& haplotype_posteriors,
                                          GenotypePosteriorMap&& genotype_posteriors)
:
haplotype_frequencies_ {std::make_shared<HaplotypePosteiorMap>(std::move(haplotype_posteriors))},
genotype_posteriors_ {std::make_shared<GenotypePosteriorMap>(std::move(genotype_posteriors))}
{}

std::unique_ptr<PopulationVariantCaller::CallerLatents>
PopulationVariantCaller::infer_latents(const std::vector<SampleIdType>& samples,
                                       const std::vector<Haplotype>& haplotypes,
                                       const HaplotypePriorMap& haplotype_priors,
                                       const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    auto model_latents = genotype_model_.infer_latents(samples, haplotypes, haplotype_priors,
                                                       haplotype_likelihoods);
    
    auto haplotype_posteriors = marginalise_haplotypes(haplotypes, model_latents.genotype_posteriors);
    
    return std::make_unique<Latents>(std::move(haplotype_posteriors),
                                     std::move(model_latents.genotype_posteriors));
}

// non member methods

namespace
{
using GM = GenotypeModel::Population;
using GenotypePosteriorMap       = GM::Latents::GenotypeProbabilityMap;
using SampleGenotypePosteriorMap = GenotypePosteriorMap::InnerMap;

struct GenotypeCall
{
    GenotypeCall() = default;
    
    template <typename T> GenotypeCall(T&& genotype, double posterior)
    :
    genotype {std::forward<T>(genotype)},
    posterior {posterior}
    {}
    
    GenomicRegion get_region() const { return phase_region; }
    
    Genotype<Allele> genotype;
    double posterior;
    GenomicRegion phase_region;
    double phase_score;
};

using GenotypeCallMap = std::unordered_map<SampleIdType, GenotypeCall>;
using GenotypeCalls   = std::vector<GenotypeCallMap>;

using AllelePosteriorMap       = ProbabilityMatrix<Allele>;
using SampleAllelePosteriorMap = AllelePosteriorMap::InnerMap;

struct VariantCall : public Mappable<VariantCall>
{
    VariantCall() = default;
    template <typename T>
    VariantCall(T&& variant, double posterior)
    :
    variant {std::forward<T>(variant)},
    posterior {posterior}
    {}
    
    GenomicRegion get_region() const { return variant.get_region(); }
    
    Variant variant;
    double posterior;
};

using VariantCalls = std::vector<VariantCall>;

struct VariantCallBlock : public Mappable<VariantCallBlock>
{
    VariantCallBlock() = default;
    template <typename T>
    VariantCallBlock(T&& variants, double posterior)
    :
    variants {std::forward<T>(variants)},
    posterior {posterior}
    {}
    
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

namespace debug
{
    void print_genotype_posteriors(const GenotypePosteriorMap& genotype_posteriors, std::size_t n = 5);
    void print_allele_posteriors(const AllelePosteriorMap& allele_posteriors, std::size_t n = 10);
    void print_variant_calls(const VariantCallBlocks& calls);
    void print_genotype_calls(const GenotypeCallMap& calls);
    void print_genotype_calls(const std::vector<GenotypeCallMap>& calls);
} // namespace debug

namespace
{

// allele posterior calculations

using AlleleBools          = std::deque<bool>; // using std::deque because std::vector<bool> is evil
using GenotypeContainments = std::vector<AlleleBools>;

auto marginalise(const SampleGenotypePosteriorMap& genotype_posteriors,
                 const AlleleBools& contained_alleles)
{
    auto result = std::inner_product(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                     std::cbegin(contained_alleles), 0.0, std::plus<void> {},
                                     [] (const auto& p, const auto is_contained) {
                                         return (is_contained) ? p.second : 0.0;
                                     });
    
    if (result > 1.0) result = 1.0; // just to account for floating point error
    
    return result;
}

auto compute_sample_allele_posteriors(const SampleGenotypePosteriorMap& genotype_posteriors,
                                      const GenotypeContainments& contained_alleles)
{
    std::vector<double> result {};
    result.reserve(contained_alleles.size());
    
    for (const auto& allele : contained_alleles) {
        result.emplace_back(marginalise(genotype_posteriors, allele));
    }
    
    return result;
}

auto find_contained_alleles(const GenotypePosteriorMap& genotype_posteriors,
                            const std::vector<Allele>& alleles)
{
    const auto num_genotypes = genotype_posteriors.size2();
    
    GenotypeContainments result {};
    
    if (num_genotypes == 0 || genotype_posteriors.empty1() || alleles.empty()) {
        return result;
    }
    
    result.reserve(alleles.size());
    
    const auto& test_sample   = genotype_posteriors.begin()->first;
    const auto genotype_begin = genotype_posteriors.begin(test_sample);
    const auto genotype_end   = genotype_posteriors.end(test_sample);
    
    for (const auto& allele : alleles) {
        result.emplace_back(num_genotypes);
        
        std::transform(genotype_begin, genotype_end, std::begin(result.back()),
                       [&] (const auto& p) {
                           return contains(p.first, allele);
                       });
    }
    
    return result;
}

AllelePosteriorMap compute_allele_posteriors(const GenotypePosteriorMap& genotype_posteriors,
                                             const std::vector<Allele>& alleles)
{
    AllelePosteriorMap result {std::cbegin(alleles), std::cend(alleles)};
    result.reserve1(genotype_posteriors.size1());
    
    const auto contained_alleles = find_contained_alleles(genotype_posteriors, alleles);
    
    for (const auto& sample_genotype_posteriors : genotype_posteriors) {
        insert_sample(sample_genotype_posteriors.first,
                      compute_sample_allele_posteriors(sample_genotype_posteriors.second, contained_alleles),
                      result);
    }
    
    return result;
}

// variant calling

double max_posterior(const Allele& allele, const AllelePosteriorMap& allele_posteriors)
{
    double result {0};
    
    for (const auto& sample_allele_posteriors : allele_posteriors) {
        const auto curr = sample_allele_posteriors.second.at(allele);
        
        if (curr > result) result = curr;
    }
    
    return result;
}

VariantCalls call_variants(const std::vector<Variant>& candidates,
                           const AllelePosteriorMap& allele_posteriors,
                           const double min_posterior)
{
    VariantCalls result {};
    result.reserve(candidates.size());
    
    for (const Variant& candidate : candidates) {
        const auto posterior = max_posterior(candidate.get_alt_allele(), allele_posteriors);
        
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
                            [] (const VariantCall& lhs, const VariantCall& rhs) {
                                return lhs.posterior < rhs.posterior;
                            })->posterior;
}

VariantCallBlocks block_variant_calls(std::vector<VariantCalls>&& segmented_calls)
{
    using std::begin; using std::end; using std::make_move_iterator; using std::for_each;
    
    VariantCallBlocks result {};
    result.reserve(segmented_calls.size());
    
    for_each(make_move_iterator(begin(segmented_calls)), make_move_iterator(end(segmented_calls)),
             [&result] (VariantCalls&& calls) {
                 std::vector<Variant> variants {};
                 variants.reserve(calls.size());
                 
                 std::transform(make_move_iterator(begin(calls)),
                                make_move_iterator(end(calls)),
                                std::back_inserter(variants),
                                [] (VariantCall&& call) {
                                    return std::move(call.variant);
                                });
                 
                 result.emplace_back(std::move(variants), max_posterior(calls));
             });
    
    return result;
}

auto segment_calls(VariantCalls&& calls)
{
    return segment_overlapped(std::make_move_iterator(std::begin(calls)),
                              std::make_move_iterator(std::end(calls)));
}

auto call_blocked_variants(const std::vector<Variant>& candidates,
                           const AllelePosteriorMap& allele_posteriors,
                           const double min_posterior)
{
    auto calls = call_variants(candidates, allele_posteriors, min_posterior);
    return block_variant_calls(segment_calls(std::move(calls)));
}

auto extract_regions(const VariantCallBlocks& variant_calls)
{
    std::vector<GenomicRegion> result {};
    result.reserve(variant_calls.size());
    
    for (const auto& segment_calls : variant_calls) {
        result.emplace_back(encompassing_region(segment_calls.variants));
    }
    
    return result;
}

void parsimonise_variant_calls(VariantCallBlocks& variant_calls, const ReferenceGenome& reference)
{
    for (auto& segment_calls : variant_calls) {
        segment_calls.variants = parsimonise_together(segment_calls.variants, reference);
    }
}

// variant genotype calling

auto call_genotype(const SampleGenotypePosteriorMap& genotype_posteriors)
{
    return *std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                             [] (const auto& lhs, const auto& rhs) {
                                 return lhs.second < rhs.second;
                             });
}

double marginalise(const Genotype<Allele>& genotype, const SampleGenotypePosteriorMap& genotype_posteriors)
{
    return std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), 0.0,
                           [&genotype] (const double curr, const auto& p) {
                               return curr + ((contains(p.first, genotype)) ? p.second : 0.0);
                           });
}

// Call cleanup

GenotypeCalls call_genotypes(const GenotypePosteriorMap& genotype_posteriors,
                             const std::vector<GenomicRegion>& variant_regions)
{
    GenotypeCalls result(variant_regions.size());
    
    for (const auto& sample_genotype_posteriors : genotype_posteriors) {
        const auto& sample_genotype_call = call_genotype(sample_genotype_posteriors.second);
        
        for (std::size_t i {0}; i < variant_regions.size(); ++i) {
            auto spliced_genotype = splice<Allele>(sample_genotype_call.first, variant_regions[i]);
            
            const auto posterior = marginalise(spliced_genotype, sample_genotype_posteriors.second);
            
            result[i].emplace(std::piecewise_construct,
                              std::forward_as_tuple(sample_genotype_posteriors.first),
                              std::forward_as_tuple(std::move(spliced_genotype), posterior));
        }
    }
    
    return result;
}

bool is_genotyped(const Variant& variant, const GenotypeCallMap& genotype_calls)
{
    return std::any_of(std::cbegin(genotype_calls), std::cend(genotype_calls),
                       [&variant] (const auto& p) {
                           return p.second.genotype.contains(variant.get_alt_allele());
                       });
}

std::size_t remove_non_genotyped_calls(VariantCallBlock& variant_calls,
                                    const GenotypeCallMap& genotype_calls)
{
    assert(!variant_calls.variants.empty());
    
    if (variant_calls.variants.size() == 1) {
        return 0;
    }
    
    auto& variants = variant_calls.variants;
    
    // TODO: to stop missing some edge cases, shouldn't need this when we figure a better way
    // to represent overlapping variants in VCF
    if (std::all_of(std::begin(variants), std::end(variants),
                    [&genotype_calls] (const Variant& call) {
                        return !is_genotyped(call, genotype_calls);
                    })) return 0;
    
    const auto it = std::remove_if(std::begin(variants), std::end(variants),
                                   [&genotype_calls] (const Variant& call) {
                                       return !is_genotyped(call, genotype_calls);
                                   });
    
    const auto result = std::distance(it, std::end(variants));
    
    variants.erase(it, std::end(variants));
    
    return result;
}

std::size_t remove_non_genotyped_calls(VariantCallBlocks& variant_calls,
                                    const GenotypeCalls& genotype_calls)
{
    std::size_t result {0};
    
    auto it = std::cbegin(genotype_calls);
    
    for (auto& block : variant_calls) {
        result += remove_non_genotyped_calls(block, *it++);
    }
    
    return result;
}

// reference genotype calling

double marginalise_reference_genotype(const Allele& reference_allele,
                                      const SampleGenotypePosteriorMap& sample_genotype_posteriors)
{
    double result {0};
    
    for (const auto& genotype_posterior : sample_genotype_posteriors) {
        if (is_homozygous(genotype_posterior.first, reference_allele)) {
            result += genotype_posterior.second;
        }
    }
    
    return result;
}

RefCalls call_reference(const GenotypePosteriorMap& genotype_posteriors,
                        const std::vector<Allele>& reference_alleles,
                        const ReadMap& reads, double min_posterior)
{
    RefCalls result {};
    
    if (reference_alleles.empty()) return result;
    
    result.reserve(reference_alleles.size());
    
    for (const auto& reference_allele : reference_alleles) {
        double min_sample_posteior {1.0};
        
        std::vector<std::pair<SampleIdType, double>> sample_posteriors {};
        sample_posteriors.reserve(genotype_posteriors.size1());
        
        for (const auto& sample_genotype_posteriors : genotype_posteriors) {
            double sample_posterior {};
            
            if (min_coverage(reads.at(sample_genotype_posteriors.first), mapped_region(reference_allele)) > 0) {
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

unsigned count_alleles(const GenotypeCallMap& genotype_calls)
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

auto count_alleles(const std::vector<Variant>& variants, const GenotypeCallMap& genotype_calls)
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

void set_vcf_genotype(const SampleIdType& sample, const GenotypeCall& genotype_call,
                      const VariantCallBlock& variant_call, VcfRecord::Builder& record)
{
    std::vector<VcfRecord::SequenceType> result {};
    result.reserve(genotype_call.genotype.ploidy());
    
    for (const auto& allele : genotype_call.genotype) {
        result.push_back(allele.get_sequence());
    }
    
    record.add_genotype(sample, result, VcfRecord::Builder::Phasing::Phased);
}

VcfRecord::Builder output_variant_call(VariantCallBlock&& block,
                                       GenotypeCallMap&& genotype_calls,
                                       const ReferenceGenome& reference,
                                       const ReadMap& reads)
{
    assert(!block.variants.empty());
    
    using std::to_string;
    
    auto result = VcfRecord::Builder {};
    
    const auto phred_quality = Maths::probability_to_phred<float>(block.posterior);
    
    const auto& reference_allele = block.variants.front().get_ref_allele();
    
    const auto& region = mapped_region(reference_allele);
    
    result.set_chromosome(contig_name(region));
    result.set_position(region_begin(region));
    result.set_ref_allele(reference_allele.get_sequence());
    result.set_alt_alleles(extract_alt_allele_sequences(block.variants));
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
        
        set_vcf_genotype(sample, sample_call.second, block, result);
        
        result.add_genotype_field(sample, "FT", "."); // TODO
        result.add_genotype_field(sample, "GQ", Octopus::to_string(Maths::probability_to_phred<float>(sample_call.second.posterior), 2));
        result.add_genotype_field(sample, "PS", to_string(region_begin(sample_call.second.phase_region) + 1));
        result.add_genotype_field(sample, "PQ", Octopus::to_string(Maths::probability_to_phred<float>(sample_call.second.phase_score), 2)); // TODO
        result.add_genotype_field(sample, "DP", to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result;
}

VcfRecord::Builder output_reference_call(RefCall call, const ReferenceGenome& reference,
                                         const ReadMap& reads, unsigned ploidy)
{
    using std::to_string;
    
    auto result = VcfRecord::Builder {};
    
    const auto phred_quality = Maths::probability_to_phred(call.posterior);
    
    const auto& region = mapped_region(call.reference_allele);
    
    result.set_chromosome(contig_name(region));
    result.set_position(region_begin(region));
    result.set_ref_allele(call.reference_allele.get_sequence().front());
    
    result.set_alt_allele("<NON_REF>");
    
    result.set_quality(phred_quality);
    
    if (region_size(region) > 1) {
        result.add_info("END", to_string(region_end(region)));
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
    
    return result;
}

void set_phasings(GenotypeCalls& variant_genotype_calls,
                  const Phaser::PhaseSet& phase_set,
                  const std::vector<GenomicRegion>& called_regions)
{
    // TODO
    for (auto& g : variant_genotype_calls) {
        for (auto& p : g) {
            const auto& phase_regions = phase_set.phase_regions.at(p.first);
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
}
} // namespace

std::vector<VcfRecord::Builder>
merge_calls(VariantCallBlocks&& variant_calls, GenotypeCalls&& variant_genotype_calls,
            RefCalls&& refcalls, const ReferenceGenome& reference, const ReadMap& reads,
            const unsigned ploidy)
{
    using std::begin; using std::end; using std::make_move_iterator; using std::move;
    
    std::vector<VcfRecord::Builder> result {};
    result.reserve(variant_calls.size() + refcalls.size());
    
    merge_transform(make_move_iterator(begin(variant_calls)), make_move_iterator(end(variant_calls)),
                    make_move_iterator(begin(variant_genotype_calls)),
                    make_move_iterator(begin(refcalls)), make_move_iterator(end(refcalls)),
                    std::back_inserter(result),
                    [&] (VariantCallBlock&& variant_call, GenotypeCallMap&& genotype_call) {
                        return output_variant_call(move(variant_call), move(genotype_call),
                                                   reference, reads);
                    },
                    [&] (RefCall&& refcall) {
                        return output_reference_call(move(refcall), reference, reads, ploidy);
                    });
    
    return result;
}

std::vector<VcfRecord::Builder>
PopulationVariantCaller::call_variants(const std::vector<Variant>& candidates,
                                       const std::vector<Allele>& callable_alleles,
                                       CallerLatents* latents,
                                       const Phaser::PhaseSet& phase_set,
                                       const ReadMap& reads) const
{
    const auto dlatents = dynamic_cast<Latents*>(latents);
    
    const auto& genotype_posteriors = *dlatents->genotype_posteriors_;
    
    //debug::print_genotype_posteriors(genotype_posteriors);
    
    const auto allele_posteriors = compute_allele_posteriors(genotype_posteriors, callable_alleles);
    
    //debug::print_allele_posteriors(allele_posteriors);
    
    auto variant_calls = call_blocked_variants(candidates, allele_posteriors, min_variant_posterior_);
    
    //debug::print_variant_calls(variant_calls);
    
    parsimonise_variant_calls(variant_calls, reference_);
    
    //debug::print_variant_calls(variant_calls);
    
    const auto called_regions = extract_regions(variant_calls);
    
    auto variant_genotype_calls = call_genotypes(genotype_posteriors, called_regions);
    
    remove_non_genotyped_calls(variant_calls, variant_genotype_calls);
    
    set_phasings(variant_genotype_calls, phase_set, called_regions); // TODO
    
    //debug::print_genotype_calls(variant_genotype_calls);
    
    auto candidate_ref_alleles = generate_candidate_reference_alleles(callable_alleles, called_regions,
                                                                      candidates, refcall_type_);
    
    auto refcalls = call_reference(genotype_posteriors, candidate_ref_alleles, reads,
                                   min_refcall_posterior_);
    
    return merge_calls(std::move(variant_calls), std::move(variant_genotype_calls),
                       std::move(refcalls), reference_, reads, ploidy_);
}

namespace debug
{
void print_genotype_posteriors(const GenotypePosteriorMap& genotype_posteriors, const std::size_t n)
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
        
        for (std::size_t i {}; i < m; ++i) {
            std::cout << "\t* ";
            print_variant_alleles(v[i].first);
            std::cout << " " << std::setprecision(20) << v[i].second << std::endl;
        }
    }
}

void print_allele_posteriors(const AllelePosteriorMap& allele_posteriors, const std::size_t n)
{
    for (const auto& sample_posteriors : allele_posteriors) {
        const auto m = std::min(n, sample_posteriors.second.size());
        std::cout << "printing top " << m << " allele posteriors for sample " << sample_posteriors.first << std::endl;
        
        std::vector<std::pair<std::reference_wrapper<const Allele>, double>> v {};
        v.reserve(sample_posteriors.second.size());
        
        std::copy(std::cbegin(sample_posteriors.second), std::cend(sample_posteriors.second),
                  std::back_inserter(v));
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [] (const auto& p) {
                          std::cout << "\t* ";
                          std::cout << p.first.get() << " " << std::setprecision(20) << p.second << '\n';
                      });
    }
}

void print_variant_calls(const VariantCallBlocks& calls)
{
    std::cout << "printing all variant calls" << std::endl;
    for (const auto& segment_calls : calls) {
        std::cout << "printing calls in segment " << encompassing_region(segment_calls.variants) << std::endl;
        std::copy(std::cbegin(segment_calls.variants), std::cend(segment_calls.variants),
                  std::ostream_iterator<Variant>(std::cout, "\n"));
    }
}

void print_genotype_calls(const GenotypeCallMap& calls)
{
    std::cout << "printing genotype calls for each sample" << std::endl;
    for (const auto& call : calls) {
        std::cout << call.first << " " << call.second.genotype << " " << call.second.posterior << std::endl;
    }
}

void print_genotype_calls(const std::vector<GenotypeCallMap>& calls)
{
    for (const auto& c : calls) {
        print_genotype_calls(c);
    }
}
} // namespace debug
    
} // namespace Octopus
