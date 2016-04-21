//
//  individual_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "individual_caller.hpp"

#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <utility>
#include <iostream>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

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
#include "probability_matrix.hpp"
#include "merge_transform.hpp"
#include "logging.hpp"

#include "coalescent_model.hpp"
#include "individual_genotype_model.hpp"

namespace Octopus
{
IndividualVariantCaller::CallerParameters::CallerParameters(double min_variant_posterior,
                                                            double min_refcall_posterior,
                                                            unsigned ploidy)
:
min_variant_posterior {min_variant_posterior},
min_refcall_posterior {min_refcall_posterior},
ploidy {ploidy}
{}

IndividualVariantCaller::IndividualVariantCaller(const ReferenceGenome& reference,
                                                 ReadPipe& read_pipe,
                                                 CandidateVariantGenerator&& candidate_generator,
                                                 VariantCaller::CallerParameters general_parameters,
                                                 CallerParameters specific_parameters)
:
VariantCaller {reference, read_pipe, std::move(candidate_generator), std::move(general_parameters)},
sample_ {read_pipe.get_samples().front()},
ploidy_ {specific_parameters.ploidy},
min_variant_log_posterior_ {std::log(specific_parameters.min_variant_posterior)},
min_refcall_log_posterior_ {std::log(specific_parameters.min_refcall_posterior)}
{}

// IndividualVariantCaller::Latents public methods

IndividualVariantCaller::Latents::Latents(const SampleIdType& sample,
                                          const std::vector<Haplotype>& haplotypes,
                                          std::vector<Genotype<Haplotype>>&& genotypes,
                                          ModelLatents&& latents)
:
genotypes_ {std::move(genotypes)},
genotype_log_posteriors_ {std::move(latents.genotype_log_probabilities)}
{
    GenotypeProbabilityMap genotype_posteriors {std::begin(genotypes_), std::end(genotypes_)};
    insert_sample(sample, Maths::exp_copy(genotype_log_posteriors_), genotype_posteriors);
    
    genotype_posteriors_  = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
    haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>(calculate_haplotype_posteriors(haplotypes));
}

std::shared_ptr<IndividualVariantCaller::Latents::HaplotypeProbabilityMap>
IndividualVariantCaller::Latents::get_haplotype_posteriors() const noexcept
{
    return haplotype_posteriors_;
}

std::shared_ptr<IndividualVariantCaller::Latents::GenotypeProbabilityMap>
IndividualVariantCaller::Latents::get_genotype_posteriors() const noexcept
{
    return genotype_posteriors_;
}

// IndividualVariantCaller::Latents private methods

IndividualVariantCaller::Latents::HaplotypeProbabilityMap
IndividualVariantCaller::Latents::calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes)
{
    assert(genotype_posteriors_ != nullptr);
    
    HaplotypeProbabilityMap result {haplotypes.size()};
    
    for (const auto& haplotype : haplotypes) {
        result.emplace(haplotype, 0.0);
    }
    
    const auto& sample = std::cbegin(*genotype_posteriors_)->first;
    
    for (const auto& p : (*genotype_posteriors_)[sample]) {
        for (const auto& haplotype : p.first.copy_unique_ref()) {
            result.at(haplotype) += p.second;
        }
    }
    
    return result;
}

std::unique_ptr<IndividualVariantCaller::CallerLatents>
IndividualVariantCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                                       const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    CoalescentModel prior_model {Haplotype {mapped_region(haplotypes.front()), reference_}};
    
    GenotypeModel::Individual model {ploidy_, prior_model};
    
    auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "There are " << genotypes.size() << " candidate genotypes";
    }
    
    const auto& sample = samples_.front();
    
    auto inferences = model.infer_latents(sample, genotypes, haplotype_likelihoods);
    
//    // TEST
//    GenotypeModel::Individual dummy_model {ploidy_ + 1, prior_model};
//    auto dummy_genotypes = generate_all_genotypes(haplotypes, ploidy_ + 1);
//    auto dummy_inferences = dummy_model.infer_latents(sample, dummy_genotypes, haplotype_likelihoods);
//    auto lp1 = std::log(0.99) + inferences.log_evidence;
//    auto lp2 = std::log(1.0 - 0.99) + dummy_inferences.log_evidence;
//    auto norm = Maths::log_sum_exp(lp1, lp2);
//    lp1 -= norm; lp2 -= norm;
//    if (lp1 < lp2) {
//        std::cout << inferences.log_evidence << " " << dummy_inferences.log_evidence << '\n';
//        auto it = std::max_element(std::cbegin(dummy_inferences.posteriors.genotype_probabilities),
//                                   std::cend(dummy_inferences.posteriors.genotype_probabilities));
//        auto itd = std::distance(std::cbegin(dummy_inferences.posteriors.genotype_probabilities), it);
//        ::debug::print_variant_alleles(dummy_genotypes[itd]);
//        std::cout << " " << *it << '\n';
//        std::cout << haplotypes.front().get_region() << std::endl;
//    }
//    // END TEST
    
    return std::make_unique<Latents>(sample, haplotypes, std::move(genotypes),
                                     std::move(inferences.posteriors));
}

namespace
{
using GM = GenotypeModel::Individual;

using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<Haplotype>>::InnerMap;

using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
using GenotypeLogProbabilityVector = std::vector<std::pair<GenotypeReference, double>>;

GenotypeLogProbabilityVector
zip(const std::vector<Genotype<Haplotype>>& genotypes, const std::vector<double>& log_probabilities)
{
    GenotypeLogProbabilityVector result {};
    result.reserve(genotypes.size());
    
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(log_probabilities),
                   std::back_inserter(result),
                   [] (const auto& g, const auto& p) {
                       return std::make_pair(std::cref(g), p);
                   });
    
    return result;
}

struct GenotypeCall
{
    GenotypeCall() = default;
    
    template <typename T> GenotypeCall(T&& genotype, double posterior)
    :
    genotype {std::forward<T>(genotype)},
    log_posterior {log_posterior}
    {}
    
    GenomicRegion get_region() const { return phase_region; }
    
    Genotype<Allele> genotype;
    double log_posterior;
    GenomicRegion phase_region;
    double phase_score;
};

using GenotypeCalls = std::vector<GenotypeCall>;

using AlleleLogPosteriorMap = std::unordered_map<Allele, double>;

struct VariantCall : public Mappable<VariantCall>
{
    VariantCall() = default;
    template <typename T>
    VariantCall(T&& variant, double log_posterior)
    :
    variant {std::forward<T>(variant)},
    log_posterior {log_posterior}
    {}
    
    GenomicRegion get_region() const { return variant.get_region(); }
    
    Variant variant;
    double log_posterior;
};

using VariantCalls = std::vector<VariantCall>;

struct VariantCallBlock : public Mappable<VariantCallBlock>
{
    VariantCallBlock() = default;
    template <typename T>
    VariantCallBlock(T&& variants, double posterior)
    :
    variants {std::forward<T>(variants)},
    log_posterior {log_posterior}
    {}
    
    const GenomicRegion& get_region() const noexcept { return variants.front().get_region(); }
    
    std::vector<Variant> variants;
    double log_posterior;
};

using VariantCallBlocks = std::vector<VariantCallBlock>;

struct RefCall : public Mappable<RefCall>
{
    RefCall() = default;
    
    template <typename A>
    RefCall(A&& reference_allele, double posterior)
    :
    reference_allele {std::forward<A>(reference_allele)},
    log_posterior {log_posterior}
    {}
    
    const GenomicRegion& get_region() const noexcept { return reference_allele.get_region(); }
    
    Allele reference_allele;
    double log_posterior;
};

using RefCalls = std::vector<RefCall>;

} // namespace

namespace debug
{
    template <typename S>
    void print_genotype_posteriors(S&& stream, const GenotypeLogProbabilityVector& genotype_log_posteriors,
                                   std::size_t n = 5);
    void print_genotype_posteriors(const GenotypeLogProbabilityVector& genotype_log_posteriors,
                                   std::size_t n = 5);
    template <typename S>
    void print_allele_posteriors(S&& stream, const AlleleLogPosteriorMap& allele_posteriors,
                                 std::size_t n = 10);
    void print_allele_posteriors(const AlleleLogPosteriorMap& allele_posteriors,
                                 std::size_t n = 10);
//        void print_variant_calls(const VariantCallBlocks& calls);
} // namespace debug

namespace
{
// allele posterior calculations

auto marginalise(const Allele& allele, const GenotypeLogProbabilityVector& genotype_log_posteriors)
{
    std::vector<double> logs {};
    logs.reserve(genotype_log_posteriors.size());
    
    for (const auto& p : genotype_log_posteriors) {
        if (contains(p.first, allele)) {
            logs.push_back(p.second);
        }
    }
    
    return Maths::log_sum_exp(logs);
}

AlleleLogPosteriorMap compute_allele_posteriors(const GenotypeLogProbabilityVector& genotype_log_posteriors,
                                             const std::vector<Allele>& alleles)
{
    AlleleLogPosteriorMap result {};
    result.reserve(alleles.size());
    
    for (const auto& allele : alleles) {
        result.emplace(allele, marginalise(allele, genotype_log_posteriors));
    }
    
    return result;
}

// variant calling

VariantCalls call_variants(const std::vector<Variant>& candidates,
                           const AlleleLogPosteriorMap& allele_log_posteriors,
                           const double min_log_posterior)
{
    assert(!allele_log_posteriors.empty());
    
    VariantCalls result {};
    result.reserve(candidates.size());
    
    for (const Variant& candidate : candidates) {
        const auto log_posterior = allele_log_posteriors.at(candidate.get_alt_allele());
        
        if (log_posterior >= min_log_posterior) {
            result.emplace_back(candidate, log_posterior);
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

double max_posterior(const VariantCalls& calls)
{
    return std::max_element(std::cbegin(calls), std::cend(calls),
                            [] (const VariantCall& lhs, const VariantCall& rhs) {
                                return lhs.log_posterior < rhs.log_posterior;
                            })->log_posterior;
}

VariantCallBlocks block_variant_calls(std::vector<VariantCalls>&& segmented_calls)
{
    VariantCallBlocks result {};
    result.reserve(segmented_calls.size());
    
    std::for_each(std::make_move_iterator(std::begin(segmented_calls)),
                  std::make_move_iterator(std::end(segmented_calls)),
                  [&result] (VariantCalls&& calls) {
                      std::vector<Variant> variants {};
                      variants.reserve(calls.size());
                      
                      std::transform(std::make_move_iterator(std::begin(calls)),
                                     std::make_move_iterator(std::end(calls)),
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
                           const AlleleLogPosteriorMap& allele_posteriors,
                           const Genotype<Haplotype>& genotype_call,
                           const double min_posterior)
{
    auto calls = call_variants(candidates, allele_posteriors, min_posterior);
    
    auto it = std::remove_if(std::begin(calls), std::end(calls),
                             [&genotype_call] (const VariantCall& call) {
                                 return !contains_exact(genotype_call, call.variant.get_alt_allele());
                             });
    calls.erase(it, std::end(calls));
    
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

// variant genotype calling

auto call_genotype(const GenotypeLogProbabilityVector& genotype_posteriors)
{
    return std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                            [] (const auto& lhs, const auto& rhs) {
                                return lhs.second < rhs.second;
                            })->first;
}

double marginalise(const Genotype<Allele>& genotype,
                   const GenotypeLogProbabilityVector& genotype_log_posteriors)
{
    std::vector<double> logs {};
    logs.reserve(genotype_log_posteriors.size());
    
    for (const auto& p : genotype_log_posteriors) {
        if (contains(p.first.get(), genotype)) {
            logs.push_back(p.second);
        }
    }
    
    return Maths::log_sum_exp(logs);
}

GenotypeCalls call_genotypes(const Genotype<Haplotype>& genotype_call,
                             const GenotypeLogProbabilityVector& genotype_log_posteriors,
                             const std::vector<GenomicRegion>& variant_regions)
{
    GenotypeCalls result {};
    result.reserve(variant_regions.size());
    
    for (const auto& region : variant_regions) {
        auto spliced_genotype = splice<Allele>(genotype_call, region);
        
        const auto posterior = marginalise(spliced_genotype, genotype_log_posteriors);
        
        result.emplace_back(std::move(spliced_genotype), posterior);
    }
    
    return result;
}

// Call cleanup

void remove_nongenotyped_calls(VariantCallBlock& block,
                               const GenotypeCall& genotype_call)
{
    const auto& genotype = genotype_call.genotype;
    auto it = std::remove_if(std::begin(block.variants), std::end(block.variants),
                             [&genotype] (const Variant& call) {
                                 return !genotype.contains(call.get_alt_allele());
                             });
    block.variants.erase(it, std::end(block.variants));
}

void remove_nongenotyped_calls(VariantCallBlocks& variant_calls,
                               const GenotypeCalls& genotype_calls)
{
    auto it = std::cbegin(genotype_calls);
    for (auto& block : variant_calls) {
        if (block.variants.size() > 1) {
            remove_nongenotyped_calls(block, *it++);
        }
    }
}

void parsimonise_variant_calls(VariantCallBlocks& variant_calls, const ReferenceGenome& reference)
{
    for (auto& segment_calls : variant_calls) {
        segment_calls.variants = parsimonise_together(segment_calls.variants, reference);
    }
}

    void parsimonise(VariantCallBlocks& variant_calls, GenotypeCalls& genotype_calls,
                     const ReferenceGenome& reference)
    {
        
    }

// reference genotype calling

double marginalise_reference_genotype(const Allele& reference_allele,
                                      const GenotypeLogProbabilityVector& genotype_log_posteriors)
{
    std::vector<double> logs {};
    logs.reserve(genotype_log_posteriors.size());
    
    for (const auto& p : genotype_log_posteriors) {
        if (is_homozygous(p.first.get(), reference_allele)) {
            logs.push_back(p.second);
        }
    }
    
    return Maths::log_sum_exp(logs);
}

RefCalls call_reference(const GenotypeLogProbabilityVector& genotype_log_posteriors,
                        const std::vector<Allele>& reference_alleles,
                        const ReadMap::mapped_type& reads, const double min_call_log_posterior)
{
    RefCalls result {};
    
    if (reference_alleles.empty()) return result;
    
    result.reserve(reference_alleles.size());
    
    for (const auto& reference_allele : reference_alleles) {
        double posterior {0};
        
        if (has_coverage(reads, mapped_region(reference_allele))) {
            posterior = marginalise_reference_genotype(reference_allele,
                                                       genotype_log_posteriors);
        }
        
        if (posterior >= min_call_log_posterior) {
            result.emplace_back(reference_allele, posterior);
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

unsigned count_alleles(const GenotypeCall& genotype_call)
{
    std::unordered_set<Allele> unique_alleles {};
    
    for (const auto& allele : genotype_call.genotype) {
        unique_alleles.emplace(allele);
    }
    
    return static_cast<unsigned>(unique_alleles.size());
}

auto count_alleles(const std::vector<Variant>& variants, const GenotypeCall& genotype_call)
{
    std::unordered_map<Allele, unsigned> allele_counts {};
    allele_counts.reserve(variants.size());
    
    for (const auto& allele : genotype_call.genotype) {
        ++allele_counts[allele];
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
    
float to_phred(const double log_posterior)
{
    return Maths::round(-log_posterior / Maths::ln_10_div_10, 2);
}

VcfRecord::Builder output_variant_call(const SampleIdType& sample,
                                       VariantCallBlock&& block,
                                       GenotypeCall&& genotype_call,
                                       const ReferenceGenome& reference,
                                       const ReadMap& reads,
                                       const bool sites_only)
{
    assert(!block.variants.empty());
    
    using std::to_string;
    
    auto result = VcfRecord::Builder {};
    
    const auto phred_quality = to_phred(block.log_posterior);
    
    const auto& reference_allele = block.variants.front().get_ref_allele();
    
    const auto& region = mapped_region(reference_allele);
    
    result.set_chromosome(contig_name(region));
    result.set_position(region_begin(region));
    result.set_ref_allele(reference_allele.get_sequence());
    result.set_alt_alleles(extract_alt_allele_sequences(block.variants));
    result.set_quality(phred_quality);
    
    result.add_info("AC",  to_strings(count_alleles(block.variants, genotype_call)));
    result.add_info("AN",  to_string(count_alleles(genotype_call)));
    result.add_info("NS",  to_string(count_samples_with_coverage(reads, region)));
    result.add_info("DP",  to_string(sum_max_coverages(reads, region)));
    result.add_info("SB",  Octopus::to_string(strand_bias(reads, region), 2));
    result.add_info("BQ",  to_string(static_cast<unsigned>(rmq_base_quality(reads, region))));
    result.add_info("MQ",  to_string(static_cast<unsigned>(rmq_mapping_quality(reads, region))));
    result.add_info("MQ0", to_string(count_mapq_zero(reads, region)));
    
    if (!sites_only) {
        result.set_format({"GT", "FT", "GQ", "PS", "PQ", "DP", "BQ", "MQ"});
        
        set_vcf_genotype(sample, genotype_call, block, result);
        
        result.add_genotype_field(sample, "FT", "."); // TODO
        result.add_genotype_field(sample, "GQ", Octopus::to_string(to_phred(genotype_call.log_posterior), 2));
        result.add_genotype_field(sample, "PS", to_string(region_begin(genotype_call.phase_region) + 1));
        result.add_genotype_field(sample, "PQ", Octopus::to_string(Maths::probability_to_phred<float>(genotype_call.phase_score), 2)); // TODO
        result.add_genotype_field(sample, "DP", to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result;
}

VcfRecord::Builder output_reference_call(const SampleIdType& sample,
                                         RefCall call, const ReferenceGenome& reference,
                                         const ReadMap& reads, unsigned ploidy,
                                         const bool sites_only)
{
    using std::to_string;
    
    auto result = VcfRecord::Builder {};
    
    const auto phred_quality = to_phred(call.log_posterior);
    
    const auto& region = mapped_region(call.reference_allele);
    
    result.set_chromosome(contig_name(region));
    result.set_position(region_begin(region));
    result.set_ref_allele(call.reference_allele.get_sequence().front());
    
    result.set_refcall();
    
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
    
    if (!sites_only) {
        result.set_format({"GT", "GQ", "DP", "BQ", "MQ"});
        
        result.add_homozygous_ref_genotype(sample, ploidy);
        result.add_genotype_field(sample, "GQ", Octopus::to_string(to_phred(call.log_posterior), 2));
        result.add_genotype_field(sample, "DP", to_string(max_coverage(reads.at(sample), region)));
        result.add_genotype_field(sample, "BQ", to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
        result.add_genotype_field(sample, "MQ", to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
    }
    
    return result;
}

void set_phasings(GenotypeCalls& variant_genotype_calls,
                  const Phaser::PhaseSet::SamplePhaseRegions& phase_set,
                  const std::vector<GenomicRegion>& called_regions)
{
    auto region_itr = std::cbegin(called_regions);
    
    for (auto& g : variant_genotype_calls) {
        const auto& call_region = *region_itr;
        
        const auto& phase = find_phase_region(phase_set, call_region);
        
        if (phase) {
            const auto overlapped = overlap_range(called_regions, phase->get().get_region(),
                                                  BidirectionallySortedTag {});
            
            assert(!overlapped.empty());
            
            g.phase_region = overlapped.front();
            g.phase_score  = phase->get().score;
        } else {
            g.phase_region = call_region;
            g.phase_score  = 0;
        }
        
        ++region_itr;
    }
}

std::vector<VcfRecord::Builder>
merge_calls(const SampleIdType& sample,
            VariantCallBlocks&& variant_calls, GenotypeCalls&& variant_genotype_calls,
            RefCalls&& refcalls, const ReferenceGenome& reference, const ReadMap& reads,
            const unsigned ploidy, const bool sites_only)
{
    using std::begin; using std::end; using std::make_move_iterator; using std::move;
    
    std::vector<VcfRecord::Builder> result {};
    result.reserve(variant_calls.size() + refcalls.size());
    
    merge_transform(make_move_iterator(begin(variant_calls)), make_move_iterator(end(variant_calls)),
                    make_move_iterator(begin(variant_genotype_calls)),
                    make_move_iterator(begin(refcalls)), make_move_iterator(end(refcalls)),
                    std::back_inserter(result),
                    [&] (VariantCallBlock&& variant_call, GenotypeCall&& genotype_call) {
                        return output_variant_call(sample, move(variant_call), move(genotype_call),
                                                   reference, reads, sites_only);
                    },
                    [&] (RefCall&& refcall) {
                        return output_reference_call(sample, move(refcall), reference, reads,
                                                     ploidy, sites_only);
                    });
    
    return result;
}
} // namespace

std::vector<VcfRecord::Builder>
IndividualVariantCaller::call_variants(const std::vector<Variant>& candidates,
                                       const std::vector<Allele>& callable_alleles,
                                       CallerLatents* latents,
                                       const Phaser::PhaseSet& phase_set,
                                       const ReadMap& reads) const
{
    assert(!callable_alleles.empty());
    
    const auto dlatents = dynamic_cast<Latents*>(latents);
    
    const auto genotype_log_posteriors = zip(dlatents->genotypes_, dlatents->genotype_log_posteriors_);
    
    if (TRACE_MODE) {
        Logging::TraceLogger log {};
        debug::print_genotype_posteriors(stream(log), genotype_log_posteriors, -1);
    } else if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        debug::print_genotype_posteriors(stream(log), genotype_log_posteriors);
    }
    
    const auto allele_posteriors = compute_allele_posteriors(genotype_log_posteriors, callable_alleles);
    
    if (TRACE_MODE) {
        Logging::TraceLogger log {};
        debug::print_allele_posteriors(stream(log), allele_posteriors, -1);
    } else if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        debug::print_allele_posteriors(stream(log), allele_posteriors);
    }
    
    const auto genotype_call = call_genotype(genotype_log_posteriors);
    
    auto variant_calls = call_blocked_variants(candidates, allele_posteriors,
                                               genotype_call, min_variant_log_posterior_);
    
    //debug::print_variant_calls(variant_calls);
    
    parsimonise_variant_calls(variant_calls, reference_);
    
    auto called_regions = extract_regions(variant_calls);
    
    auto variant_genotype_calls = call_genotypes(genotype_call, genotype_log_posteriors, called_regions);
    
    set_phasings(variant_genotype_calls, phase_set.phase_regions.at(sample_), called_regions); // TODO
    
    //debug::print_genotype_calls(variant_genotype_calls);
    
    auto candidate_ref_alleles = generate_candidate_reference_alleles(callable_alleles, called_regions,
                                                                      candidates, refcall_type_);
    
    auto refcalls = call_reference(genotype_log_posteriors, candidate_ref_alleles, reads.at(sample_),
                                   min_refcall_log_posterior_);
    
    return merge_calls(sample_, std::move(variant_calls), std::move(variant_genotype_calls),
                       std::move(refcalls), reference_, reads, ploidy_, call_sites_only_);
}

namespace debug
{
    template <typename S>
    void print_genotype_posteriors(S&& stream,
                                   const GenotypeLogProbabilityVector& genotype_posteriors,
                                   const std::size_t n)
    {
        const auto m = std::min(n, genotype_posteriors.size());
        
        if (m == genotype_posteriors.size()) {
            stream << "Printing all genotype log posteriors " << '\n';
        } else {
            stream << "Printing top " << m << " genotype log posteriors " << '\n';
        }
        
        using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
        
        std::vector<std::pair<GenotypeReference, double>> v {};
        v.reserve(genotype_posteriors.size());
        
        std::copy(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                  std::back_inserter(v));
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&] (const auto& p) {
                          ::debug::print_variant_alleles(stream, p.first);
                          stream << " " << p.second << '\n';
                      });
    }
    
    void print_genotype_posteriors(const GenotypeLogProbabilityVector& genotype_posteriors,
                                   const std::size_t n)
    {
        print_genotype_posteriors(std::cout, genotype_posteriors, n);
    }
    
    template <typename S>
    void print_allele_posteriors(S&& stream, const AlleleLogPosteriorMap& allele_posteriors,
                                 const std::size_t n)
    {
        const auto m = std::min(n, allele_posteriors.size());
        
        if (m == allele_posteriors.size()) {
            stream << "Printing all allele posteriors " << '\n';
        } else {
            stream << "Printing top " << m << " allele posteriors " << '\n';
        }
        
        using AlleleReference = std::reference_wrapper<const Allele>;
        
        std::vector<std::pair<AlleleReference, double>> v {};
        v.reserve(allele_posteriors.size());
        
        std::copy(std::cbegin(allele_posteriors), std::cend(allele_posteriors),
                  std::back_inserter(v));
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&] (const auto& p) {
                          stream << p.first.get() << " " << p.second << '\n';
                      });
    }
    
    void print_allele_posteriors(const AlleleLogPosteriorMap& allele_posteriors,
                                 const std::size_t n)
    {
        print_allele_posteriors(std::cout, allele_posteriors, n);
    }
} // namespace debug
} // namespace Octopus
