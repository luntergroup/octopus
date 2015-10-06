//
//  cancer_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 16/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_caller.hpp"

#include <unordered_map>
#include <algorithm> // std::max_element, std::any_of
#include <numeric>

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
#include "maths.hpp"
#include "cancer_genotype.hpp"
#include "genotype_model.hpp"
#include "cancer_genotype_model.hpp"
#include "read_utils.hpp"
#include "string_utils.hpp"

#include "haplotype_prior_model.hpp"

#include <iostream> // TEST

namespace Octopus
{
// public methods
    
    CancerVariantCaller::CancerVariantCaller(ReferenceGenome& reference, CandidateVariantGenerator& candidate_generator,
                                             RefCallType refcall_type, double min_variant_posterior,
                                             double min_refcall_posterior, const SampleIdType& normal_sample)
    :
    VariantCaller {reference, candidate_generator, refcall_type},
    normal_sample_ {normal_sample},
    min_variant_posterior_ {min_variant_posterior},
    min_refcall_posterior_ {min_refcall_posterior}
    {}
    
    std::string CancerVariantCaller::do_get_details() const
    {
        return "cancer caller. normal sample = " + normal_sample_;
    }

    GenomicRegion CancerVariantCaller::get_init_region(const GenomicRegion& region, const ReadMap& reads,
                                                       const std::vector<Variant>& candidates)
    {
        return region;
    }
    
    GenomicRegion CancerVariantCaller::get_next_region(const GenomicRegion& current_region, const ReadMap& reads,
                                                       const std::vector<Variant>& candidates)
    {
        return GenomicRegion {"TEST", 0, 0};
    }
    
    // private methods
    
    using GenotypeCalls = std::unordered_map<Octopus::SampleIdType, std::pair<CancerGenotype<Allele>, double>>;
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
    
    std::ostream& operator<<(std::ostream& os, const std::array<double, 3>& arr)
    {
        os << arr[0] << " " << arr[1] << " " << arr[2];
        return os;
    }
    
    std::ostream& operator<<(std::ostream& os, const std::unordered_map<Octopus::SampleIdType, std::array<double, 3>>& m)
    {
        for (const auto& p : m) os << p.first << ": " << p.second << "\n";
        return os;
    }
    
    static std::vector<GenomicRegion> get_segment_regions(const std::vector<std::vector<Allele>>& segments)
    {
        std::vector<GenomicRegion> result {};
        result.reserve(segments.size());
        for (const auto& segment : segments) result.push_back(segment.front().get_region());
        return result;
    }
    
    auto find_map_genotype(const GenotypeModel::Cancer::GenotypeProbabilities& genotype_posteriors)
    {
        return *std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                 [] (const auto& lhs, const auto& rhs) {
                                     return lhs.second < rhs.second;
                                 });
    }
    
    bool is_cancer(const GenotypeModel::Cancer::Latents& latents)
    {
        if (find_map_genotype(latents.genotype_posteriors).second > 0.95) {
            for (const auto& sample_genotype_weights : latents.genotype_weights) {
                if (sample_genotype_weights.second[2] > 0.2) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    bool is_homozygous_reference(const CancerGenotype<Haplotype>& genotype, const Haplotype& reference)
    {
        return genotype.get_normal_genotype().count(reference) == genotype.ploidy();
    }
    
    bool has_cancer_sample(const GenotypeModel::Cancer::GenotypeWeights& genotype_weights, double min_weight = 0.1)
    {
        return std::any_of(std::cbegin(genotype_weights), std::cend(genotype_weights),
                           [min_weight] (const auto& p) { return p.second[2] >= min_weight; });
    }
    
    unsigned count_alleles(const GenotypeCalls& genotype_calls)
    {
        if (genotype_calls.empty()) return 0;
        
        std::unordered_set<Allele> unique_alleles {};
        
        for (const auto& sample_call : genotype_calls) {
            for (const auto& allele : sample_call.second.first.get_normal_genotype()) {
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
            for (const auto& allele : sample_call.second.first.get_normal_genotype()) {
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
    
    std::vector<VcfRecord::SequenceType> to_vcf_genotype(const CancerGenotype<Allele>& genotype)
    {
        std::vector<VcfRecord::SequenceType> result {};
        result.reserve(genotype.ploidy());
        for (const auto& allele : genotype.get_normal_genotype()) result.push_back(allele.get_sequence());
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
    
    VcfRecord output_reference_call(RefCall call, ReferenceGenome& reference, const ReadMap& reads)
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
            result.add_homozygous_ref_genotype(sample, 2);
            result.add_genotype_field(sample, "GP", std::to_string(Maths::probability_to_phred(sample_posteior.second)));
            result.add_genotype_field(sample, "DP", std::to_string(max_coverage(reads.at(sample), region)));
            result.add_genotype_field(sample, "BQ", std::to_string(static_cast<unsigned>(rmq_base_quality(reads.at(sample), region))));
            result.add_genotype_field(sample, "MQ", std::to_string(static_cast<unsigned>(rmq_mapping_quality(reads.at(sample), region))));
        }
        
        return result.build_once();
    }
    
    std::vector<VcfRecord>
    CancerVariantCaller::call_variants(const GenomicRegion& region, const std::vector<Variant>& candidates,
                                       const ReadMap& reads)
    {
        std::vector<VcfRecord> result {};
        
        Octopus::HaplotypeTree tree {reference_};
        extend_tree(candidates, tree);
        
        auto haplotypes = tree.get_haplotypes(region);
        
        std::cout << "there are " << haplotypes.size() << " haplotypes" << std::endl;
        
        GenotypeModel::Cancer genotype_model {normal_sample_};
        
        auto latents = genotype_model.evaluate(haplotypes, reads, reference_);
        
        auto map_genotype = find_map_genotype(latents.genotype_posteriors);
        
        if (map_genotype.second >= min_variant_posterior_) {
            Haplotype reference_haplotype {reference_, region};
            
            auto has_variant = !is_homozygous_reference(map_genotype.first, reference_haplotype);
            auto has_cancer  = has_cancer_sample(latents.genotype_weights);
            
            if (has_variant || has_cancer) {
                if (has_variant && has_cancer) {
                    std::cout << "somatic mutation somewhere" << std::endl;
                } else if (has_variant && !has_cancer) {
                    std::cout << "probable normal variation" << std::endl;
                } else {
                    std::cout << "wierdness" << std::endl;
                }
            }
            
//            auto alleles  = decompose(candidates);
//            auto segments = segment(alleles);
//            auto regions  = get_segment_regions(segments);
        }
        
        print_variant_alleles(map_genotype.first);
        std::cout << std::endl << std::endl;
        std::cout << latents.genotype_weights << std::endl;
        std::cout << std::boolalpha << is_cancer(latents) << std::endl;
        
        return result;
    }

} // namespace Octopus
