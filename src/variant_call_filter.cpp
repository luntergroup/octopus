//
//  variant_call_filter.cpp
//  Octopus
//
//  Created by Daniel Cooke on 31/05/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "variant_call_filter.hpp"

#include <unordered_map>
#include <map>
#include <limits>
#include <numeric>
#include <cmath>

#include "common.hpp"

#include "vcf_reader.hpp"
#include "vcf_writer.hpp"
#include "vcf_header.hpp"
#include "vcf_record.hpp"
#include "variant.hpp"

#include "genomic_region.hpp"
#include "mappable_flat_set.hpp"
#include "mappable_map.hpp"
#include "aligned_read.hpp"
#include "read_utils.hpp"

#include "genotype_reader.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "maths.hpp"
#include "string_utils.hpp"

#include <boost/math/distributions/hypergeometric.hpp>

namespace Octopus
{
VariantCallFilter::VariantCallFilter(const ReferenceGenome& reference, const ReadPipe& read_pipe)
:
reference_ {reference},
read_pipe_ {read_pipe}
{}

namespace
{
    auto mapped_region(const VcfRecord& call)
    {
        using SizeType = GenomicRegion::SizeType;
        const auto begin = call.pos() - 1;
        return GenomicRegion {call.chrom(), begin, begin + static_cast<SizeType>(call.ref().size())};
    }
    
    auto mapped_regions(const std::vector<VcfRecord>& calls)
    {
        std::vector<GenomicRegion> result {};
        result.reserve(calls.size());
        
        std::transform(std::cbegin(calls), std::cend(calls), std::back_inserter(result),
                       [] (const auto& call) { return mapped_region(call); });
        
        std::sort(std::begin(result), std::end(result));
        
        return result;
    }
    
    auto encompassing_region(const std::vector<VcfRecord>& calls)
    {
        return encompassing_region(mapped_regions(calls));
    }
    
    auto fetch_reads(const std::vector<VcfRecord>& calls, const ReadPipe& read_pipe)
    {
        return read_pipe.fetch_reads(encompassing_region(calls));
    }
    
    using HaplotypeSupportMap = std::multimap<Haplotype, AlignedRead>;
    
    using HaplotypeLikelihoods = std::vector<std::vector<double>>;
    
    auto max_posterior_haplotypes(const Genotype<Haplotype>& genotype, const unsigned read,
                                  const HaplotypeLikelihoods& likelihoods)
    {
        std::vector<unsigned> result {};
        result.reserve(genotype.ploidy());
        
        auto max_likelihood = std::numeric_limits<double>::lowest();
        
        for (unsigned k {0}; k < genotype.ploidy(); ++k) {
            const auto curr = likelihoods[k][read];
            
            if (Maths::almost_equal(curr, max_likelihood)) {
                result.push_back(k);
            } else if (curr > max_likelihood) {
                result.assign({k});
                max_likelihood = curr;
            }
        }
        
        result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
        
        return result;
    }
    
    auto calculate_support(const Genotype<Haplotype>& genotype, const ReadContainer& reads,
                           const HaplotypeLikelihoods& likelihoods)
    {
        HaplotypeSupportMap result {};
        
        //std::cout << mapped_region(genotype) << std::endl;
        for (unsigned i {0}; i < reads.size(); ++i) {
            //std::cout << reads[i].mapped_region() << " " << reads[i].cigar_string() << std::endl;
            
            const auto top = max_posterior_haplotypes(genotype, i, likelihoods);
            
            if (top.size() == 1) {
                result.emplace(genotype[top.front()], reads[i]);
            }
        }
        
        return result;
    }
    
    auto expand(const Genotype<Haplotype>& genotype, Haplotype::SizeType n)
    {
        Genotype<Haplotype> result {genotype.ploidy()};
        
        for (const auto& haplotype : genotype) {
            result.emplace(expand(haplotype, n));
        }
        
        return result;
    }
    
    auto calculate_likelihoods(const Genotype<Haplotype>& genotype, const ReadContainer& reads)
    {
        ReadMap tmp {};
        
        tmp["sample"] = reads;
        
        const auto& genotype_region = mapped_region(genotype);
        
        const auto reads_region = encompassing_region(reads);
        
        unsigned min_lhs_expansion {0};
        
        if (begins_before(reads_region, genotype_region)) {
            min_lhs_expansion += begin_distance(reads_region, genotype_region);
        }
        
        unsigned min_rhs_expansion {0};
        
        if (ends_before(genotype_region, reads_region)) {
            min_rhs_expansion += end_distance(genotype_region, reads_region);
        }
        
        const auto min_expansion = std::max({min_lhs_expansion, min_rhs_expansion, 20u});
        
        std::map<Haplotype, Haplotype> expanded_haplotypes {};
        
        for (const auto& haplotype : genotype) {
            expanded_haplotypes.emplace(haplotype, expand(haplotype, min_expansion));
        }
        
        const auto expanded_genotype = expand(genotype, min_expansion);
        
        const auto haplotypes = expanded_genotype.copy_unique();
        
//        for (const auto& haplotype : haplotypes) {
//            std::cout << haplotype.mapped_region() << std::endl;
//            ::debug::print_variant_alleles(haplotype);
//            std::cout << '\n';
//        }
//        
//        std::cout << reads[48].mapped_region() << " " << reads[48].cigar_string() << std::endl;
        
        HaplotypeLikelihoodCache likelihoods {static_cast<unsigned>(haplotypes.size()), {"sample"}};
        
        likelihoods.populate(tmp, haplotypes);
        
        HaplotypeLikelihoods result(genotype.ploidy());
        
        std::transform(std::cbegin(genotype), std::cend(genotype), std::begin(result),
                       [&likelihoods, &expanded_haplotypes] (const auto& haplotype) {
                           return likelihoods.log_likelihoods("sample", expanded_haplotypes.at(haplotype));
                       });
        
        return result;
    }
    
    auto calculate_support(const Genotype<Haplotype>& genotype, const ReadContainer& reads)
    {
        const auto likelihoods = calculate_likelihoods(genotype, reads);
        return calculate_support(genotype, reads, likelihoods);
    }
    
    decltype(auto) find_genotype(const GenotypeMap& genotypes, const SampleIdType& sample,
                                 const GenomicRegion& region)
    {
        return genotypes.at(sample).overlap_range(region).front();
    }
    
    using VariantMap = std::multimap<SampleIdType, Variant>;
    
    auto extract_variants(const VcfRecord& call, const std::vector<SampleIdType>& samples)
    {
        VariantMap result {};
        
        const auto region = mapped_region(call);
        
        const auto& alt_alleles = call.alt();
        
        for (const auto& sample : samples) {
            auto genotype = call.get_sample_value(sample, "GT");
            
            std::sort(std::begin(genotype), std::end(genotype));
            
            genotype.erase(std::unique(std::begin(genotype), std::end(genotype)), std::end(genotype));
            
            for (const auto& allele : genotype) {
                const auto it = std::find(std::cbegin(alt_alleles), std::cend(alt_alleles), allele);
                
                if (it != std::cend(alt_alleles)) {
                    result.emplace(std::piecewise_construct,
                                   std::forward_as_tuple(sample),
                                   std::forward_as_tuple(region, call.ref(), *it));
                }
            }
        }
        
        return result;
    }
    
    using VariantSupportMap = std::multimap<Variant, AlignedRead>;
    
    template <typename Range>
    auto calculate_variant_support(const Genotype<Haplotype>& genotype, const ReadContainer& reads,
                                   Range variants)
    {
        const auto haplotype_support = calculate_support(genotype, reads);
        
        VariantSupportMap result {};
        
        std::for_each(variants.first, variants.second,
                      [&haplotype_support, &result] (const auto& p) {
                          for (const auto& hp : haplotype_support) {
                              if (overlaps(p.second, hp.second)
                                  && hp.first.contains(p.second.alt_allele())) {
                                        result.emplace(p.second, hp.second);
                              }
                          }
                      });
        
        return result;
    }
    
    auto calculate_alt_frequency_credible_region(const VariantSupportMap& support, const Variant& v,
                                                 const ReadContainer& reads)
    {
        const auto er = support.equal_range(v);
        
        const auto num_alt = std::distance(er.first, er.second);
        
        const auto num_ref = count_overlapped(reads, mapped_region(v)) - num_alt;
        
        return Maths::beta_hdi(num_alt + 0.5, num_ref + 0.5, 0.9999);
    }
    
    struct ReadDirectionCounts
    {
        unsigned forward, reverse;
    };
    
    auto calculate_direction_counts(const VariantSupportMap& support, const Variant& v)
    {
        const auto er = support.equal_range(v);
        
        const auto c = std::count_if(er.first, er.second,
                                     [] (const auto& p) {
                                         return p.second.is_marked_reverse_mapped();
                                     });
        
        return ReadDirectionCounts {
            static_cast<unsigned>(std::distance(er.first, er.second) - c),
            static_cast<unsigned>(c)
        };
    }
    
    double fisher_exact(const unsigned a, const unsigned b,
                        const unsigned c, const unsigned d)
    {
        if (a + b == 0 || c + d == 0) return 1.0;
        
        const auto r = a + c;
        const auto n = c + d;
        
        if (r == 0 || n == 0) return 1.0;
        
        const auto N = a + b + c + d;
        
        const boost::math::hypergeometric_distribution<> dist {r, n, N};
        
        const auto cutoff = pdf(dist, c);
        
        const auto min_k = static_cast<unsigned>(std::max(0, static_cast<int>(r + n - N)));
        
        double result {0};
        
        for (auto k = min_k, max_k = std::min(r, n); k <= max_k; ++k) {
            const auto p = pdf(dist, k);
            if (p <= cutoff) result += p;
        }
        
        return result;
    }
    
    auto calculate_strand_bias(const VariantSupportMap& support, const Variant& v,
                               const ReadContainer& reads)
    {
        const auto alt_counts = calculate_direction_counts(support, v);
        
        ReadDirectionCounts other_counts {
            static_cast<unsigned>(count_forward(reads, mapped_region(v))) - alt_counts.forward,
            static_cast<unsigned>(count_reverse(reads, mapped_region(v))) - alt_counts.reverse
        };
        
        return fisher_exact(alt_counts.forward, alt_counts.reverse,
                            other_counts.forward, other_counts.reverse);
    }
    
    auto convert_to_probabilities(const std::vector<unsigned>& counts)
    {
        const auto norm = std::accumulate(std::cbegin(counts), std::cend(counts), 0.0);
        
        std::vector<double> result(counts.size());
        
        std::transform(std::cbegin(counts), std::cend(counts), std::begin(result),
                       [norm] (const auto c) { return static_cast<double>(c) / norm; });
        
        return result;
    }
    
    auto kl_divergence(const std::vector<double>& p, const std::vector<double>& q)
    {
        return std::inner_product(std::cbegin(p), std::cend(p), std::cbegin(q), 0.0,
                                  std::plus<> {}, [] (const auto a, const auto b) {
                                      return a * std::log(a / b);
                                  });
    }
    
    auto symmetric_kl_divergence(const std::vector<double>& p, const std::vector<double>& q)
    {
        return kl_divergence(p, q) + kl_divergence(q, p);
    }
    
    auto calculate_mq_bias(const VariantSupportMap& support, const Variant& v,
                           const ReadContainer& reads)
    {
        using Q = AlignedRead::QualityType;
        
        static constexpr unsigned Max_qualities {std::numeric_limits<Q>::max() - std::numeric_limits<Q>::min()};
        
        std::vector<unsigned> variant_mq_histogram(Max_qualities, 0);
        
        const auto er = support.equal_range(v);
        
        std::for_each(er.first, er.second,
                      [&variant_mq_histogram] (const auto& p) {
                          ++variant_mq_histogram[p.second.mapping_quality()];
                      });
        
        std::vector<unsigned> other_mq_histogram(Max_qualities, 0);
        
        for (const auto& read : reads) {
            ++other_mq_histogram[read.mapping_quality()];
        }
        
        std::transform(std::cbegin(other_mq_histogram), std::cend(other_mq_histogram),
                       std::cbegin(variant_mq_histogram), std::begin(other_mq_histogram),
                       [] (const auto a, const auto b) { return a - b; });
        
        if (std::accumulate(std::cbegin(variant_mq_histogram), std::cend(variant_mq_histogram), 0) < 5) {
            return 0.0;
        }
        
        // Prevent 0 counts (hence 0 probability), so KL is well defined
        for (auto& e : variant_mq_histogram) ++e;
        for (auto& e : other_mq_histogram) ++e;
        
        const auto variant_dist = convert_to_probabilities(variant_mq_histogram);
        const auto other_dist   = convert_to_probabilities(other_mq_histogram);
        
        return symmetric_kl_divergence(variant_dist, other_dist);
    }
    
    auto get_dummy_model_posterior(const VcfRecord& call)
    {
        return std::stod(call.info_value("DMBF").front());
    }
} // namespace

void VariantCallFilter::filter(const VcfReader& source, VcfWriter& dest, const RegionMap& regions) const
{
    const auto header = source.fetch_header();
    
    if (!dest.is_header_written()) {
        VcfHeader::Builder hb {header};
        
        hb.add_filter("AB", "Alt allele frequency is less than expected");
        hb.add_filter("MODEL", "The caller specific model filter failed");
        
        dest << hb.build_once();
    }
    
    const auto samples = header.samples();
    
    for (const auto& p : regions) {
        for (const auto& region : p.second) {
            auto calls = source.fetch_records(region);
            
            if (calls.empty()) continue;
            
            const auto reads = fetch_reads(calls, read_pipe_);
            
            const auto reads_region = encompassing_region(reads);
            
            const auto genotypes = extract_genotypes(calls, header, reference_, reads_region);
            
            for (auto& call : calls) {
                //std::cout << call << std::endl;
                
                const auto call_region = mapped_region(call);
                
                const auto call_reads = copy_overlapped(reads, call_region);
                
                VcfRecord::Builder cb {call};
                
                cb.set_info("NS",  count_samples_with_coverage(call_reads));
                cb.set_info("DP",  sum_max_coverages(call_reads));
                cb.set_info("SB",  Octopus::to_string(strand_bias(call_reads), 2));
                cb.set_info("BQ",  static_cast<unsigned>(rmq_base_quality(call_reads)));
                cb.set_info("MQ",  static_cast<unsigned>(rmq_mapping_quality(call_reads)));
                cb.set_info("MQ0", count_mapq_zero(call_reads));
                
                bool filtered {false};
                
                const auto dummy_model_posterior = get_dummy_model_posterior(call);
                
                cb.clear_info("DMBF");
                
                if (dummy_model_posterior > 0.1) {
                    cb.add_filter("MODEL");
                    filtered = true;
                }
                
                if (call.qual() && *call.qual() < 10) {
                    cb.add_filter("q10");
                    filtered = true;
                }
                
                const auto rmq = std::ceil(rmq_mapping_quality(call_reads));
                
                if (rmq < 40) {
                    cb.add_filter("MQ");
                    filtered = true;
                }
                
                const auto variants = extract_variants(call, samples);
                
                bool strand_biased {true};
                bool sample_rmq_failed {false};
                bool sample_kl_failed {false};
                bool sample_allele_biased {false};
                bool all_homozygous {true};
                
                for (const auto& sample : samples) {
                    //std::cout << sample << std::endl;
                    
                    const auto& sample_call_reads = call_reads.at(sample);
                    
                    if (rmq >= 40) {
                        const auto sample_rmq = std::ceil(rmq_mapping_quality(sample_call_reads));
                        
                        if (sample_rmq < 40) {
                            sample_rmq_failed = true;
                        }
                    }
                    
                    const auto& genotype = find_genotype(genotypes, sample, call_region);
                    
                    const auto sample_variants = variants.equal_range(sample);
                    
                    if (sample_variants.first == sample_variants.second) continue;
                    
                    std::vector<ReadDirectionCounts> count {};
                    
                    if (call.is_heterozygous(sample)) {
                        all_homozygous = false;
                        
                        const auto variant_support = calculate_variant_support(genotype, sample_call_reads,
                                                                               sample_variants);
                        
//                        const auto cr = calculate_alt_frequency_credible_region(variant_support, sample_variants.first->second,
//                                                                                sample_call_reads);
//                        
//                        if (cr.second < 0.5) {
//                            sample_allele_biased = true;
//                        }
                        
                        const auto pval = calculate_strand_bias(variant_support, sample_variants.first->second,
                                                                sample_call_reads);
                        
                        if (pval > 0.00005) {
                            strand_biased = false;
                        }
                        
                        const auto mq_bias = calculate_mq_bias(variant_support, sample_variants.first->second,
                                                               sample_call_reads);
                        
                        if (mq_bias > 0.5) {
                            sample_kl_failed = true;
                        }
                    }
                    
                    cb.set_format(sample, "DP", max_coverage(call_reads.at(sample)));
                    cb.set_format(sample, "BQ", static_cast<unsigned>(rmq_base_quality(call_reads.at(sample))));
                    cb.set_format(sample, "MQ", static_cast<unsigned>(rmq_mapping_quality(call_reads.at(sample))));
                }
                
                if (sample_rmq_failed) {
                    cb.add_filter("MQ");
                    filtered = true;
                }
                
                if (!all_homozygous) {
                    if (sample_allele_biased) {
                        cb.add_filter("AB");
                        filtered = true;
                    }
                    
                    if (sample_kl_failed) {
                        cb.add_filter("KL");
                        filtered = true;
                    }
                    
                    if (strand_biased) {
                        cb.add_filter("SB");
                        filtered = true;
                    }
                }
                
                if (!filtered) {
                    cb.set_passed();
                }
                
                call = cb.build();
            }
            
//            const auto it = std::remove_if(std::begin(calls), std::end(calls),
//                                           [] (const auto& call) {
//                                               return !is_somatic(call);
//                                           });
//            calls.erase(it, std::end(calls));
            
            write(calls, dest);
        }
    }
}
} // namespace Octopus
