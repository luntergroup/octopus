//
//  population_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_genotype_model.hpp"

#include <algorithm> // std::transform
#include <cmath>

#include "read_model.hpp"
#include "maths.hpp"
#include "read_utils.hpp"

#include <iostream> // TEST

namespace Octopus
{
    namespace GenotypeModel
    {
    // public methods
    
    Population::Population(unsigned ploidy, unsigned max_em_iterations, double em_epsilon)
    :
    ploidy_ {ploidy},
    max_em_iterations_ {max_em_iterations},
    em_epsilon_ {em_epsilon}
    {}
    
    // non member methods
    
    struct GenotypeLogProbability
    {
        GenotypeLogProbability() = delete;
        GenotypeLogProbability(const Genotype<Haplotype>& genotype, double log_probability)
        : genotype {genotype}, log_probability {log_probability} {}
        
        const Genotype<Haplotype>& genotype;
        double log_probability;
    };
        
    using GenotypeLogMarginals = std::vector<GenotypeLogProbability>;
    
    using GenotypeLogLikelihood = double;
    using SampleGenotypeLogLikelihoods = std::vector<GenotypeLogLikelihood>;
    using GenotypeLogLikelihoods = std::unordered_map<SampleIdType, SampleGenotypeLogLikelihoods>;
    
    using GenotypePosterior = double;
    using SampleGenotypePosteriors = std::vector<GenotypePosterior>;
    using GenotypePosteriors = std::unordered_map<SampleIdType, SampleGenotypePosteriors>;
    
    // DEBUG methods
    
    inline void print(const HaplotypePriorCounts& prior_counts)
    {
        for (const auto& h : prior_counts) {
            print_variant_alleles(h.first);
            std::cout << " " << h.second << std::endl;
        }
    }
    
//    inline void print(const GenotypeLikelihoods& likelihoods, size_t n = 5)
//    {
//        std::cout << "top " << n << " genotype likelihoods for each sample" << std::endl;
//        for (const auto& sample_likelihoods : likelihoods) {
//            std::cout << sample_likelihoods.first << ":" << std::endl;
//            std::vector<std::pair<Genotype<Haplotype>, double>> v {};
//            v.reserve(sample_likelihoods.second.size());
//            std::copy(std::cbegin(sample_likelihoods.second), std::cend(sample_likelihoods.second),
//                      std::back_inserter(v));
//            std::sort(std::begin(v), std::end(v), [] (const auto& lhs, const auto& rhs) {
//                return lhs.second > rhs.second;
//            });
//            for (unsigned i {}; i < std::min(n, v.size()); ++i) {
//                print_variant_alleles(v[i].first);
//                std::cout << " " << v[i].second << std::endl;
//            }
//        }
//    }
//
//    void print_top_haplotypes(const std::vector<Haplotype>& haplotypes, const ReadMap& reads, size_t n = 3)
//    {
//        std::cout << "top " << n << " haplotype likelihoods for each read in each sample" << std::endl;
//        auto m = std::min(n, haplotypes.size());
//        for (const auto& sample_reads : reads) {
//            std::cout << sample_reads.first << ":" << std::endl;
//            for (const auto& read : sample_reads.second) {
//                std::cout << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
//                std::vector<std::pair<Haplotype, double>> top {};
//                top.reserve(haplotypes.size());
//                for (const auto& haplotype : haplotypes) {
//                    top.emplace_back(haplotype, ReadModel(2).log_probability(read, haplotype));
//                }
//                std::sort(std::begin(top), std::end(top), [] (const auto& lhs, const auto& rhs) {
//                    return lhs.second > rhs.second;
//                });
//                for (unsigned i {}; i < m; ++i) {
//                    std::cout << "\t* ";
//                    print_variant_alleles(top[i].first);
//                    std::cout << " " << top[i].second << std::endl;
//                }
//            }
//        }
//    }
//    
//    void print_top_haplotypes(const std::vector<Haplotype>& haplotypes, const ReadMap& reads,
//                              const SampleIdType& sample, size_t n = 3)
//    {
//        std::cout << "top " << n << " haplotype likelihoods for each read in sample " << sample << std::endl;
//        auto m = std::min(n, haplotypes.size());
//        for (const auto& read : reads.at(sample)) {
//            std::cout << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
//            std::vector<std::pair<Haplotype, double>> top {};
//            top.reserve(haplotypes.size());
//            for (const auto& haplotype : haplotypes) {
//                top.emplace_back(haplotype, ReadModel(2).log_probability(read, haplotype));
//            }
//            std::sort(std::begin(top), std::end(top), [] (const auto& lhs, const auto& rhs) {
//                return lhs.second > rhs.second;
//            });
//            for (unsigned i {}; i < m; ++i) {
//                std::cout << "\t* ";
//                print_variant_alleles(top[i].first);
//                std::cout << " " << top[i].second << std::endl;
//            }
//        }
//    }
//    
//    void print_top_genotypes(const std::vector<Genotype<Haplotype>>& genotypes, const ReadMap& reads, size_t n = 3)
//    {
//        std::cout << "top " << n << " genotype likelihoods for each read in each sample" << std::endl;
//        auto m = std::min(n, genotypes.size());
//        for (const auto& sample_reads : reads) {
//            std::cout << sample_reads.first << ":" << std::endl;
//            for (const auto& read : sample_reads.second) {
//                std::cout << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
//                std::vector<std::pair<Genotype<Haplotype>, double>> top {};
//                top.reserve(genotypes.size());
//                for (const auto& genotype : genotypes) {
//                    top.emplace_back(genotype, ReadModel(2).log_probability(read, genotype));
//                }
//                std::sort(std::begin(top), std::end(top), [] (const auto& lhs, const auto& rhs) {
//                    return lhs.second > rhs.second;
//                });
//                for (unsigned i {}; i < m; ++i) {
//                    std::cout << "\t* ";
//                    print_variant_alleles(top[i].first);
//                    std::cout << " " << top[i].second << std::endl;
//                }
//            }
//        }
//    }
//    
//    void print_top_genotypes(const std::vector<Genotype<Haplotype>>& genotypes, const ReadMap& reads,
//                             const SampleIdType& sample, size_t n = 3)
//    {
//        std::cout << "top " << n << " genotype likelihoods for each read in sample " << sample << std::endl;
//        auto m = std::min(n, genotypes.size());
//        for (const auto& read : reads.at(sample)) {
//            std::cout << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
//            std::vector<std::pair<Genotype<Haplotype>, double>> top {};
//            top.reserve(genotypes.size());
//            for (const auto& genotype : genotypes) {
//                top.emplace_back(genotype, ReadModel(2).log_probability(read, genotype));
//            }
//            std::sort(std::begin(top), std::end(top), [] (const auto& lhs, const auto& rhs) {
//                return lhs.second > rhs.second;
//            });
//            for (unsigned i {}; i < m; ++i) {
//                std::cout << "\t* ";
//                print_variant_alleles(top[i].first);
//                std::cout << " " << top[i].second << std::endl;
//            }
//        }
//    }
    
    // non member methods
    
    GenotypeLogMarginals
    init_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                                const HaplotypeFrequencies& haplotype_frequencies)
    {
        GenotypeLogMarginals result {};
        result.reserve(genotypes.size());
        
        for (const auto& genotype : genotypes) {
            result.emplace_back(genotype, log_hardy_weinberg(genotype, haplotype_frequencies));
        }
        
        return result;
    }
    
    GenotypeLogLikelihoods
    compute_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                     const ReadMap& reads, ReadModel& read_model)
    {
        GenotypeLogLikelihoods result {};
        result.reserve(reads.size());
        
        for (const auto& sample_reads_p : reads) {
            const auto& sample       = sample_reads_p.first;
            const auto& sample_reads = sample_reads_p.second;
            
            SampleGenotypeLogLikelihoods sample_log_likelihoods(genotypes.size());
            
            std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(sample_log_likelihoods),
                           [&sample_reads, &read_model] (const auto& genotype) {
                               return read_model.log_probability(std::cbegin(sample_reads),
                                                                 std::cend(sample_reads), genotype);
                           });
            
            result.emplace(sample, std::move(sample_log_likelihoods));
        }
        
        return result;
    }
    
    static void normalise_exp(SampleGenotypePosteriors& unnormalised_log_posteriors)
    {
        const auto norm = Maths::log_sum_exp<double>(unnormalised_log_posteriors);
        for (auto& p : unnormalised_log_posteriors) {
            p -= norm;
            p = std::exp(p);
        }
    }
    
    GenotypePosteriors
    init_genotype_posteriors(const GenotypeLogMarginals& genotype_log_marginals,
                                 const GenotypeLogLikelihoods& genotype_log_likilhoods)
    {
        GenotypePosteriors result {};
        result.reserve(genotype_log_likilhoods.size());
        
        for (const auto& sample_genotype_log_likilhoods : genotype_log_likilhoods) {
            SampleGenotypePosteriors sample_result(genotype_log_marginals.size());
            
            std::transform(std::cbegin(genotype_log_marginals), std::cend(genotype_log_marginals),
                           std::cbegin(sample_genotype_log_likilhoods.second), std::begin(sample_result),
                           [] (const auto& genotype_log_marginal, const auto& genotype_log_likilhood) {
                               return genotype_log_marginal.log_probability + genotype_log_likilhood;
                           });
            
            normalise_exp(sample_result);
            
            result.emplace(sample_genotype_log_likilhoods.first, std::move(sample_result));
        }
        
        return result;
    }
    
    void update_genotype_log_marginals(GenotypeLogMarginals& marginals,
                                       const HaplotypeFrequencies& haplotype_frequencies)
    {
        for (auto& p : marginals) {
            p.log_probability = log_hardy_weinberg(p.genotype, haplotype_frequencies);
        }
    }
    
    void update_genotype_posteriors(GenotypePosteriors& genotype_posteriors,
                                    const HaplotypeFrequencies& haplotype_frequencies,
                                    const GenotypeLogMarginals& genotype_log_marginals,
                                    const GenotypeLogLikelihoods& genotype_log_likilhoods)
    {
        for (auto& sample_genotype_posteriors : genotype_posteriors) {
            
            const auto& sample_genotype_log_likilhoods = genotype_log_likilhoods.at(sample_genotype_posteriors.first);
            
            std::transform(std::cbegin(genotype_log_marginals), std::cend(genotype_log_marginals),
                           std::cbegin(sample_genotype_log_likilhoods),
                           std::begin(sample_genotype_posteriors.second),
                           [] (const auto& genotype_log_marginal, const auto& genotype_log_likilhood) {
                               return genotype_log_marginal.log_probability + genotype_log_likilhood;
                           });
            
            normalise_exp(sample_genotype_posteriors.second);
        }
    }
    
    HaplotypeFrequencies
    compute_haplotype_frequencies(const HaplotypePriorCounts& haplotype_prior_counts,
                                  const std::vector<Genotype<Haplotype>>& genotypes,
                                  const GenotypePosteriors& genotype_posteriors)
    {
        HaplotypeFrequencies result {};
        result.reserve(haplotype_prior_counts.size());
        
        const double norm = genotype_posteriors.size() * genotypes.front().ploidy();
        
        for (const auto& haplotype_count : haplotype_prior_counts) {
            const auto& haplotype = haplotype_count.first;
            
            double haplotype_frequency = haplotype_count.second;
            
            for (size_t k {}; k < genotypes.size(); ++k) {
                auto n = genotypes[k].count(haplotype);
                
                if (n > 0) {
                    double t {};
                    
                    for (const auto& sample_genotype_log_posteriors : genotype_posteriors) {
                        t += sample_genotype_log_posteriors.second[k];
                    }
                    
                    haplotype_frequency += n * t;
                }
            }
            
            result.emplace(haplotype, haplotype_frequency / norm);
        }
        
        return result;
    }
    
    double do_em_iteration(const std::vector<Haplotype>& haplotypes,
                           const std::vector<Genotype<Haplotype>>& genotypes,
                           GenotypePosteriors& genotype_posteriors,
                           HaplotypeFrequencies& haplotype_frequencies,
                           GenotypeLogMarginals& genotype_log_marginals,
                           const GenotypeLogLikelihoods& genotype_log_likilhoods,
                           const HaplotypePriorCounts& haplotype_prior_counts)
    {
        auto new_frequencies = compute_haplotype_frequencies(haplotype_prior_counts, genotypes,
                                                             genotype_posteriors);
        
        auto max_change = max_haplotype_frequency_change(haplotype_frequencies, new_frequencies);
        
        update_genotype_log_marginals(genotype_log_marginals, new_frequencies);
        update_genotype_posteriors(genotype_posteriors, new_frequencies,
                                   genotype_log_marginals, genotype_log_likilhoods);
        
        haplotype_frequencies = new_frequencies;
        
        return max_change;
    }
    
    Population::Latents
    make_latents(std::vector<Genotype<Haplotype>>&& genotypes,
                 GenotypePosteriors&& genotype_posteriors,
                 HaplotypeFrequencies&& haplotype_frequencies)
    {
        Population::GenotypeProbabilities result {};
        for (auto&& sample_genotype_log_posteriors : genotype_posteriors) {
            Population::SampleGenotypeProbabilities sample_result {};
            sample_result.reserve(genotypes.size());
            
            auto it = std::cbegin(sample_genotype_log_posteriors.second);
            
            for (const auto& genotype : genotypes) {
                sample_result.emplace(genotype, *it);
                ++it;
            }
            
            result.emplace(std::move(sample_genotype_log_posteriors.first), std::move(sample_result));
        }
        
        return Population::Latents {std::move(result), std::move(haplotype_frequencies)};
    }
        
    // private methods
        
    Population::Latents
    Population::evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads, ReferenceGenome& reference)
    {
        auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
        
        ReadModel read_model {ploidy_, max_sample_read_count(reads), haplotypes.size()};
        
        const auto genotype_log_likilhoods = compute_genotype_log_likelihoods(genotypes, reads, read_model);
        
        read_model.clear_cache();
        
        auto haplotype_prior_counts = compute_haplotype_prior_counts(haplotypes, reference, haplotype_prior_model_);
        
        auto haplotype_frequencies  = init_haplotype_frequencies(haplotype_prior_counts);
        auto genotype_log_marginals = init_genotype_log_marginals(genotypes, haplotype_frequencies);
        
        auto genotype_posteriors = init_genotype_posteriors(genotype_log_marginals, genotype_log_likilhoods);
        
        for (unsigned n {}; n < max_em_iterations_; ++n) {
            //std::cout << "EM iteration " << n << std::endl;
            double c = do_em_iteration(haplotypes, genotypes, genotype_posteriors, haplotype_frequencies,
                                       genotype_log_marginals, genotype_log_likilhoods, haplotype_prior_counts);
            if (c < em_epsilon_) break;
        }
        
        return make_latents(std::move(genotypes), std::move(genotype_posteriors), std::move(haplotype_frequencies));
    }
    
    } // namespace GenotypeModel
} // namespace Octopus

