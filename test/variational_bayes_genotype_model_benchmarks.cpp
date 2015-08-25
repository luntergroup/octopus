//
//  variational_bayes_genotype_model_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <random>

#include "test_common.h"
#include "test_utils.h"
#include "benchmark_utils.h"
#include "reference_genome.h"
#include "read_manager.h"
#include "allele.h"
#include "variant.h"
#include "variant_utils.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"
#include "bayesian_genotype_model.h"
#include "variational_bayes_genotype_model.h"
#include "maths.h"
#include "read_filter.h"
#include "read_filters.h"
#include "haplotype_tree.h"

using Octopus::HaplotypeTree;
using Octopus::ReadModel;
using Octopus::BayesianGenotypeModel::VariationalBayesGenotypeModel;

//BOOST_AUTO_TEST_CASE(responsability calculation)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    auto a_region = parse_region("1:10000000-10001000", human);
//    
//    ReadManager a_read_manager {std::vector<std::string> {human_1000g_bam1}};
//    
//    auto sample_id = a_read_manager.get_sample_ids().at(0);
//    
//    auto reads = a_read_manager.fetch_reads(sample_id, a_region);
//    
//    unsigned ploidy {2};
//    ReadModel read_model {ploidy};
//    VariationalBayesGenotypeModel the_model {read_model, ploidy};
//    
//    Haplotype reference_haplotype {human, a_region};
//    
//    std::vector<Haplotype> haplotypes {reference_haplotype};
//    
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//    
//    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
//    for (const auto& haplotype : haplotypes) {
//        pseudo_counts[haplotype] = 1.0;
//    }
//    
//    auto f_genotype_resp = [&sample_id, &the_model, &reads, &genotypes, &pseudo_counts] () {
//        the_model.genotype_responsabilities(genotypes, reads.cbegin(), reads.cend(), pseudo_counts, sample_id);
//    };
//    
//    auto time = benchmark<std::chrono::microseconds>(f_genotype_resp, 1).count();
//    
//    std::cout << "genotype responsabilities time " << time << std::endl;
//}

//BOOST_AUTO_TEST_CASE(where is the time being taken?)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager {std::vector<std::string> {human_1000g_bam1}};
//    
//    auto a_region = parse_region("11:67503118-67503253", human);
//    
//    auto samples = a_read_manager.get_sample_ids();
//    
//    auto reads = a_read_manager.fetch_reads(samples, a_region);
//    
//    CandidateVariantGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//    
//    candidate_generator.add_reads(reads["HG00101"].cbegin(), reads["HG00101"].cend());
//    
//    auto candidates = candidate_generator.get_candidates(a_region);
//    
//    HaplotypeTree haplotype_tree {human};
//    
//    for (const auto& candidate : candidates) {
//        haplotype_tree.extend(candidate.get_reference_allele());
//        haplotype_tree.extend(candidate.get_alternative_allele());
//    }
//    
//    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
//    
//    unique_least_complex(haplotypes);
//    
//    unsigned ploidy {2};
//    
//    ReadModel a_read_model {ploidy};
//    
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//    
//    auto f_read_probs = [&samples, &genotypes, &reads, &a_read_model] () {
//        for (const auto& sample : samples) {
//            for (const auto& genotype : genotypes) {
//                a_read_model.log_probability(reads[sample].cbegin(), reads[sample].cend(), genotype, sample);
//            }
//        }
//        a_read_model.clear_cache();
//    };
//    
//    ReadModel a_new_read_model {ploidy};
//    
//    VariationalBayesGenotypeModel the_model {a_new_read_model, ploidy};
//    
//    Octopus::BayesianGenotypeModel::HaplotypePseudoCounts<double> pseudo_counts {};
//    for (const auto& haplotype : haplotypes) {
//        pseudo_counts[haplotype] = 1.0;
//    }
//    
//    Octopus::BayesianGenotypeModel::GenotypeProbabilities<std::string, double> responsabilities(samples.size());
//    
//    auto f_genotype_resp = [&the_model, &reads, &genotypes, &pseudo_counts, &responsabilities] () {
//        for (const auto& sample : reads) {
//            responsabilities[sample.first] = the_model.genotype_responsabilities(genotypes,
//                                                                                 sample.second.cbegin(),
//                                                                                 sample.second.cend(),
//                                                                                 pseudo_counts,
//                                                                                 sample.first);
//        }
//    };
//    
//    ReadModel another_new_read_model {ploidy};
//    
//    VariationalBayesGenotypeModel a_new_model {another_new_read_model, ploidy};
//    
//    SamplesReads the_reads {};
//    for (const auto& sample_reads : reads) {
//        the_reads.emplace(sample_reads.first, std::make_pair(sample_reads.second.cbegin(),
//                                                             sample_reads.second.cend()));
//    }
//    
//    auto f_update = [&a_new_model, &the_reads, &genotypes, &pseudo_counts] () {
//        update_parameters(a_new_model, genotypes, pseudo_counts, the_reads, 3);
//    };
//    
//    auto read_model_time     = benchmark<std::chrono::milliseconds>(f_read_probs, 1).count();
//    auto genotype_model_time = benchmark<std::chrono::milliseconds>(f_genotype_resp, 1).count();
//    auto model_update_time   = benchmark<std::chrono::milliseconds>(f_update, 1).count();
//    
//    std::cout << "read_model_time " << read_model_time << std::endl;
//    std::cout << "genotype_model_time " << genotype_model_time << std::endl;
//    std::cout << "model_update_time " << model_update_time << std::endl;
//}

//BOOST_AUTO_TEST_CASE(prior calculation)
//{
//    std::size_t size {100};
//    std::vector<double> probs(size);
//    std::vector<double> alphas(size);
//    
//    std::default_random_engine generator;
//    std::uniform_real_distribution<double> distribution(0.0, 1.0);
//    
//    std::generate_n(probs.begin(), size, [&distribution, &generator] () { return distribution(generator); });
//    
//    auto f_prior = [&probs, &alphas] () {
//        for (std::size_t i {}; i < probs.size(); ++i) {
//            alphas[i] = std::log(digamma_inv(probs[i]));
//        }
//    };
//    
//    auto f_prior_time = benchmark<std::chrono::microseconds>(f_prior, 10000).count();
//    
//    std::cout << "f_prior_time " << f_prior_time << std::endl;
//}
