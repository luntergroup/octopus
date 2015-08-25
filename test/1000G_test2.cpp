//
//  1000G_test2.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <cstddef>
#include <vector>
#include <unordered_map>

#include "test_common.h"
#include "test_utils.h"
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
#include "haplotype_prior_model.h"
#include "variational_bayes_genotype_model.h"
#include "maths.h"
#include "read_filter.h"
#include "read_filters.h"
#include "haplotype_tree.h"

using std::cout;
using std::endl;

//BOOST_AUTO_TEST_CASE(1000G test 2: 11:27282193-27282290)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager {std::vector<std::string> {human_1000g_bam1}};
//    
//    auto a_region = parse_region("11:27282193-27282290", human);
//    
//    auto samples = a_read_manager.get_sample_ids();
//    
//    cout << "Calling variants in region " << a_region << " for samples ";
//    for (const auto& sample : samples) cout << sample << " ";
//    cout << endl;
//    
//    auto reads = a_read_manager.fetch_reads(samples, a_region);
//    
//    CandidateVariantGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//    
//    for (const auto& sample_reads : reads) {
//        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
//    }
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
//    cout << "there are " << haplotypes.size() << " haplotypes" << endl;
//    
//    Haplotype reference_haplotype {human, a_region};
//    
//    Haplotype true_haplotype {human, a_region};
//    true_haplotype.push_back(Allele {parse_region("11:27282267-27282268", human), "G"});
//    
//    Haplotype false_haplotype1 {human, a_region};
//    false_haplotype1.push_back(Allele {parse_region("11:27282258-27282259", human), "C"});
//    false_haplotype1.push_back(Allele {parse_region("11:27282267-27282268", human), "G"});
//    
//    Haplotype false_haplotype2 {human, a_region};
//    false_haplotype2.push_back(Allele {parse_region("11:27282199-27282200", human), "C"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282207-27282208", human), "C"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282218-27282219", human), "C"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282239-27282240", human), "T"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282267-27282268", human), "G"});
//    
//    Haplotype false_haplotype3 {human, a_region};
//    false_haplotype3.push_back(Allele {parse_region("11:27282239-27282240", human), "T"});
//    
//    Haplotype false_haplotype4 {human, a_region};
//    false_haplotype4.push_back(Allele {parse_region("11:27282199-27282200", human), "C"});
//    false_haplotype4.push_back(Allele {parse_region("11:27282207-27282208", human), "C"});
//    false_haplotype4.push_back(Allele {parse_region("11:27282218-27282219", human), "C"});
//    false_haplotype4.push_back(Allele {parse_region("11:27282239-27282240", human), "T"});
//    false_haplotype4.push_back(Allele {parse_region("11:27282258-27282259", human), "C"});
//    false_haplotype4.push_back(Allele {parse_region("11:27282267-27282268", human), "G"});
//    false_haplotype4.push_back(Allele {parse_region("11:27282272-27282273", human), "G"});
//    
//    std::unordered_map<Haplotype, std::string> haplotype_labels {};
//    haplotype_labels.emplace(reference_haplotype, "ref");
//    haplotype_labels.emplace(true_haplotype, "true");
//    haplotype_labels.emplace(false_haplotype1, "false1");
//    haplotype_labels.emplace(false_haplotype2, "false2");
//    haplotype_labels.emplace(false_haplotype3, "false3");
//    haplotype_labels.emplace(false_haplotype4, "false4");
//    
////    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
////    
////    for (const auto& haplotype : haplotypes) {
////        pseudo_counts[haplotype] = 0.028;
////    }
////    
////    pseudo_counts[reference_haplotype] = 3.0;
////    pseudo_counts[true_haplotype]      = 0.37;
////    pseudo_counts[false_haplotype3]    = 0.37;
//    
//    auto haplotype_priors = get_haplotype_prior_probabilities<double>(haplotypes, candidates.cbegin(),
//                                                                      candidates.cend());
//    
//    auto pseudo_counts = get_prior_pseudo_counts(haplotype_priors, reference_haplotype, 3.0);
//    
////    cout << pseudo_counts[reference_haplotype] << endl;
////    cout << pseudo_counts[true_haplotype] << endl;
////    cout << pseudo_counts[false_haplotype1] << endl;
////    cout << pseudo_counts[false_haplotype2] << endl;
////    cout << pseudo_counts[false_haplotype3] << endl;
////    cout << pseudo_counts[false_haplotype4] << endl;
////    
////    exit(0);
//    
//    unsigned ploidy {2};
//    
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//    
//    cout << "there are " << genotypes.size() << " genotypes" << endl;
//    
//    Genotype ref_ref {reference_haplotype, reference_haplotype};
//    Genotype ref_true {reference_haplotype, true_haplotype};
//    Genotype ref_false1 {reference_haplotype, false_haplotype1};
//    Genotype ref_false2 {reference_haplotype, false_haplotype2};
//    Genotype ref_false3 {reference_haplotype, false_haplotype3};
//    Genotype true_true {true_haplotype, true_haplotype};
//    Genotype true_false1 {true_haplotype, false_haplotype2};
//    Genotype true_false2 {true_haplotype, false_haplotype2};
//    Genotype true_false3 {true_haplotype, false_haplotype3};
//    Genotype true_false4 {true_haplotype, false_haplotype4};
//    Genotype false1_false1 {false_haplotype1, false_haplotype1};
//    Genotype false1_false2 {false_haplotype1, false_haplotype2};
//    Genotype false1_false3 {false_haplotype1, false_haplotype3};
//    Genotype false2_false2 {false_haplotype2, false_haplotype2};
//    Genotype false2_false3 {false_haplotype2, false_haplotype3};
//    Genotype false3_false3 {false_haplotype3, false_haplotype3};
//    
//    ReadModel a_read_model {ploidy};
//    
//    cout << endl;
//    for (const auto& sample_reads : reads) {
//        cout << "genotype likelihoods for sample " << sample_reads.first << endl;
//        cout << "ref,ref     = " << a_read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), ref_ref, sample_reads.first) << endl;
//        cout << "ref,true    = " << a_read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), ref_true, sample_reads.first) << endl;
//        cout << "ref,false2  = " << a_read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), ref_false2, sample_reads.first) << endl;
//        cout << "ref,false3  = " << a_read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), ref_false3, sample_reads.first) << endl;
//        cout << "true,false2 = " << a_read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), true_false2, sample_reads.first) << endl;
//        cout << "true,false3 = " << a_read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), true_false3, sample_reads.first) << endl;
//        cout << "true,false4 = " << a_read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), true_false4, sample_reads.first) << endl;
//        cout << endl;
//    }
//    
//    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
//    
//    SamplesReads the_reads {};
//    for (const auto& sample_reads : reads) {
//        the_reads.emplace(sample_reads.first, std::make_pair(sample_reads.second.cbegin(),
//                                                             sample_reads.second.cend()));
//    }
//    
//    auto results = update_parameters(the_model, genotypes, pseudo_counts, the_reads, 3);
//    auto responsabilities        = results.first;
//    auto posterior_pseudo_counts = results.second;
//    
//    cout << "haplotype posterior counts:" << endl;
//    cout << "ref    = " << posterior_pseudo_counts[reference_haplotype] << endl;
//    cout << "true   = " << posterior_pseudo_counts[true_haplotype] << endl;
//    cout << "false1 = " << posterior_pseudo_counts[false_haplotype1] << endl;
//    cout << "false2 = " << posterior_pseudo_counts[false_haplotype2] << endl;
//    cout << "false3 = " << posterior_pseudo_counts[false_haplotype3] << endl;
//    cout << "false4 = " << posterior_pseudo_counts[false_haplotype4] << endl;
//    cout << endl;
//    
//    for (const auto& sample : samples) {
//        cout << "**** Sample " << sample << " ****" << endl << endl;
//        
//        auto& sample_responsabilities = responsabilities[sample];
//        
//        std::sort(genotypes.begin(), genotypes.end(), [&sample_responsabilities] (const auto& g1, const auto& g2) {
//            return sample_responsabilities[g1] > sample_responsabilities[g2];
//        });
//        
//        cout << "top genotype posteriors:" << endl;
//        for (unsigned i {}; i < 5; ++i) {
//            cout << get_label(haplotype_labels, genotypes.at(i).at(0)) << ","
//            << get_label(haplotype_labels, genotypes.at(i).at(1)) << " "
//            << sample_responsabilities.at(genotypes.at(i)) << endl;
//        }
//        cout << endl;
//        
//        cout << "allele posteriors:" << endl;
//        for (const auto& variant : candidates) {
//            cout << variant << " "
//            << the_model.posterior_probability_allele_in_sample(variant.get_reference_allele(), haplotypes,
//                                                                sample_responsabilities, genotypes)
//            << " "
//            << the_model.posterior_probability_allele_in_sample(variant.get_alternative_allele(), haplotypes,
//                                                                sample_responsabilities, genotypes)
//            << endl;
//        }
//        cout << endl << endl;
//    }
//}
