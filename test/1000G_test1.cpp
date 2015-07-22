//
//  1000G_test1.cpp
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
#include "reference_genome_factory.h"
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

//BOOST_AUTO_TEST_CASE(1000G test 1: 11:67503118-67503253)
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
//    cout << "Calling variants in region " << a_region << " for samples ";
//    for (const auto& sample : samples) cout << sample << " ";
//    cout << endl;
//    
//    auto reads = a_read_manager.fetch_reads(samples, a_region);
//    
//    CandidateVariantGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 10));
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
//    true_haplotype.push_back(Allele {parse_region("11:67503147-67503148", human), "A"});
//    true_haplotype.push_back(Allele {parse_region("11:67503214-67503215", human), "A"});
//    
//    std::unordered_map<Haplotype, std::string> haplotype_labels {};
//    haplotype_labels.emplace(reference_haplotype, "ref");
//    haplotype_labels.emplace(true_haplotype, "true");
//    
//    unsigned ploidy {2};
//    
//    ReadModel a_read_model {ploidy};
//    
//    Genotype hom_ref {reference_haplotype, reference_haplotype};
//    Genotype het_alt {reference_haplotype, true_haplotype};
//    Genotype hom_alt {true_haplotype, true_haplotype};
//    
//    for (const auto& sample_reads : reads) {
//        cout << "genotype likelihoods for sample " << sample_reads.first << endl;
//        cout << "ref,ref  = " << a_read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), hom_ref, sample_reads.first) << endl;
//        cout << "ref,alt  = " << a_read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), het_alt, sample_reads.first) << endl;
//        cout << "alt,alt  = " << a_read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), hom_alt, sample_reads.first) << endl;
//        cout << endl;
//    }
//    
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//    //std::vector<Genotype> genotypes {hom_ref, het_alt, hom_alt};
//    
//    cout << "there are " << genotypes.size() << " genotypes" << endl;
//    
//    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
//    
////    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
////    for (const auto& haplotype : haplotypes) {
////        pseudo_counts[haplotype] = 0.196;
////    }
////    
////    pseudo_counts[reference_haplotype] = 134;
//    
//    auto haplotype_priors = get_haplotype_prior_probabilities<double>(haplotypes, candidates.cbegin(),
//                                                                      candidates.cend());
//    
//    auto pseudo_counts = get_prior_pseudo_counts(haplotype_priors, reference_haplotype, 3);
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
//    cout << "ref  = "  << posterior_pseudo_counts[reference_haplotype] << endl;
//    cout << "true = "  << posterior_pseudo_counts[true_haplotype] << endl;
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
//        for (unsigned i {}; i < 3; ++i) {
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
