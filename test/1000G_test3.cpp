//
//  1000G_test3.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

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
#include "variant_candidate_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"
#include "variational_bayes_genotype_model.h"
#include "maths.h"
#include "read_filter.h"
#include "read_filters.h"
#include "haplotype_tree.h"

using std::cout;
using std::endl;

//TEST_CASE("1000G test 3: 2:104142870-104142884", "[validation]")
//{
//    unsigned ploidy {2};
//
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human(a_factory.make(human_reference_fasta));
//
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1, human_1000g_bam2});
//
//    VariantCandidateGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//
//    auto a_region = parse_region("2:104142870-104142884", human);
//
//    auto reference_sequence = human.get_sequence(a_region);
//
//    auto sample_ids = a_read_manager.get_sample_ids();
//
//    auto some_reads = a_read_manager.fetch_reads(sample_ids, a_region);
//
//    candidate_generator.add_reads(some_reads[sample_ids[0]].cbegin(), some_reads[sample_ids[0]].cend());
//    candidate_generator.add_reads(some_reads[sample_ids[1]].cbegin(), some_reads[sample_ids[1]].cend());
//
//    auto variants = candidate_generator.get_candidates(a_region);
//
//    REQUIRE(variants.size() == 3);
//
//    Haplotype reference_haplotype {human, a_region}; // there are no reads completely supporting the reference
//
//    Haplotype hap1 {human, a_region};
//    add_to_back(variants[0], hap1); // well supported insert
//    add_to_back(variants[2], hap1); // well supported snp
//
//    Haplotype hap2 {human, a_region};
//    add_to_back(variants[1], hap2); // this is a low quality snp
//
//    Haplotype hap3 {human, a_region};
//    add_to_back(variants[0], hap3);
//    add_to_back(variants[1], hap3);
//    add_to_back(variants[2], hap3);
//
//    unsigned num_haplotypes {4};
//    std::vector<Haplotype> haplotypes {reference_haplotype, hap1, hap2, hap3};
//
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//
//    REQUIRE(genotypes.size() == num_genotypes(num_haplotypes, ploidy));
//
//    ReadModel a_read_model {ploidy};
//
//    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
//
//    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
//    pseudo_counts[reference_haplotype] = 1;
//    pseudo_counts[hap1]                = 1;
//    pseudo_counts[hap2]                = 1;
//    pseudo_counts[hap3]                = 1;
//
//    SamplesReads the_reads {};
//    the_reads.push_back({some_reads[sample_ids[0]].cbegin(), some_reads[sample_ids[0]].cend()});
//    the_reads.push_back({some_reads[sample_ids[1]].cbegin(), some_reads[sample_ids[1]].cend()});
//
//    auto results = update_parameters(the_model, genotypes, pseudo_counts, the_reads, 10);
//    auto responsabilities        = results.first;
//    auto posterior_pseudo_counts = results.second;
//
////    for (auto& p : responsabilities[0]) {
////        cout << p.first << " " << p.second << endl;
////    }
////    cout << endl;
////    for (auto& p : responsabilities[1]) {
////        cout << p.first << " " << p.second << endl;
////    }
////    cout << endl;
////    for (const auto& c : posterior_pseudo_counts) {
////        cout << c.first << " " << c.second << endl;
////    }
//}
