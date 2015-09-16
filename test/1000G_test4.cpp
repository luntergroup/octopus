//
//  1000G_test4.cpp
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
#include "maths.h"
#include "read_filter.h"
#include "read_filters.h"
#include "haplotype_tree.h"

using std::cout;
using std::endl;

//BOOST_AUTO_TEST_CASE(1000G test 4: 16:9299900-9300110)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human(a_factory.make(human_reference_fasta));
//
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
//
//    CandidateVariantGenerator candidate_generator {};
//
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//
//    auto sample_ids = a_read_manager.get_sample_ids();
//    auto the_sample_id = sample_ids.at(0);
//
//    auto a_region = parse_region("16:9300000-9300100", human);
//
//    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
//
//    using ReadIterator = std::vector<AlignedRead>::const_iterator;
//    ReadFilter<ReadIterator> a_read_filter {};
//    a_read_filter.register_filter(is_mapped);
//
//    std::vector<AlignedRead> good_reads {}, bad_reads {};
//    good_reads.reserve(reads.size());
//    bad_reads.reserve(reads.size());
//    a_read_filter.filter_reads(std::make_move_iterator(reads.begin()), std::make_move_iterator(reads.end()),
//                               ContextBackInserter(good_reads), ContextBackInserter(bad_reads));
//
//    candidate_generator.add_reads(good_reads.cbegin(), good_reads.cend());
//
//    auto candidates = candidate_generator.get_candidates(a_region);
//
//    std::vector<Variant> aligned_candidates {};
//    std::transform(candidates.begin(), candidates.end(),
//                   std::back_inserter(aligned_candidates),
//                   [&human] (const auto& a_variant) {
//                       return left_align(a_variant, human);
//                   });
//
//    Haplotype reference_haplotype {human, a_region};
//
//    //for (auto& c : candidates) { cout << c << endl; }
//
//    // 2 snps
//    Haplotype hap1 {human, a_region};
//    add_to_back(candidates[3], hap1);
//    add_to_back(candidates[5], hap1);
//
//    // smaller deletion
//    Haplotype hap2 {human, a_region};
//    add_to_back(candidates[2], hap2);
//
//    // larger deletion
//    Haplotype hap3 {human, a_region};
//    add_to_back(candidates[1], hap3);
//
//    // smaller insertion
//    Haplotype hap4 {human, a_region};
//    add_to_back(candidates[0], hap4);
//
//    // smaller insertion & smaller deletion
//    Haplotype hap5 {human, a_region};
//    add_to_back(candidates[0], hap5);
//    add_to_back(candidates[2], hap5);
//
//    // bigger insertion
//    Haplotype hap6 {human, a_region};
//    add_to_back(candidates[6], hap6);
//
//    // snps & smaller insertion
//    Haplotype hap7 {human, a_region};
//    add_to_back(candidates[0], hap7);
//    add_to_back(candidates[3], hap7);
//    add_to_back(candidates[5], hap7);
//
//    // snps & bigger insertion
//    Haplotype hap8 {human, a_region};
//    add_to_back(candidates[3], hap8);
//    add_to_back(candidates[5], hap8);
//    add_to_back(candidates[6], hap8);
//
//    // smaller insertion & bigger insertion
//    Haplotype hap9 {human, a_region};
//    add_to_back(candidates[0], hap9);
//    add_to_back(candidates[6], hap9);
//
//    std::vector<Haplotype> haplotypes {reference_haplotype, hap1, hap2, hap3, hap4, hap5, hap6, hap7, hap8, hap9};
//
//    unsigned ploidy {2};
//    ReadModel a_read_model {ploidy};
//    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
//    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
//    for (const auto& h : haplotypes) { pseudo_counts[h] = 1; }
//
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//
//    SamplesReads the_reads {};
//    the_reads.push_back({good_reads.cbegin(), good_reads.cend()});
//
//    auto results = update_parameters(the_model, genotypes, pseudo_counts, the_reads, 10);
//    auto responsabilities        = results.first;
//    auto posterior_pseudo_counts = results.second;
//
////    for (const auto& h : haplotypes) { cout << h << endl; }
////
////    for (auto& p : responsabilities[0]) {
////        if (p.second > 0.2) {
////            cout << p.first << " " << p.second << endl;
////        }
////    }
//}

