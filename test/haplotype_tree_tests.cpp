////
////  haplotype_tree_tests.cpp
////  Octopus
////
////  Created by Daniel Cooke on 22/03/2015.
////  Copyright (c) 2015 Oxford University. All rights reserved.
////
//
//#define BOOST_TEST_DYN_LINK
//
//#include <boost/test/unit_test.hpp>
//
//#include <iostream>
//#include <vector>
//#include <string>
//#include <cstddef>
//#include <algorithm>
//
//#include "test_common.hpp"
//#include "reference_genome.hpp"
//#include "read_manager.hpp"
//#include "variant.hpp"
//#include "candidate_variant_generator.hpp"
//#include "candidate_generator_builder.hpp"
//#include "haplotype.hpp"
//#include "haplotype_tree.hpp"
//#include "mappable_algorithms.hpp"
//
//using std::cout;
//using std::endl;
//
//using Octopus::HaplotypeTree;
//using Octopus::CandidateGeneratorBuilder;
//
//BOOST_AUTO_TEST_SUITE(Components)
//BOOST_AUTO_TEST_SUITE(HaplotypeTreeTests)
//
//static void sort(std::vector<Haplotype>& haplotypes)
//{
//    std::sort(std::begin(haplotypes), std::end(haplotypes));
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_tree_splits_overlapping_snps_into_different_branches)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
//    const Allele allele2 {parse_region("4:1000000-1000001", human), "C"};
//    const Allele allele3 {parse_region("4:1000000-1000001", human), "G"};
//    const Allele allele4 {parse_region("4:1000001-1000002", human), "G"};
//    const Allele allele5 {parse_region("4:1000001-1000002", human), "C"};
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    
//    haplotype_tree.extend(allele1);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 1);
//    
//    haplotype_tree.extend(allele2);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
//    
//    haplotype_tree.extend(allele3);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
//    
//    haplotype_tree.extend(allele4);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
//    
//    haplotype_tree.extend(allele5);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 6);
//}
//
//BOOST_AUTO_TEST_CASE(clear_leaves_the_tree_empty)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
//    const Allele allele2 {parse_region("4:1000000-1000001", human), "C"};
//    const Allele allele3 {parse_region("4:1000000-1000001", human), "G"};
//    const Allele allele4 {parse_region("4:1000001-1000002", human), "G"};
//    const Allele allele5 {parse_region("4:1000001-1000002", human), "C"};
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    haplotype_tree.extend(allele1).extend(allele2).extend(allele3).extend(allele4).extend(allele5);
//    
//    BOOST_REQUIRE(haplotype_tree.num_haplotypes() == 6);
//    
//    haplotype_tree.clear();
//    
//    BOOST_CHECK(haplotype_tree.empty());
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 0);
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_tree_ignores_duplicate_alleles_coming_from_same_allele)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
//    const Allele allele2 {parse_region("4:1000000-1000001", human), "C"};
//    const Allele allele3 {parse_region("4:1000000-1000001", human), "A"};
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    haplotype_tree.extend(allele1).extend(allele2).extend(allele3);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
//    
//    const Allele allele4 {parse_region("4:1000001-1000001", human), "A"};
//    const Allele allele5 {parse_region("4:1000001-1000001", human), "C"};
//    const Allele allele6 {parse_region("4:1000001-1000001", human), "C"};
//    
//    haplotype_tree.extend(allele4).extend(allele5).extend(allele6);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 4);
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_tree_ignores_insertions_followed_immediatly_by_deletions_and_vice_versa)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("16:9300037-9300037", human), "TG"};
//    const Allele allele2 {parse_region("16:9300037-9300051", human), ""};
//    
//    HaplotypeTree haplotype_tree {"16", human};
//    haplotype_tree.extend(allele1).extend(allele2);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 1);
//    
//    const auto region = parse_region("16:9300037-9300037", human);
//    
//    const auto haplotypes = haplotype_tree.extract_haplotypes(region);
//    
//    BOOST_CHECK(haplotypes[0].contains(allele1));
//    BOOST_CHECK(!haplotypes[0].contains(allele2));
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_tree_does_not_bifurcate_on_alleles_positioned_past_the_leading_alleles)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
//    const Allele allele2 {parse_region("4:1000001-1000002", human), "C"};
//    const Allele allele3 {parse_region("4:1000002-1000002", human), "GC"};
//    const Allele allele4 {parse_region("4:1000005-1000007", human), ""};
//    const Allele allele5 {parse_region("4:1000007-1000008", human), "G"};
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    haplotype_tree.extend(allele1).extend(allele2).extend(allele3).extend(allele4).extend(allele5);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 1);
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_tree_can_generate_haplotypes_in_a_region)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
//    const Allele allele2 {parse_region("4:1000002-1000003", human), "C"};
//    const Allele allele3 {parse_region("4:1000002-1000003", human), "G"};
//    const Allele allele4 {parse_region("4:1000004-1000005", human), "T"};
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    haplotype_tree.extend(allele1).extend(allele2).extend(allele3).extend(allele4);
//    
//    auto region = parse_region("4:1000000-1000005", human);
//    
//    auto haplotypes = haplotype_tree.extract_haplotypes(region);
//    
//    BOOST_CHECK(haplotypes.size() == 2);
//    
//    sort(haplotypes);
//    
//    BOOST_CHECK(haplotypes[0].get_sequence() == "ATCCT");
//    BOOST_CHECK(haplotypes[1].get_sequence() == "ATGCT");
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_tree_can_generate_haplotypes_ending_in_different_regions)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
//    const Allele allele2 {parse_region("4:1000002-1000006", human), ""};
//    const Allele allele3 {parse_region("4:1000002-1000003", human), "G"};
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    
//    haplotype_tree.extend(allele1).extend(allele2).extend(allele3);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
//    
//    const auto region1 = parse_region("4:1000000-1000006", human);
//    
//    auto haplotypes = haplotype_tree.extract_haplotypes(region1);
//    
//    BOOST_CHECK(haplotypes.size() == 2);
//    
//    const auto region2 = parse_region("4:1000000-1000003", human);
//    
//    haplotypes = haplotype_tree.extract_haplotypes(region2);
//    
//    BOOST_CHECK(haplotypes.size() == 2);
//    
//    sort(haplotypes);
//    
//    BOOST_CHECK(haplotypes[0].get_sequence() == "ATG");
//    BOOST_CHECK(haplotypes[1].get_sequence() == "AT");
//}
//
//BOOST_AUTO_TEST_CASE(leading_haplotypes_can_be_removed_from_the_tree)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
//    const Allele allele2 {parse_region("4:1000002-1000003", human), "C"};
//    const Allele allele3 {parse_region("4:1000002-1000003", human), "G"};
//    const Allele allele4 {parse_region("4:1000004-1000005", human), "T"};
//    const Allele allele5 {parse_region("4:1000004-1000005", human), "C"};
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    
//    haplotype_tree.extend(allele1).extend(allele2).extend(allele3).extend(allele4).extend(allele5);
//    
//    const auto region = encompassing_region(allele1, allele5);
//    
//    auto haplotypes = haplotype_tree.extract_haplotypes(region);
//    
//    BOOST_CHECK(haplotypes.size() == 4);
//    
//    sort(haplotypes);
//    
//    haplotype_tree.prune_all(haplotypes[0]); // ATCCC
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
//    
//    haplotype_tree.prune_all(haplotypes[1]); // ATCCT
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
//    
//    haplotypes = haplotype_tree.extract_haplotypes(region);
//    
//    sort(haplotypes);
//    
//    BOOST_CHECK(haplotypes[0].get_sequence() == "ATGCC");
//    BOOST_CHECK(haplotypes[1].get_sequence() == "ATGCT");
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_tree_only_contains_haplotypes_with_added_alleles)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "C"};
//    const Allele allele2 {parse_region("4:1000001-1000002", human), "T"};
//    const Allele allele3 {parse_region("4:1000001-1000002", human), "G"};
//    const Allele allele4 {parse_region("4:1000002-1000003", human), "C"};
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    haplotype_tree.extend(allele1).extend(allele2).extend(allele3).extend(allele4);
//    
//    const auto region = encompassing_region(allele1, allele4); // reference = CTC
//    
//    Haplotype hap1 {region, human};
//    hap1.push_back(allele1);
//    hap1.push_back(allele2);
//    hap1.push_back(allele4);
//    
//    BOOST_REQUIRE(hap1.get_sequence() == "CTC");
//    BOOST_CHECK(haplotype_tree.contains(hap1));
//    
//    Haplotype hap2 {region, human};
//    hap2.push_back(allele1);
//    hap2.push_back(allele3);
//    hap2.push_back(allele4);
//    
//    BOOST_REQUIRE(hap2.get_sequence() == "CGC");
//    BOOST_CHECK(haplotype_tree.contains(hap2));
//    
//    const Allele allele5 {parse_region("4:1000000-1000001", human), "G"};
//    
//    Haplotype hap3 {region, human};
//    hap3.push_back(allele5);
//    hap3.push_back(allele2);
//    hap3.push_back(allele4);
//    
//    BOOST_REQUIRE(hap3.get_sequence() == "GTC");
//    BOOST_CHECK(!haplotype_tree.contains(hap3));
//    
//    Haplotype hap4 {region, human};
//    hap4.push_back(allele5);
//    hap4.push_back(allele3);
//    hap4.push_back(allele4);
//    
//    BOOST_REQUIRE(hap4.get_sequence() == "GGC");
//    BOOST_CHECK(!haplotype_tree.contains(hap4));
//    
//    const Allele allele6 {parse_region("4:1000001-1000002", human), "C"};
//    
//    Haplotype hap5 {region, human};
//    hap5.push_back(allele1);
//    hap5.push_back(allele6);
//    hap5.push_back(allele4);
//    
//    BOOST_REQUIRE(hap5.get_sequence() == "CCC");
//    BOOST_CHECK(!haplotype_tree.contains(hap5));
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_tree_contains_haplotypes_with_implicit_reference_alleles)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "C"};
//    const Allele allele2 {parse_region("4:1000001-1000002", human), "T"};
//    const Allele allele3 {parse_region("4:1000001-1000002", human), "G"};
//    const Allele allele4 {parse_region("4:1000002-1000003", human), "C"};
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    haplotype_tree.extend(allele1).extend(allele2).extend(allele3).extend(allele4);
//    
//    const auto region = encompassing_region(allele1, allele4); // reference = CTC
//    
//    Haplotype hap1 {region, human};
//    
//    BOOST_REQUIRE(hap1.get_sequence() == "CTC");
//    BOOST_CHECK(haplotype_tree.contains(hap1));
//    
//    Haplotype hap2 {region, human};
//    hap2.push_back(allele2);
//    
//    BOOST_REQUIRE(hap2.get_sequence() == "CTC");
//    BOOST_CHECK(haplotype_tree.contains(hap2));
//    
//    Haplotype hap3 {region, human};
//    hap3.push_back(allele3);
//    
//    BOOST_REQUIRE(hap3.get_sequence() == "CGC");
//    BOOST_CHECK(haplotype_tree.contains(hap3));
//    
//    const Allele allele5 {parse_region("4:1000000-1000001", human), "G"};
//    
//    Haplotype hap4 {region, human};
//    hap4.push_back(allele5);
//    
//    BOOST_REQUIRE(hap4.get_sequence() == "GTC");
//    BOOST_CHECK(!haplotype_tree.contains(hap4));
//    
//    const Allele allele6 {parse_region("4:1000001-1000002", human), "C"};
//    
//    Haplotype hap5 {region, human};
//    hap5.push_back(allele6);
//    
//    BOOST_REQUIRE(hap5.get_sequence() == "CCC");
//    BOOST_CHECK(!haplotype_tree.contains(hap5));
//}
//
//BOOST_AUTO_TEST_CASE(prune_all_gets_haplotypes_with_implicit_reference_alleles)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "C"};
//    const Allele allele2 {parse_region("4:1000001-1000002", human), "T"};
//    const Allele allele3 {parse_region("4:1000001-1000002", human), "G"};
//    const Allele allele4 {parse_region("4:1000002-1000003", human), "C"};
//    
//    const auto region = encompassing_region(allele1, allele4);
//    
//    Haplotype hap {region, human};
//    hap.push_back(allele2);
//    
//    BOOST_REQUIRE(hap.get_sequence() == "CTC");
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    
//    haplotype_tree.extend(allele1).extend(allele2).extend(allele3).extend(allele4);
//    
//    haplotype_tree.prune_all(hap);
//    
//    BOOST_REQUIRE(haplotype_tree.num_haplotypes() == 1);
//    BOOST_CHECK(haplotype_tree.extract_haplotypes().front().get_sequence() == "CGC");
//}
//
//BOOST_AUTO_TEST_CASE(pruned_branches_can_still_be_extended)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "C"};
//    const Allele allele2 {parse_region("4:1000001-1000002", human), "T"};
//    const Allele allele3 {parse_region("4:1000001-1000002", human), "G"};
//    const Allele allele4 {parse_region("4:1000002-1000003", human), "C"};
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    haplotype_tree.extend(allele1).extend(allele2).extend(allele3).extend(allele4);
//    
//    const auto region = encompassing_region(allele1, allele4);
//    
//    Haplotype hap {region, human};
//    hap.push_back(allele2);
//    
//    haplotype_tree.prune_all(hap);
//    
//    BOOST_REQUIRE(haplotype_tree.num_haplotypes() == 1);
//    
//    const Allele allele5 {parse_region("4:1000002-1000003", human), "T"};
//    
//    haplotype_tree.extend(allele5);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
//}
//
//BOOST_AUTO_TEST_CASE(extending_on_mnps_results_in_backtracked_bifurification)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("16:9300039-9300051", human), "TGTGTGTGCGTT"};
//    const Allele allele2 {parse_region("16:9300039-9300051", human), ""};
//    
//    const Allele allele3 {parse_region("16:9300047-9300048", human), "C"};
//    const Allele allele4 {parse_region("16:9300047-9300048", human), "T"};
//    const Allele allele5 {parse_region("16:9300050-9300051", human), "T"};
//    const Allele allele6 {parse_region("16:9300050-9300051", human), "G"};
//    
//    HaplotypeTree haplotype_tree {"16", human};
//    
//    haplotype_tree.extend(allele1).extend(allele2);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
//    
//    haplotype_tree.extend(allele3);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
//    
//    haplotype_tree.extend(allele4);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 4);
//    
//    haplotype_tree.extend(allele5);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 5);
//    
//    haplotype_tree.extend(allele6);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 8);
//    
//    auto a_region = parse_region("16:9300039-9300051", human);
//    
//    auto haplotypes = haplotype_tree.extract_haplotypes(a_region);
//    
//    BOOST_CHECK(haplotypes.size() == 8);
//    
//    sort(haplotypes);
//    
//    auto it = std::unique(haplotypes.begin(), haplotypes.end());
//    
//    haplotypes.erase(it, haplotypes.end());
//    
//    BOOST_CHECK(haplotypes.size() == 5);
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_tree_can_selectively_extend_branches)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
//    const Allele allele2 {parse_region("4:1000000-1000003", human), ""};
//    
//    haplotype_tree.extend(allele1).extend(allele2);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
//    
//    haplotype_tree.extend(allele1).extend(allele2);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
//    
//    const Allele allele3 {parse_region("4:1000001-1000002", human), "C"};
//    const Allele allele4 {parse_region("4:1000002-1000003", human), "G"};
//    
//    haplotype_tree.extend(allele3);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
//    
//    haplotype_tree.extend(allele4);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 4);
//    
//    const Allele allele5 {parse_region("4:1000003-1000004", human), "T"};
//    
//    haplotype_tree.extend(allele5);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 4);
//    
//    const Allele allele6 {parse_region("4:1000003-1000004", human), "A"};
//    
//    haplotype_tree.extend(allele6);
//    
//    BOOST_CHECK(haplotype_tree.num_haplotypes() == 8);
//}
//
//BOOST_AUTO_TEST_CASE(prune_unqiue_leaves_a_single_haplotype_which_contains_the_same_alleles_as_the_given_haplotype)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    ReadManager read_manager {NA12878_low_coverage};
//    
//    const auto region = parse_region("16:9200000-9202500", human);
//    
//    const auto sample = read_manager.get_samples().front();
//    
//    auto reads = read_manager.fetch_reads(sample, region);
//    
//    auto candidate_generator = CandidateGeneratorBuilder().set_reference(human)
//    .add_generator(CandidateGeneratorBuilder::Generator::Alignment)
//    .build();
//    
//    add_reads(reads, candidate_generator);
//    
//    auto candidates = candidate_generator.generate_candidates(region);
//    
//    BOOST_REQUIRE(candidates.size() < 20);
//    
//    HaplotypeTree haplotype_tree {"16", human};
//    
//    extend_tree(candidates, haplotype_tree);
//    
//    auto haplotypes = haplotype_tree.extract_haplotypes(region);
//    
//    sort(haplotypes);
//    
//    const auto duplicate_itr = std::adjacent_find(std::cbegin(haplotypes), std::cend(haplotypes));
//    
//    BOOST_REQUIRE(duplicate_itr != std::cend(haplotypes));  // otherwise the test is pointless
//    
//    const auto er = std::equal_range(std::begin(haplotypes), std::end(haplotypes), *duplicate_itr);
//    
//    const auto haplotype_to_prune= *er.first;
//    
//    haplotype_tree.prune_unique(haplotype_to_prune);
//    
//    auto new_haplotypes = haplotype_tree.extract_haplotypes(region);
//    
//    sort(new_haplotypes);
//    
//    const auto er2 = std::equal_range(new_haplotypes.begin(), new_haplotypes.end(), *duplicate_itr);
//    
//    BOOST_CHECK(std::distance(er2.first, er2.second) == 1);
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_tree_survives_serious_pruning)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    ReadManager read_manager {NA12878_low_coverage};
//    
//    CandidateGeneratorBuilder builder {};
//    builder.set_reference(human);
//    builder.set_max_variant_size(10);
//    builder.add_generator(CandidateGeneratorBuilder::Generator::Alignment);
//    auto candidate_generator = builder.build();
//    
//    const auto sample = read_manager.get_samples().front();
//    
//    const auto region = parse_region("16:9300000-9302500", human);
//    
//    auto reads = read_manager.fetch_reads(sample, region);
//    add_reads(reads, candidate_generator);
//    
//    const auto candidates = candidate_generator.generate_candidates(region);
//    
//    BOOST_REQUIRE(candidates.size() >= 15); // to make the test interesting
//    
//    HaplotypeTree haplotype_tree {"16", human};
//    
//    extend_tree(candidates, haplotype_tree);
//    
//    auto haplotypes = haplotype_tree.extract_haplotypes(region);
//    
//    sort(haplotypes);
//    
//    const Haplotype reference_haplotype {region, human};
//    
//    // remove everything other than the reference
//    const auto er = std::equal_range(haplotypes.begin(), haplotypes.end(), reference_haplotype);
//    haplotypes.erase(er.first, er.second);
//    for (const auto& haplotype : haplotypes) {
//        haplotype_tree.prune_all(haplotype);
//    }
//    
//    BOOST_REQUIRE(haplotype_tree.num_haplotypes() > 0);
//    
//    auto pruned_haplotypes = haplotype_tree.extract_haplotypes(region);
//    
//    BOOST_CHECK(std::all_of(std::cbegin(pruned_haplotypes), std::cend(pruned_haplotypes),
//                            [&] (const auto& haplotype) { return haplotype == reference_haplotype; }));
//    
//    std::sort(std::begin(pruned_haplotypes), std::end(pruned_haplotypes));
//    pruned_haplotypes.erase(std::unique(std::begin(pruned_haplotypes), std::end(pruned_haplotypes)),
//                            std::end(pruned_haplotypes));
//    
//    BOOST_REQUIRE(pruned_haplotypes.size() == 1);
//    
//    haplotype_tree.prune_unique(pruned_haplotypes.front());
//    
//    BOOST_REQUIRE(haplotype_tree.num_haplotypes() == 1);
//    
//    const auto last_haplotype = haplotype_tree.extract_haplotypes(region).front();
//    
//    BOOST_CHECK(last_haplotype == reference_haplotype);
//    
//    haplotype_tree.prune_all(last_haplotype);
//    
//    BOOST_CHECK(haplotype_tree.empty());
//}
//
//BOOST_AUTO_TEST_CASE(contains_returns_true_if_the_given_haplotype_is_in_the_tree_in_any_form)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
//    const Allele allele2 {parse_region("4:1000000-1000001", human), "C"};
//    const Allele allele3 {parse_region("4:1000000-1000001", human), "G"};
//    const Allele allele4 {parse_region("4:1000001-1000002", human), "G"};
//    const Allele allele5 {parse_region("4:1000001-1000002", human), "C"};
//    
//    const auto region = encompassing_region(allele1, allele5);
//    
//    HaplotypeTree haplotype_tree {"4", human};
//    
//    haplotype_tree.extend(allele1);
//    haplotype_tree.extend(allele2);
//    haplotype_tree.extend(allele3);
//    haplotype_tree.extend(allele4);
//    haplotype_tree.extend(allele5);
//    
//    Haplotype haplotype1 {region, human};
//    haplotype1.push_back(allele1);
//    haplotype1.push_back(allele4);
//    
//    BOOST_CHECK(haplotype_tree.contains(haplotype1));
//    
//    Haplotype haplotype2 {region, human};
//    haplotype2.push_back(allele1);
//    haplotype2.push_back(Allele {parse_region("4:1000001-1000002", human), "A"});
//    
//    BOOST_CHECK(!haplotype_tree.contains(haplotype2));
//}
//
//BOOST_AUTO_TEST_CASE(is_unique_return_true_if_the_given_haplotype_occurs_extactly_once_in_the_tree)
//{
//    // TODO
//}
//
//BOOST_AUTO_TEST_CASE(remove_can_clear_specific_regions_from_the_tree)
//{
//    // TODO
//}
//
//BOOST_AUTO_TEST_SUITE_END()
//BOOST_AUTO_TEST_SUITE_END()
