////
////  haplotype_tests.cpp
////  Octopus
////
////  Created by Daniel Cooke on 02/04/2015.
////  Copyright (c) 2015 Oxford University. All rights reserved.
////
//
//#define BOOST_TEST_DYN_LINK
//
//#include <boost/test/unit_test.hpp>
//
//#include <iostream>
//#include <string>
//
//#include "test_common.hpp"
//#include "reference_genome.hpp"
//#include "read_manager.hpp"
//#include "variant.hpp"
//#include "haplotype.hpp"
//#include "mappable_algorithms.hpp"
//
//using std::cout;
//using std::endl;
//
//BOOST_AUTO_TEST_SUITE(Components)
//BOOST_AUTO_TEST_SUITE(Haplotypes)
//
//BOOST_AUTO_TEST_CASE(empty_haplotypes_are_the_unaltered_reference_sequence)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("3:1000000-1000010", human);
//    
//    const auto ref_sequence = human.fetch_sequence(region); // CCAACAAGCA
//    
//    const Haplotype ref_haplotype {region, human};
//    
//    BOOST_REQUIRE(ref_haplotype.sequence() == ref_sequence);
//    BOOST_REQUIRE(ref_haplotype.sequence(region) == ref_sequence);
//    BOOST_REQUIRE(ref_haplotype.sequence(region.contig_region()) == ref_sequence);
//}
//
//BOOST_AUTO_TEST_CASE(alleles_can_be_added_to_front_and_back_of_haplotypes)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const Allele allele1 {"3", 1000004, "A"};
//    
//    const auto region1 = allele1.mapped_region();
//    
//    Haplotype haplotype1 {region1, human};
//    haplotype1.push_back(allele1);
//    
//    BOOST_CHECK(haplotype1.sequence(region1) == allele1.sequence());
//    BOOST_CHECK(haplotype1.sequence() == allele1.sequence());
//    
//    Haplotype haplotype2 {region1, human};
//    haplotype2.push_front(allele1);
//    
//    BOOST_CHECK(haplotype2.sequence(region1) == haplotype1.sequence());
//    BOOST_CHECK(haplotype2.sequence() == haplotype1.sequence());
//    
//    const Allele allele2 {"3", 1000005, "C"};
//    const Allele allele3 {"3", 1000006, "G"};
//    
//    const auto region2 = encompassing_region(allele1, allele3);
//    
//    Haplotype haplotype3 {region2, human};
//    haplotype3.push_back(allele1);
//    haplotype3.push_back(allele2);
//    haplotype3.push_back(allele3);
//    
//    Haplotype haplotype4 {region2, human};
//    haplotype4.push_front(allele3);
//    haplotype4.push_front(allele2);
//    haplotype4.push_front(allele1);
//    
//    BOOST_CHECK(haplotype3.sequence() == "ACG");
//    BOOST_CHECK(haplotype3.sequence(region2) == "ACG");
//    
//    BOOST_CHECK(haplotype3 == haplotype4);
//}
//
//BOOST_AUTO_TEST_CASE(regions_within_haplotype_region_flanking_explicitly_added_alleles_are_reference)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("3:1000000-1000010", human);
//    
//    const Allele allele1 {"3", 1000004, "A"};
//    
//    Haplotype haplotype1 {region, human};
//    haplotype1.push_front(allele1);
//    
//    BOOST_CHECK(haplotype1.sequence() == "CCAAAAAGCA");
//    BOOST_CHECK(haplotype1.sequence(region) == "CCAAAAAGCA");
//    
//    const Allele allele2 {"3", 1000003, "C"};
//    const Allele allele3 {"3", 1000005, "G"};
//    
//    auto haplotype2 = haplotype1;
//    haplotype2.push_front(allele2);
//    haplotype2.push_back(allele3);
//    
//    BOOST_CHECK(haplotype2.sequence() == "CCACAGAGCA");
//    BOOST_CHECK(haplotype2.sequence(region) == "CCACAGAGCA");
//    
//    const Allele allele4 {parse_region("3:1000002-1000004", human), ""};
//    const Allele allele5 {parse_region("3:1000005-1000008", human), ""};
//    
//    auto haplotype3 = haplotype1;
//    haplotype3.push_front(allele4);
//    haplotype3.push_back(allele5);
//    
//    BOOST_CHECK(haplotype3.sequence() == "CCACA");
//    BOOST_CHECK(haplotype3.sequence(region) == "CCACA");
//}
//
//BOOST_AUTO_TEST_CASE(alleles_not_explicitly_added_to_haplotypes_are_assumed_reference)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const GenomicRegion region {"7", 1000000, 1000100};
//    
//    const Haplotype reference_haplotype {region, human};
//    
//    BOOST_CHECK(reference_haplotype.contains(make_reference_allele(region, human)));
//    
//    const GenomicRegion a_sub_region {"7", 1000010, 1000090};
//    
//    BOOST_CHECK(reference_haplotype.contains(make_reference_allele(a_sub_region, human)));
//    
//    const GenomicRegion a_left_overlapping_region {"7", 999999, 1000090};
//    
//    BOOST_CHECK(!reference_haplotype.contains(make_reference_allele(a_left_overlapping_region, human)));
//    
//    const GenomicRegion a_right_overlapping_region {"7", 1000090, 1000101};
//    
//    BOOST_CHECK(!reference_haplotype.contains(make_reference_allele(a_right_overlapping_region, human)));
//}
//
//BOOST_AUTO_TEST_CASE(alleles_explicitly_added_to_haplotypes_should_be_contained)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const GenomicRegion region {"3", 1000000, 1000020}; // CCAACAAGCATTGGTGTGGC
//    
//    // Variants the haplotype will contain
//    const Variant variant1 {"3", 1000002, "A", "T"};
//    const Variant variant2 {"3", 1000004, "CA", ""};
//    const Variant variant3 {"3", 1000008, "", "C"};
//    const Variant variant4 {"3", 1000010, "TT", "GG"};
//    const Variant variant5 {"3", 1000014, "T", "C"};
//    const Variant variant6 {"3", 1000018, "G", "A"};
//    
//    // Parts of the haplotype which remain reference
//    const GenomicRegion ref_part1 {"3", 1000003, 1000004};
//    const GenomicRegion ref_part2 {"3", 1000012, 1000014};
//    const GenomicRegion ref_part3 {"3", 1000015, 1000018};
//    
//    // Variants which the haplotype will not contain
//    const Variant false_variant1 {"3", 1000002, "A", "C"};
//    const Variant false_variant2 {"3", 1000008, "", "T"};
//    const Variant false_variant3 {"3", 1000010, "TT", "AC"};
//    const Variant false_variant4 {"3", 1000014, "T", "A"};
//    
//    Haplotype haplotype {region, human};
//    
//    BOOST_CHECK(haplotype.sequence() == human.fetch_sequence(region));
//    
//    add_alt_to_back(variant1, haplotype);
//    add_alt_to_back(variant2, haplotype);
//    add_alt_to_back(variant3, haplotype);
//    add_alt_to_back(variant4, haplotype);
//    add_alt_to_back(variant5, haplotype);
//    add_alt_to_back(variant6, haplotype);
//    
//    BOOST_CHECK(haplotype.sequence() == "CCTAAGCCAGGGGCGTGAC");
//    
//    // The variant regions should be contained
//    BOOST_CHECK(contains(haplotype, variant1));
//    BOOST_CHECK(contains(haplotype, variant2));
//    BOOST_CHECK(contains(haplotype, variant3));
//    BOOST_CHECK(contains(haplotype, variant4));
//    BOOST_CHECK(contains(haplotype, variant5));
//    BOOST_CHECK(contains(haplotype, variant6));
//    
//    // But not the reference alleles (alts were added)
//    BOOST_CHECK(!haplotype.contains(variant1.ref_allele()));
//    BOOST_CHECK(!haplotype.contains(variant2.ref_allele()));
//    BOOST_CHECK(!haplotype.contains(variant3.ref_allele()));
//    BOOST_CHECK(!haplotype.contains(variant4.ref_allele()));
//    BOOST_CHECK(!haplotype.contains(variant5.ref_allele()));
//    BOOST_CHECK(!haplotype.contains(variant6.ref_allele()));
//    
//    BOOST_CHECK(!haplotype.contains(false_variant1.alt_allele()));
//    BOOST_CHECK(!haplotype.contains(false_variant2.alt_allele()));
//    BOOST_CHECK(!haplotype.contains(false_variant3.alt_allele()));
//    BOOST_CHECK(!haplotype.contains(false_variant4.alt_allele()));
//    
//    const GenomicRegion ref_begin_bit {"3", 1000000, 1000002};
//    const GenomicRegion ref_end_bit {"3", 1000019, 1000020};
//    
//    BOOST_CHECK(haplotype.contains(make_reference_allele(ref_begin_bit, human)));
//    BOOST_CHECK(haplotype.contains(make_reference_allele(ref_part1,     human)));
//    BOOST_CHECK(haplotype.contains(make_reference_allele(ref_part2,     human)));
//    BOOST_CHECK(haplotype.contains(make_reference_allele(ref_part3,     human)));
//    BOOST_CHECK(haplotype.contains(make_reference_allele(ref_end_bit,   human)));
//}
//
////BOOST_AUTO_TEST_CASE(if_a_haplotype_contains_an_insertion_then_it_should_also_contain_the_parsimonious_version)
////{
////    auto human = make_reference(human_reference_fasta);
////    
////    auto region = parse_region("3:1000000-1000020", human);
////    
////    Allele a1 {parse_region("3:1000005-1000005", human), "ACGT"};
////    Allele a2 {parse_region("3:1000004-1000005", human), "CACGT"}; // parsimonious version
////    
////    Haplotype hap1 {human, region};
////    Haplotype hap2 {human, region};
////    hap2.push_back(a1);
////    Haplotype hap3 {human, region};
////    hap3.push_back(a2);
////    
////    BOOST_CHECK(!contains(hap1, a1));
////    BOOST_CHECK(!contains(hap1, a2));
////    BOOST_CHECK(contains(hap2, a1));
////    BOOST_CHECK(contains(hap2, a1));
////    BOOST_CHECK(contains(hap3, a1));
////    BOOST_CHECK(contains(hap3, a1));
////}
//
//BOOST_AUTO_TEST_CASE(mnps_decompose)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("3:1000000-1000020", human);
//    
//    const Allele an_allele {parse_region("3:1000010-1000012", human), "GG"};
//    const Allele a_sub_allele {parse_region("3:1000010-1000011", human), "G"};
//    const Allele another_sub_allele {parse_region("3:1000011-1000012", human), "G"};
//    const Allele not_a_sub_allele {parse_region("3:1000010-1000011", human), "C"};
//    
//    Haplotype hap {region, human};
//    hap.push_back(an_allele);
//    
//    BOOST_CHECK(hap.contains(an_allele));
//    BOOST_CHECK(hap.contains(a_sub_allele));
//    BOOST_CHECK(hap.contains(another_sub_allele));
//    BOOST_CHECK(!hap.contains(not_a_sub_allele));
//}
//
//BOOST_AUTO_TEST_CASE(deletions_decompose)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("3:1000000-1000020", human);
//    
//    const Allele allele {parse_region("3:1000010-1000012", human), ""};
//    const Allele sub_allele {parse_region("3:1000010-1000011", human), ""};
//    const Allele another_sub_allele {parse_region("3:1000011-1000012", human), ""};
//    const Allele not_a_sub_allele1 {parse_region("3:1000010-1000011", human), "C"};
//    const Allele not_a_sub_allele2 {parse_region("3:1000010-1000013", human), ""};
//    
//    Haplotype hap {region, human};
//    hap.push_back(allele);
//    
//    BOOST_CHECK(hap.contains(allele));
//    BOOST_CHECK(hap.contains(sub_allele));
//    BOOST_CHECK(hap.contains(another_sub_allele));
//    BOOST_CHECK(!hap.contains(not_a_sub_allele1));
//    BOOST_CHECK(!hap.contains(not_a_sub_allele2));
//}
//
//BOOST_AUTO_TEST_CASE(insertions_decompose)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("3:1000000-1000020", human);
//    
//    const Allele allele {parse_region("3:1000010-1000010", human), "AT"};
//    const Allele sub_allele {parse_region("3:1000010-1000010", human), "A"};
//    const Allele another_sub_allele {parse_region("3:1000010-1000010", human), "T"};
//    const Allele not_a_sub_allele1 {parse_region("3:1000010-1000010", human), "C"};
//    const Allele not_a_sub_allele2 {parse_region("3:1000010-1000011", human), "A"};
//    
//    Haplotype hap {region, human};
//    hap.push_back(allele);
//    
//    BOOST_CHECK(hap.contains(allele));
//    BOOST_CHECK(hap.contains(sub_allele));
//    BOOST_CHECK(hap.contains(another_sub_allele));
//    BOOST_CHECK(!hap.contains(not_a_sub_allele1));
//    BOOST_CHECK(!hap.contains(not_a_sub_allele2));
//}
//
//BOOST_AUTO_TEST_CASE(haplotype_equate_when_alleles_infer_same_sequence)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("16:9300000-9300100", human);
//    
//    const Allele allele1 {parse_region("16:9300037-9300037", human), "TG"};
//    const Allele allele2 {parse_region("16:9300039-9300051", human), ""};
//    const Allele allele3 {parse_region("16:9300041-9300051", human), ""};
//    
//    Haplotype hap1 {region, human};
//    hap1.push_back(allele3);
//    
//    Haplotype hap2 {region, human};
//    hap2.push_back(allele1);
//    hap2.push_back(allele2);
//    
//    BOOST_CHECK(hap1.sequence() == hap2.sequence());
//    BOOST_CHECK(hap1 == hap2);
//    
//    const Allele allele4 {parse_region("16:9300037-9300038", human), "T"};
//    const Allele allele5 {parse_region("16:9300038-9300039", human), "C"};
//    const Allele allele6 {parse_region("16:9300037-9300039", human), "TC"};
//    
//    Haplotype hap3 {region, human};
//    hap3.push_back(allele4);
//    hap3.push_back(allele5);
//    
//    Haplotype hap4 {region, human};
//    hap4.push_back(allele6);
//    
//    BOOST_CHECK(hap3.sequence() == hap4.sequence());
//    BOOST_CHECK(hap3 == hap4);
//}
//
//BOOST_AUTO_TEST_CASE(haplotypes_behave_at_boundries)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("16:9299940-9300100", human);
//    
//    const Allele allele1 {parse_region("16:9299945-9299946", human), "T"};
//    const Allele allele2 {parse_region("16:9299946-9299957", human), "CGCATTACAAC"};
//    const Allele allele3 {parse_region("16:9299957-9299958", human), "C"};
//    const Allele allele4 {make_reference_allele("16:9299958-9300037", human)};
//    const Allele allele5 {parse_region("16:9300037-9300037", human), ""};
//    const Allele allele6 {parse_region("16:9300037-9300039", human), "TG"};
//    const Allele allele7 {parse_region("16:9300039-9300051", human), "TGTGTGTGCGTT"};
//    const Allele allele8 {parse_region("16:9300051-9300061", human), "TGTGTGTGTG"};
//    const Allele allele9 {parse_region("16:9300061-9300062", human), "G"};
//    const Allele allele10 {parse_region("16:9300062-9300072", human), "GTGTGTGTGT"};
//    const Allele allele11 {parse_region("16:9300072-9300073", human), "G"};
//    const Allele allele12 {parse_region("16:9300073-9300074", human), "G"};
//    const Allele allele13 {parse_region("16:9300074-9300075", human), "G"};
//    
//    Haplotype haplotype {region, human};
//    
//    haplotype.push_back(allele1);
//    haplotype.push_back(allele2);
//    haplotype.push_back(allele3);
//    haplotype.push_back(allele4);
//    haplotype.push_back(allele5);
//    haplotype.push_back(allele6);
//    haplotype.push_back(allele7);
//    haplotype.push_back(allele8);
//    haplotype.push_back(allele9);
//    haplotype.push_back(allele10);
//    haplotype.push_back(allele11);
//    haplotype.push_back(allele12);
//    haplotype.push_back(allele13);
//    
//    const Allele test_allele1 {parse_region("16:9300037-9300050", human), ""};
//    const Allele test_allele2 {parse_region("16:9300037-9300051", human), ""};
//    const Allele test_allele3 {parse_region("16:9300037-9300052", human), ""};
//    
//    BOOST_CHECK(!haplotype.contains(test_allele1));
//    BOOST_CHECK(!haplotype.contains(test_allele2));
//    BOOST_CHECK(!haplotype.contains(test_allele3));
//}
//
//BOOST_AUTO_TEST_CASE(haplotypes_can_be_copied_and_moved)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("16:9299940-9300100", human);
//    
//    const Allele allele1 {parse_region("16:9300037-9300037", human), "TG"};
//    const Allele allele2 {parse_region("16:9300039-9300051", human), ""};
//    
//    Haplotype hap {region, human};
//    hap.push_back(allele1);
//    hap.push_back(allele2);
//    
//    const auto hap_copy = hap;
//    
//    BOOST_CHECK(hap_copy.contains(allele1));
//    BOOST_CHECK(hap_copy.contains(allele2));
//    
//    const auto moved_hap = std::move(hap);
//    
//    BOOST_CHECK(moved_hap.contains(allele1));
//    BOOST_CHECK(moved_hap.contains(allele2));
//}
//
//BOOST_AUTO_TEST_CASE(get_sequence_works_at_boundries)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    Haplotype haplotype {GenomicRegion {"22", 16231802, 16232338}, human};
//    haplotype.push_back(ContigAllele {ContigRegion {16232051, 16232074}, "AAAGAGAAAGAAAGAAAGAAAGA"});
//    
//    const GenomicRegion splice_region {"22", 16232073, 16232077};
//    
//    const auto splice_haplotype = haplotype.sequence(splice_region);
//    
//    exit(0);
//}
//
//BOOST_AUTO_TEST_CASE(Haplotype_can_be_spliced)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//                                     
//    const auto allele1  = make_allele("16:9299945-9299946", "T", human);
//    const auto allele2  = make_allele("16:9299946-9299957", "CGCATTACAAC", human);
//    const auto allele3  = make_allele("16:9299957-9299958", "C", human);
//    const auto allele4  = make_reference_allele("16:9299958-9300037", human);
//    const auto allele5  = make_allele("16:9300037-9300037", "", human);
//    const auto allele6  = make_allele("16:9300037-9300039", "TG", human);
//    const auto allele7  = make_allele("16:9300039-9300051", "TGTGTGTGCGTT", human);
//    const auto allele8  = make_allele("16:9300051-9300061", "TGTGTGTGTG", human);
//    const auto allele9  = make_allele("16:9300061-9300062", "G", human);
//    const auto allele10 = make_allele("16:9300062-9300072", "GTGTGTGTGT", human);
//    const auto allele11 = make_allele("16:9300072-9300073", "G", human);
//    const auto allele12 = make_allele("16:9300073-9300074", "G", human);
//    const auto allele13 = make_allele("16:9300074-9300075", "G", human);
//    
//    std::vector<Allele> alleles {
//        allele1, allele2, allele3, allele4, allele5, allele6, allele7, allele8, allele9, allele10,
//        allele11, allele12, allele13
//    };
//    
//    const auto region = parse_region("16:9299940-9300100", human);
//    
//    Haplotype haplotype {region, human};
//    
//    for (const auto& allele : alleles) haplotype.push_back(allele);
//    
//    const auto splice_region = parse_region("16:9299964-9300083", human);
//    
//    const auto spliced = splice<Haplotype>(haplotype, splice_region);
//    
//    BOOST_CHECK(is_same_region(spliced, splice_region));
//    BOOST_CHECK(contains(haplotype, spliced));
//    
//    BOOST_CHECK(splice<Haplotype>(haplotype, region) == haplotype);
//    
//    const auto allele_region = encompassing_region(alleles);
//    
//    const auto allele_region_splice = splice<Haplotype>(haplotype, allele_region);
//    
//    BOOST_CHECK(std::all_of(std::cbegin(alleles), std::cend(alleles),
//                            [&] (const auto& allele) {
//                                return allele_region_splice.contains_exact(allele); })
//                );
//    
//    const auto partial_allele_region = parse_region("16:9299947-9300074", human);
//    
//    const auto partial_allele_region_splice = splice<Haplotype>(haplotype, partial_allele_region);
//    
//    BOOST_CHECK(!partial_allele_region_splice.contains(allele1));
//    BOOST_CHECK(!partial_allele_region_splice.contains_exact(allele2));
//    BOOST_CHECK(!partial_allele_region_splice.contains(allele13));
//    BOOST_CHECK(std::all_of(std::next(std::cbegin(alleles), 2), std::prev(std::cend(alleles)),
//                            [&] (const auto& allele) {
//                                return allele_region_splice.contains_exact(allele); })
//                );
//    
//    const auto insertion_region = mapped_region(allele5);
//    
//    const auto insertion_region_haplotype_splice = splice<Haplotype>(haplotype, insertion_region);
//    
//    BOOST_CHECK(is_same_region(insertion_region_haplotype_splice, insertion_region));
//    BOOST_CHECK(insertion_region_haplotype_splice.sequence() == "");
//    BOOST_CHECK(insertion_region_haplotype_splice.contains_exact(allele5));
//    
//    const auto insertion_region_allele_splice = splice<Allele>(haplotype, insertion_region);
//    
//    BOOST_CHECK(insertion_region_allele_splice == allele5);
//    
//    const auto insertion_before_snp_region = head_region(allele9);
//    
//    const auto insertion_before_snp_haplotype_splice = splice<Haplotype>(haplotype, insertion_before_snp_region);
//    
//    BOOST_CHECK(is_same_region(insertion_before_snp_haplotype_splice, insertion_before_snp_region));
//    BOOST_CHECK(insertion_before_snp_haplotype_splice.sequence() == "");
//    BOOST_CHECK(is_reference(insertion_before_snp_haplotype_splice));
//    
//    const auto insertion_before_snp_allele_splice = splice<Allele>(haplotype, insertion_before_snp_region);
//    
//    BOOST_CHECK(is_empty_region(insertion_before_snp_allele_splice));
//    BOOST_CHECK(is_reference(insertion_before_snp_allele_splice, human));
//}
//
//BOOST_AUTO_TEST_CASE(splicing_around_insertions_works)
//{
//    using std::cbegin; using std::cend; using std::for_each; using std::next;
//    
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto variant0 = make_variant("6:31235411-31235411", "A", human);
//    const auto variant1 = make_variant("6:31235411-31235412", "A", human);
//    const auto variant2 = make_variant("6:31235412-31235413", "",  human);
//    const auto variant3 = make_variant("6:31235413-31235414", "A", human);
//    const auto variant4 = make_variant("6:31235414-31235414", "A", human);
//    const auto variant5 = make_variant("6:31235414-31235415", "A", human);
//    const auto variant6 = make_variant("6:31235415-31235415", "A", human);
//    
//    const std::vector<Variant> variants {variant0, variant1, variant2, variant3, variant4, variant5, variant6};
//    
//    const auto region = expand(encompassing_region(variants), 10);
//    
//    Haplotype ref_haplotype {region, human};
//    for (const auto& variant : variants) add_ref_to_back(variant, ref_haplotype);
//    
//    Haplotype alt_haplotype {region, human};
//    for (const auto& variant : variants) add_alt_to_back(variant, alt_haplotype);
//    
//    // front
//    
//    // TODO: before insertion
//    
//    // on insertion
//    
//    auto splice_hap_ref0 = splice<Haplotype>(ref_haplotype, mapped_region(variant0));
//    auto splice_hap_alt0 = splice<Haplotype>(alt_haplotype, mapped_region(variant0));
//    
//    BOOST_CHECK(is_same_region(splice_hap_ref0, variant0));
//    BOOST_CHECK(splice_hap_ref0.sequence() == "");
//    BOOST_CHECK(is_same_region(splice_hap_alt0, variant0));
//    BOOST_CHECK(splice_hap_alt0.sequence() == "A");
//    BOOST_CHECK(splice_hap_ref0.difference(splice_hap_alt0).size() == 1);
//    
//    // TODO: after insertion
//    
//    // middle
//    
//    // TODO: before insertion
//    
//    // on insertion
//    
//    auto splice_hap_ref4 = splice<Haplotype>(ref_haplotype, mapped_region(variant4));
//    auto splice_hap_alt4 = splice<Haplotype>(alt_haplotype, mapped_region(variant4));
//    
//    BOOST_CHECK(is_same_region(splice_hap_ref4, variant4));
//    BOOST_CHECK(splice_hap_ref4.sequence() == "");
//    BOOST_CHECK(is_same_region(splice_hap_alt4, variant4));
//    BOOST_CHECK(splice_hap_alt4.sequence() == "A");
//    BOOST_CHECK(splice_hap_ref4.difference(splice_hap_alt4).size() == 1);
//    
//    // TODO: after insertion
//    
//    // back
//    
//    // TODO: before insertion
//    
//    // on insertion
//    
//    auto splice_hap_ref6 = splice<Haplotype>(ref_haplotype, mapped_region(variant6));
//    auto splice_hap_alt6 = splice<Haplotype>(alt_haplotype, mapped_region(variant6));
//    
//    BOOST_CHECK(is_same_region(splice_hap_ref6, variant6));
//    BOOST_CHECK(splice_hap_ref6.sequence() == "");
//    BOOST_CHECK(is_same_region(splice_hap_alt6, variant6));
//    BOOST_CHECK(splice_hap_alt6.sequence() == "A");
//    BOOST_CHECK(splice_hap_ref6.difference(splice_hap_alt6).size() == 1);
//    
//    // TODO: after insertion
//    
//    // TODO
//}
//
//BOOST_AUTO_TEST_SUITE_END() // Haplotypes
//BOOST_AUTO_TEST_SUITE_END() // Components
