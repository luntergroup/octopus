//
//  haplotype_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 02/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <cstddef>
#include <set>

#include "test_common.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "variant.hpp"
#include "variant_utils.hpp"
#include "candidate_variant_generator.hpp"
#include "alignment_candidate_variant_generator.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "mappable_algorithms.hpp"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(Components)

static void add_alt_to_front(const Variant& variant, Haplotype& haplotype)
{
    haplotype.push_front(variant.get_alt_allele());
}

static void add_alt_to_back(const Variant& variant, Haplotype& haplotype)
{
    haplotype.push_back(variant.get_alt_allele());
}

BOOST_AUTO_TEST_CASE(alleles_can_be_added_to_front_and_back_of_haplotypes)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    auto region = *parse_region("3:1000000-1000010", human);
    
    auto ref_sequence = human.get_sequence(region); // CCAACAAGCA
    
    Haplotype ref_haplotype {region, human};
    
    BOOST_CHECK(ref_haplotype.get_sequence() == ref_sequence);
    
    Allele allele1 {"3", 1000004, "A"};
    
    Haplotype haplotype2 {region, human};
    haplotype2.push_front(allele1);
    
    BOOST_CHECK(haplotype2.get_sequence(region) == "CCAAAAAGCA");
    
    Haplotype haplotype3 {region, human};
    haplotype3.push_back(allele1);
    
    BOOST_CHECK(haplotype2.get_sequence(region) == haplotype3.get_sequence(region));
    BOOST_CHECK(haplotype2.get_sequence() == haplotype3.get_sequence());
    
    Variant variant2 {"3", 1000004, "CA", ""};
    Variant variant3 {"3", 1000008, "", "C"};
    
    Haplotype haplotype4 {region, human};
    add_alt_to_back(variant2, haplotype4);
    add_alt_to_back(variant3, haplotype4);
    
    BOOST_CHECK(haplotype4.get_sequence() == "CCAAAGCCA");
    
    Haplotype haplotype5 {region, human};
    add_alt_to_back(variant2, haplotype5);
    add_alt_to_back(variant3, haplotype5);
    
    BOOST_CHECK(haplotype5.get_sequence() == "AGC");
    
    Variant variant4 {"3", 1000004, "CA", "GG"};
    
    Haplotype haplotype6 {region, human};
    add_alt_to_back(variant4, haplotype6);
    
    BOOST_CHECK(haplotype6.get_sequence() == "CCAAGGAGCA");
    
    Variant variant5 {"3", 1000004, "C", "G"};
    Variant variant6 {"3", 1000005, "A", "G"};
    
    Haplotype haplotype7 {region, human};
    add_alt_to_back(variant6, haplotype7);
    add_alt_to_front(variant5, haplotype7);
    
    BOOST_CHECK(haplotype7.get_sequence() == haplotype6.get_sequence());
}

BOOST_AUTO_TEST_CASE(haplotypes_work_with_real_data)
{
    BOOST_REQUIRE(test_file_exists(ecoli_reference_fasta));
    
    auto ecoli = make_reference(ecoli_reference_fasta);
    
    ReadManager read_manager {ecoli_bam};
    
    Octopus::CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<Octopus::AlignmentCandidateVariantGenerator>(ecoli, 0));
    
    auto region = *parse_region("R00000042:99640-99745", ecoli);
    
    auto reference_sequence = ecoli.get_sequence(region);
    
    const auto sample = read_manager.get_samples().front();
    
    auto some_reads = read_manager.fetch_reads(sample, region);
    
    add_reads(some_reads, candidate_generator);
    
    auto variants = candidate_generator.get_candidates(region);
    
    BOOST_CHECK(variants.size() == 12);
    
    Haplotype haplotype1 {region, ecoli};
    for (const auto& variant : variants) {
        if (is_snp(variant)) {
            add_alt_to_back(variant, haplotype1);
        }
    }
    
    BOOST_CHECK(haplotype1.get_sequence() == "AGCGTGGGTAAACAAAGCCATGCTATCAGCACCGCCAGCGGCGTTGGCGAACA"
            "TTTTGCTGATAAAACTGCGTTAATTACGCGTCTTAAATTACTGATTGCTGAG");
}

BOOST_AUTO_TEST_CASE(alleles_not_explicitly_added_to_haplotypes_are_assumed_reference)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    GenomicRegion region {"7", 1000000, 1000100};
    
    Haplotype reference_haplotype {region, human};
    
    BOOST_CHECK(reference_haplotype.contains(get_reference_allele(region, human)));
    
    GenomicRegion a_sub_region {"7", 1000010, 1000090};
    
    BOOST_CHECK(reference_haplotype.contains(get_reference_allele(a_sub_region, human)));
    
    GenomicRegion a_left_overlapping_region {"7", 999999, 1000090};
    
    BOOST_CHECK(!reference_haplotype.contains(get_reference_allele(a_left_overlapping_region, human)));
    
    GenomicRegion a_right_overlapping_region {"7", 1000090, 1000101};
    
    BOOST_CHECK(!reference_haplotype.contains(get_reference_allele(a_right_overlapping_region, human)));
}

BOOST_AUTO_TEST_CASE(alleles_explicitly_added_to_haplotypes_should_be_contained)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    GenomicRegion region {"3", 1000000, 1000020}; // CCAACAAGCATTGGTGTGGC
    
    // Variants the haplotype will contain
    Variant variant1 {"3", 1000002, "A", "T"};
    Variant variant2 {"3", 1000004, "CA", ""};
    Variant variant3 {"3", 1000008, "", "C"};
    Variant variant4 {"3", 1000010, "TT", "GG"};
    Variant variant5 {"3", 1000014, "T", "C"};
    Variant variant6 {"3", 1000018, "G", "A"};
    
    // Parts of the haplotype which remain reference
    GenomicRegion ref_part1 {"3", 1000003, 1000004};
    GenomicRegion ref_part2 {"3", 1000012, 1000014};
    GenomicRegion ref_part3 {"3", 1000015, 1000018};
    
    // Variants which the haplotype does not contain
    Variant false_variant1 {"3", 1000002, "A", "C"};
    Variant false_variant2 {"3", 1000008, "", "T"};
    Variant false_variant3 {"3", 1000010, "TT", "AC"};
    Variant false_variant4 {"3", 1000014, "T", "A"};
    
    Haplotype haplotype_unbounded {region, human};
    add_alt_to_back(variant1, haplotype_unbounded);
    add_alt_to_back(variant2, haplotype_unbounded);
    add_alt_to_back(variant3, haplotype_unbounded);
    add_alt_to_back(variant4, haplotype_unbounded);
    add_alt_to_back(variant5, haplotype_unbounded);
    add_alt_to_back(variant6, haplotype_unbounded);
    
    BOOST_CHECK(haplotype_unbounded.get_sequence() == "TAAGCCAGGGGCGTGA");
    
    BOOST_CHECK(contains(haplotype_unbounded, variant1));
    BOOST_CHECK(contains(haplotype_unbounded, variant2));
    BOOST_CHECK(contains(haplotype_unbounded, variant3));
    BOOST_CHECK(contains(haplotype_unbounded, variant4));
    BOOST_CHECK(contains(haplotype_unbounded, variant5));
    BOOST_CHECK(contains(haplotype_unbounded, variant6));
    
    BOOST_CHECK(!haplotype_unbounded.contains(variant1.get_ref_allele()));
    BOOST_CHECK(!haplotype_unbounded.contains(variant2.get_ref_allele()));
    BOOST_CHECK(!haplotype_unbounded.contains(variant3.get_ref_allele()));
    BOOST_CHECK(!haplotype_unbounded.contains(variant4.get_ref_allele()));
    BOOST_CHECK(!haplotype_unbounded.contains(variant5.get_ref_allele()));
    BOOST_CHECK(!haplotype_unbounded.contains(variant6.get_ref_allele()));
    
    BOOST_CHECK(!contains(haplotype_unbounded, false_variant1));
    BOOST_CHECK(!contains(haplotype_unbounded, false_variant2));
    BOOST_CHECK(!contains(haplotype_unbounded, false_variant3));
    BOOST_CHECK(!contains(haplotype_unbounded, false_variant4));
    
    BOOST_CHECK(haplotype_unbounded.contains(get_reference_allele(ref_part1, human)));
    BOOST_CHECK(haplotype_unbounded.contains(get_reference_allele(ref_part2, human)));
    BOOST_CHECK(haplotype_unbounded.contains(get_reference_allele(ref_part3, human)));
    
    Haplotype haplotype_bounded {region, human};
    
    BOOST_CHECK(haplotype_bounded.get_sequence() == human.get_sequence(region));
    
    add_alt_to_back(variant1, haplotype_bounded);
    add_alt_to_back(variant2, haplotype_bounded);
    add_alt_to_back(variant3, haplotype_bounded);
    add_alt_to_back(variant4, haplotype_bounded);
    add_alt_to_back(variant5, haplotype_bounded);
    add_alt_to_back(variant6, haplotype_bounded);
    
    BOOST_CHECK(haplotype_bounded.get_sequence() == "CCTAAGCCAGGGGCGTGAC");
    
    BOOST_CHECK(contains(haplotype_bounded, variant1));
    BOOST_CHECK(contains(haplotype_bounded, variant2));
    BOOST_CHECK(contains(haplotype_bounded, variant3));
    BOOST_CHECK(contains(haplotype_bounded, variant4));
    BOOST_CHECK(contains(haplotype_bounded, variant5));
    BOOST_CHECK(contains(haplotype_bounded, variant6));
    
    BOOST_CHECK(!haplotype_bounded.contains(variant1.get_ref_allele()));
    BOOST_CHECK(!haplotype_bounded.contains(variant2.get_ref_allele()));
    BOOST_CHECK(!haplotype_bounded.contains(variant3.get_ref_allele()));
    BOOST_CHECK(!haplotype_bounded.contains(variant4.get_ref_allele()));
    BOOST_CHECK(!haplotype_bounded.contains(variant5.get_ref_allele()));
    BOOST_CHECK(!haplotype_bounded.contains(variant6.get_ref_allele()));
    
    BOOST_CHECK(!contains(haplotype_bounded, false_variant1));
    BOOST_CHECK(!contains(haplotype_bounded, false_variant2));
    BOOST_CHECK(!contains(haplotype_bounded, false_variant3));
    BOOST_CHECK(!contains(haplotype_bounded, false_variant4));
    
    GenomicRegion ref_begin_bit {"3", 1000000, 1000002};
    GenomicRegion ref_end_bit {"3", 1000019, 1000020};
    
    BOOST_CHECK(haplotype_bounded.contains(get_reference_allele(ref_begin_bit, human)));
    BOOST_CHECK(haplotype_bounded.contains(get_reference_allele(ref_part1,     human)));
    BOOST_CHECK(haplotype_bounded.contains(get_reference_allele(ref_part2,     human)));
    BOOST_CHECK(haplotype_bounded.contains(get_reference_allele(ref_part3,     human)));
    BOOST_CHECK(haplotype_bounded.contains(get_reference_allele(ref_end_bit,   human)));
}

//BOOST_AUTO_TEST_CASE(if_a_haplotype_contains_an_insertion_then_it_should_also_contain_the_parsimonious_version)
//{
//    auto human = make_reference(human_reference_fasta);
//    
//    auto region = parse_region("3:1000000-1000020", human);
//    
//    Allele a1 {*parse_region("3:1000005-1000005", human), "ACGT"};
//    Allele a2 {*parse_region("3:1000004-1000005", human), "CACGT"}; // parsimonious version
//    
//    Haplotype hap1 {human, region};
//    Haplotype hap2 {human, region};
//    hap2.push_back(a1);
//    Haplotype hap3 {human, region};
//    hap3.push_back(a2);
//    
//    BOOST_CHECK(!contains(hap1, a1));
//    BOOST_CHECK(!contains(hap1, a2));
//    BOOST_CHECK(contains(hap2, a1));
//    BOOST_CHECK(contains(hap2, a1));
//    BOOST_CHECK(contains(hap3, a1));
//    BOOST_CHECK(contains(hap3, a1));
//}

BOOST_AUTO_TEST_CASE(mnps_decompose)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    auto region = *parse_region("3:1000000-1000020", human);
    
    Allele an_allele {*parse_region("3:1000010-1000012", human), "GG"};
    
    Allele a_sub_allele {*parse_region("3:1000010-1000011", human), "G"};
    Allele another_sub_allele {*parse_region("3:1000011-1000012", human), "G"};
    Allele not_a_sub_allele {*parse_region("3:1000010-1000011", human), "C"};
    
    Haplotype hap {region, human};
    hap.push_back(an_allele);
    
    BOOST_CHECK(hap.contains(an_allele));
    BOOST_CHECK(hap.contains(a_sub_allele));
    BOOST_CHECK(hap.contains(another_sub_allele));
    BOOST_CHECK(!hap.contains(not_a_sub_allele));
}

BOOST_AUTO_TEST_CASE(deletions_decompose)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    auto region = *parse_region("3:1000000-1000020", human);
    
    Allele allele {*parse_region("3:1000010-1000012", human), ""};
    
    Allele sub_allele {*parse_region("3:1000010-1000011", human), ""};
    Allele another_sub_allele {*parse_region("3:1000011-1000012", human), ""};
    Allele not_a_sub_allele1 {*parse_region("3:1000010-1000011", human), "C"};
    Allele not_a_sub_allele2 {*parse_region("3:1000010-1000013", human), ""};
    
    Haplotype hap {region, human};
    hap.push_back(allele);
    
    BOOST_CHECK(hap.contains(allele));
    BOOST_CHECK(hap.contains(sub_allele));
    BOOST_CHECK(hap.contains(another_sub_allele));
    BOOST_CHECK(!hap.contains(not_a_sub_allele1));
    BOOST_CHECK(!hap.contains(not_a_sub_allele2));
}

BOOST_AUTO_TEST_CASE(insertions_decompose)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    auto region = *parse_region("3:1000000-1000020", human);
    
    Allele allele {*parse_region("3:1000010-1000010", human), "AT"};
    
    Allele sub_allele {*parse_region("3:1000010-1000010", human), "A"};
    Allele another_sub_allele {*parse_region("3:1000010-1000010", human), "T"};
    Allele not_a_sub_allele1 {*parse_region("3:1000010-1000010", human), "C"};
    Allele not_a_sub_allele2 {*parse_region("3:1000010-1000011", human), "A"};
    
    Haplotype hap {region, human};
    hap.push_back(allele);
    
    BOOST_CHECK(hap.contains(allele));
    BOOST_CHECK(hap.contains(sub_allele));
    BOOST_CHECK(hap.contains(another_sub_allele));
    BOOST_CHECK(!hap.contains(not_a_sub_allele1));
    BOOST_CHECK(!hap.contains(not_a_sub_allele2));
}

BOOST_AUTO_TEST_CASE(haplotype_equate_when_alleles_infer_same_sequence)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    auto region = *parse_region("16:9300000-9300100", human);
    
    Allele allele1 {*parse_region("16:9300037-9300037", human), "TG"};
    Allele allele2 {*parse_region("16:9300039-9300051", human), ""};
    Allele allele3 {*parse_region("16:9300041-9300051", human), ""};
    
    Haplotype hap1 {region, human};
    hap1.push_back(allele3);
    
    Haplotype hap2 {region, human};
    hap2.push_back(allele1);
    hap2.push_back(allele2);
    
    BOOST_CHECK(hap1.get_sequence() == hap2.get_sequence());
    BOOST_CHECK(hap1 == hap2);
    
    Allele allele4 {*parse_region("16:9300037-9300038", human), "T"};
    Allele allele5 {*parse_region("16:9300038-9300039", human), "C"};
    Allele allele6 {*parse_region("16:9300037-9300039", human), "TC"};
    
    Haplotype hap3 {region, human};
    hap3.push_back(allele4);
    hap3.push_back(allele5);
    
    Haplotype hap4 {region, human};
    hap4.push_back(allele6);
    
    BOOST_CHECK(hap3.get_sequence() == hap4.get_sequence());
    BOOST_CHECK(hap3 == hap4);
}

BOOST_AUTO_TEST_CASE(haplotypes_can_be_compared_for_structural_complexity)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    auto region = *parse_region("16:9300000-9300100", human);
    
    Allele allele1 {*parse_region("16:9300037-9300037", human), "TG"};
    Allele allele2 {*parse_region("16:9300039-9300051", human), ""};
    Allele allele3 {*parse_region("16:9300041-9300051", human), ""};
    
    Haplotype hap1 {region, human};
    hap1.push_back(allele3);
    
    Haplotype hap2 {region, human};
    hap2.push_back(allele1);
    hap2.push_back(allele2);
    
    BOOST_CHECK(IsLessComplex()(hap1, hap2));
}

BOOST_AUTO_TEST_CASE(haplotypes_behave_at_boundries)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    auto region = *parse_region("16:9299940-9300100", human);
    
    Allele allele1 {*parse_region("16:9299945-9299946", human), "T"};
    Allele allele2 {*parse_region("16:9299946-9299957", human), "CGCATTACAAC"};
    Allele allele3 {*parse_region("16:9299957-9299958", human), "C"};
    Allele allele4 {get_reference_allele(*parse_region("16:9299958-9300037", human), human)};
    Allele allele5 {*parse_region("16:9300037-9300037", human), ""};
    Allele allele6 {*parse_region("16:9300037-9300039", human), "TG"};
    Allele allele7 {*parse_region("16:9300039-9300051", human), "TGTGTGTGCGTT"};
    Allele allele8 {*parse_region("16:9300051-9300061", human), "TGTGTGTGTG"};
    Allele allele9 {*parse_region("16:9300061-9300062", human), "G"};
    Allele allele10 {*parse_region("16:9300062-9300072", human), "GTGTGTGTGT"};
    Allele allele11 {*parse_region("16:9300072-9300073", human), "G"};
    Allele allele12 {*parse_region("16:9300073-9300074", human), "G"};
    Allele allele13 {*parse_region("16:9300074-9300075", human), "G"};
    
    Haplotype haplotype {region, human};
    
    haplotype.push_back(allele1);
    haplotype.push_back(allele2);
    haplotype.push_back(allele3);
    haplotype.push_back(allele4);
    haplotype.push_back(allele5);
    haplotype.push_back(allele6);
    haplotype.push_back(allele7);
    haplotype.push_back(allele8);
    haplotype.push_back(allele9);
    haplotype.push_back(allele10);
    haplotype.push_back(allele11);
    haplotype.push_back(allele12);
    haplotype.push_back(allele13);
    
    Allele test_allele1 {*parse_region("16:9300037-9300050", human), ""};
    Allele test_allele2 {*parse_region("16:9300037-9300051", human), ""};
    Allele test_allele3 {*parse_region("16:9300037-9300052", human), ""};
    
    BOOST_CHECK(!haplotype.contains(test_allele1));
    BOOST_CHECK(!haplotype.contains(test_allele2));
    BOOST_CHECK(!haplotype.contains(test_allele3));
}

BOOST_AUTO_TEST_CASE(haplotypes_can_be_copied_and_moved)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    auto region = *parse_region("16:9299940-9300100", human);
    
    Allele allele1 {*parse_region("16:9300037-9300037", human), "TG"};
    Allele allele2 {*parse_region("16:9300039-9300051", human), ""};
    
    Haplotype hap {region, human};
    hap.push_back(allele1);
    hap.push_back(allele2);
    
    auto hap_copy = hap;
    
    BOOST_CHECK(hap_copy.contains(allele1));
    BOOST_CHECK(hap_copy.contains(allele2));
    
    auto moved_hap = std::move(hap);
    
    BOOST_CHECK(moved_hap.contains(allele1));
    BOOST_CHECK(moved_hap.contains(allele2));
}

BOOST_AUTO_TEST_CASE(Haplotype_can_be_spliced)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    auto human = make_reference(human_reference_fasta);
    
    auto region = *parse_region("16:9299940-9300100", human);
    
    Allele allele1 {*parse_region("16:9299945-9299946", human), "T"};
    Allele allele2 {*parse_region("16:9299946-9299957", human), "CGCATTACAAC"};
    Allele allele3 {*parse_region("16:9299957-9299958", human), "C"};
    Allele allele4 {get_reference_allele(*parse_region("16:9299958-9300037", human), human)};
    Allele allele5 {*parse_region("16:9300037-9300037", human), ""};
    Allele allele6 {*parse_region("16:9300037-9300039", human), "TG"};
    Allele allele7 {*parse_region("16:9300039-9300051", human), "TGTGTGTGCGTT"};
    Allele allele8 {*parse_region("16:9300051-9300061", human), "TGTGTGTGTG"};
    Allele allele9 {*parse_region("16:9300061-9300062", human), "G"};
    Allele allele10 {*parse_region("16:9300062-9300072", human), "GTGTGTGTGT"};
    Allele allele11 {*parse_region("16:9300072-9300073", human), "G"};
    Allele allele12 {*parse_region("16:9300073-9300074", human), "G"};
    Allele allele13 {*parse_region("16:9300074-9300075", human), "G"};
    
    Haplotype haplotype {region, human};
    
    haplotype.push_back(allele1);
    haplotype.push_back(allele2);
    haplotype.push_back(allele3);
    haplotype.push_back(allele4);
    haplotype.push_back(allele5);
    haplotype.push_back(allele6);
    haplotype.push_back(allele7);
    haplotype.push_back(allele8);
    haplotype.push_back(allele9);
    haplotype.push_back(allele10);
    haplotype.push_back(allele11);
    haplotype.push_back(allele12);
    haplotype.push_back(allele13);
    
    auto splice_region = *parse_region("16:9299964-9300083", human);
    
    auto spliced = splice<Haplotype>(haplotype, splice_region);
    
    BOOST_CHECK(get_region(spliced) == splice_region);
    BOOST_CHECK(contains(haplotype, spliced));
}

//BOOST_AUTO_TEST_CASE(unique_least_complex removes haplotypes that infer the same sequence,"
//          "leaving the haplotype with the fewest alterations to the reference)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    auto region = *parse_region("16:9300000-9300100", human);
//    
//    Allele allele1 {*parse_region("16:9300037-9300037", human), "TG"};
//    Allele allele2 {*parse_region("16:9300037-9300051", human), ""};
//    Allele allele3 {*parse_region("16:9300051-9300051", human), ""};
//    Allele allele4 {*parse_region("16:9300051-9300061", human), "TGTGTGTGTG"};
//    Allele allele5 {*parse_region("16:9300061-9300062", human), "T"};
//    Allele allele6 {*parse_region("16:9300062-9300072", human), "GTGTGTGTGT"};
//    Allele allele7 {*parse_region("16:9300072-9300073", human), "T"};
//    Allele allele8 {*parse_region("16:9300073-9300074", human), "G"};
//    Allele allele9 {*parse_region("16:9300074-9300075", human), "G"};
//    
//    Allele allele10 {*parse_region("16:9300037-9300039", human), "TG"};
//    Allele allele11 {*parse_region("16:9300039-9300051", human), ""};
//    
//    Haplotype hap1 {human, region};
//    hap1.push_back(allele1);
//    hap1.push_back(allele2);
//    hap1.push_back(allele3);
//    hap1.push_back(allele4);
//    hap1.push_back(allele5);
//    hap1.push_back(allele6);
//    hap1.push_back(allele7);
//    hap1.push_back(allele8);
//    hap1.push_back(allele9);
//    
//    Haplotype hap2 {human, region};
//    hap2.push_back(allele10);
//    hap2.push_back(allele11);
//    hap2.push_back(allele3);
//    hap2.push_back(allele4);
//    hap2.push_back(allele5);
//    hap2.push_back(allele6);
//    hap2.push_back(allele7);
//    hap2.push_back(allele8);
//    hap2.push_back(allele9);
//    
//    BOOST_CHECK(hap1 == hap2);
//    
//    std::vector<Haplotype> haplotypes {hap1, hap2};
//    
//    unique_least_complex(haplotypes);
//    
//    BOOST_CHECK(haplotypes.size() == 1);
//    BOOST_CHECK(!haplotypes[0].contains(allele2));
//}

BOOST_AUTO_TEST_SUITE_END()
