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

BOOST_AUTO_TEST_CASE(alleles_can_be_added_to_front_and_back_of_haplotypes)
{
    auto human = make_reference(human_reference_fasta);
    
    auto region = parse_region("3:1000000-1000010", human);
    
    auto the_reference_sequence = human.get_sequence(region); // CCAACAAGCA
    
    Haplotype reference_haplotype {human, region};
    
    BOOST_CHECK(reference_haplotype.get_sequence() == the_reference_sequence);
    
    Variant variant_1 {"3", 1000004, "C", "A"};
    
    Haplotype haplotype2 {human};
    add_to_back(variant_1, haplotype2);
    
    BOOST_CHECK(haplotype2.get_sequence(region) == "CCAAAAAGCA");
    
    Haplotype haplotype3 {human, region};
    add_to_back(variant_1, haplotype3);
    
    BOOST_CHECK(haplotype2.get_sequence(region) == haplotype3.get_sequence());
    
    Variant variant_2 {"3", 1000004, "CA", ""};
    Variant variant_3 {"3", 1000008, "", "C"};
    
    Haplotype haplotype4 {human, region};
    add_to_back(variant_2, haplotype4);
    add_to_back(variant_3, haplotype4);
    
    BOOST_CHECK(haplotype4.get_sequence() == "CCAAAGCCA");
    
    Haplotype haplotype5 {human};
    add_to_back(variant_2, haplotype5);
    add_to_back(variant_3, haplotype5);
    
    BOOST_CHECK(haplotype5.get_sequence() == "AGC");
    
    Variant variant_4 {"3", 1000004, "CA", "GG"};
    
    Haplotype haplotype6 {human, region};
    add_to_back(variant_4, haplotype6);
    
    BOOST_CHECK(haplotype6.get_sequence() == "CCAAGGAGCA");
    
    Variant variant_5 {"3", 1000004, "C", "G"};
    Variant variant_6 {"3", 1000005, "A", "G"};
    
    Haplotype haplotype7 {human, region};
    add_to_back(variant_6, haplotype7);
    add_to_front(variant_5, haplotype7);
    
    BOOST_CHECK(haplotype7.get_sequence() == haplotype6.get_sequence());
}

BOOST_AUTO_TEST_CASE(haplotypes_work_with_real_data)
{
    auto ecoli = make_reference(ecoli_reference_fasta);
    
    ReadManager a_read_manager {ecoli_bam};
    
    Octopus::CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<Octopus::AlignmentCandidateVariantGenerator>(ecoli, 0));
    
    auto region = parse_region("R00000042:99640-99745", ecoli);
    
    auto reference_sequence = ecoli.get_sequence(region);
    
    auto sample_ids = a_read_manager.get_samples();
    auto the_sample_id = sample_ids.at(0);
    
    auto some_reads = a_read_manager.fetch_reads(the_sample_id, region);
    
    candidate_generator.add_reads(some_reads.begin(), some_reads.end());
    
    auto variants = candidate_generator.get_candidates(region);
    
    BOOST_CHECK(variants.size() == 12);
    
    Haplotype haplotype1 {ecoli, region};
    for (const auto& variant : variants) {
        if (is_snp(variant)) {
            add_to_back(variant, haplotype1);
        }
    }
    
    BOOST_CHECK(haplotype1.get_sequence() == "AGCGTGGGTAAACAAAGCCATGCTATCAGCACCGCCAGCGGCGTTGGCGAACA"
            "TTTTGCTGATAAAACTGCGTTAATTACGCGTCTTAAATTACTGATTGCTGAG");
}

BOOST_AUTO_TEST_CASE(alleles_not_explicitly_added_to_haplotypes_are_assumed_reference)
{
    auto human = make_reference(human_reference_fasta);
    
    GenomicRegion region {"7", 1000000, 1000100};
    
    Haplotype a_reference_haplotype {human, region};
    
    BOOST_CHECK(a_reference_haplotype.contains(get_reference_allele(region, human)));
    
    GenomicRegion a_sub_region {"7", 1000010, 1000090};
    
    BOOST_CHECK(a_reference_haplotype.contains(get_reference_allele(a_sub_region, human)));
    
    GenomicRegion a_left_overlapping_region {"7", 999999, 1000090};
    
    BOOST_CHECK(!a_reference_haplotype.contains(get_reference_allele(a_left_overlapping_region, human)));
    
    GenomicRegion a_right_overlapping_region {"7", 1000090, 1000101};
    
    BOOST_CHECK(!a_reference_haplotype.contains(get_reference_allele(a_right_overlapping_region, human)));
}

BOOST_AUTO_TEST_CASE(alleles_explicitly_added_to_haplotypes_should_be_contained)
{
    auto human = make_reference(human_reference_fasta);
    
    GenomicRegion region {"3", 1000000, 1000020}; // CCAACAAGCATTGGTGTGGC
    
    // Variants the haplotype will contain
    Variant variant_1 {"3", 1000002, "A", "T"};
    Variant variant_2 {"3", 1000004, "CA", ""};
    Variant variant_3 {"3", 1000008, "", "C"};
    Variant variant_4 {"3", 1000010, "TT", "GG"};
    Variant variant_5 {"3", 1000014, "T", "C"};
    Variant variant_6 {"3", 1000018, "G", "A"};
    
    // Parts of the haplotype which remain reference
    GenomicRegion a_reference_part1 {"3", 1000003, 1000004};
    GenomicRegion a_reference_part2 {"3", 1000012, 1000014};
    GenomicRegion a_reference_part3 {"3", 1000015, 1000018};
    
    // Variants which the haplotype does not contain
    Variant false_variant_1 {"3", 1000002, "A", "C"};
    Variant false_variant_2 {"3", 1000008, "", "T"};
    Variant false_variant_3 {"3", 1000010, "TT", "AC"};
    Variant false_variant_4 {"3", 1000014, "T", "A"};
    
    Haplotype haplotype_unbounded {human};
    add_to_back(variant_1, haplotype_unbounded);
    add_to_back(variant_2, haplotype_unbounded);
    add_to_back(variant_3, haplotype_unbounded);
    add_to_back(variant_4, haplotype_unbounded);
    add_to_back(variant_5, haplotype_unbounded);
    add_to_back(variant_6, haplotype_unbounded);
    
    BOOST_CHECK(haplotype_unbounded.get_sequence() == "TAAGCCAGGGGCGTGA");
    
    BOOST_CHECK(contains(haplotype_unbounded, variant_1));
    BOOST_CHECK(contains(haplotype_unbounded, variant_2));
    BOOST_CHECK(contains(haplotype_unbounded, variant_3));
    BOOST_CHECK(contains(haplotype_unbounded, variant_4));
    BOOST_CHECK(contains(haplotype_unbounded, variant_5));
    BOOST_CHECK(contains(haplotype_unbounded, variant_6));
    
    BOOST_CHECK(!haplotype_unbounded.contains(variant_1.get_reference_allele()));
    BOOST_CHECK(!haplotype_unbounded.contains(variant_2.get_reference_allele()));
    BOOST_CHECK(!haplotype_unbounded.contains(variant_3.get_reference_allele()));
    BOOST_CHECK(!haplotype_unbounded.contains(variant_4.get_reference_allele()));
    BOOST_CHECK(!haplotype_unbounded.contains(variant_5.get_reference_allele()));
    BOOST_CHECK(!haplotype_unbounded.contains(variant_6.get_reference_allele()));
    
    BOOST_CHECK(!contains(haplotype_unbounded, false_variant_1));
    BOOST_CHECK(!contains(haplotype_unbounded, false_variant_2));
    BOOST_CHECK(!contains(haplotype_unbounded, false_variant_3));
    BOOST_CHECK(!contains(haplotype_unbounded, false_variant_4));
    
    BOOST_CHECK(haplotype_unbounded.contains(get_reference_allele(a_reference_part1, human)));
    BOOST_CHECK(haplotype_unbounded.contains(get_reference_allele(a_reference_part2, human)));
    BOOST_CHECK(haplotype_unbounded.contains(get_reference_allele(a_reference_part3, human)));
    
    Haplotype haplotype_bounded {human, region};
    
    BOOST_CHECK(haplotype_bounded.get_sequence() == human.get_sequence(region));
    
    add_to_back(variant_1, haplotype_bounded);
    add_to_back(variant_2, haplotype_bounded);
    add_to_back(variant_3, haplotype_bounded);
    add_to_back(variant_4, haplotype_bounded);
    add_to_back(variant_5, haplotype_bounded);
    add_to_back(variant_6, haplotype_bounded);
    
    BOOST_CHECK(haplotype_bounded.get_sequence() == "CCTAAGCCAGGGGCGTGAC");
    
    BOOST_CHECK(contains(haplotype_bounded, variant_1));
    BOOST_CHECK(contains(haplotype_bounded, variant_2));
    BOOST_CHECK(contains(haplotype_bounded, variant_3));
    BOOST_CHECK(contains(haplotype_bounded, variant_4));
    BOOST_CHECK(contains(haplotype_bounded, variant_5));
    BOOST_CHECK(contains(haplotype_bounded, variant_6));
    
    BOOST_CHECK(!haplotype_bounded.contains(variant_1.get_reference_allele()));
    BOOST_CHECK(!haplotype_bounded.contains(variant_2.get_reference_allele()));
    BOOST_CHECK(!haplotype_bounded.contains(variant_3.get_reference_allele()));
    BOOST_CHECK(!haplotype_bounded.contains(variant_4.get_reference_allele()));
    BOOST_CHECK(!haplotype_bounded.contains(variant_5.get_reference_allele()));
    BOOST_CHECK(!haplotype_bounded.contains(variant_6.get_reference_allele()));
    
    BOOST_CHECK(!contains(haplotype_bounded, false_variant_1));
    BOOST_CHECK(!contains(haplotype_bounded, false_variant_2));
    BOOST_CHECK(!contains(haplotype_bounded, false_variant_3));
    BOOST_CHECK(!contains(haplotype_bounded, false_variant_4));
    
    GenomicRegion reference_begin_bit {"3", 1000000, 1000002};
    GenomicRegion reference_end_bit {"3", 1000019, 1000020};
    
    BOOST_CHECK(haplotype_bounded.contains(get_reference_allele(reference_begin_bit, human)));
    BOOST_CHECK(haplotype_bounded.contains(get_reference_allele(a_reference_part1, human)));
    BOOST_CHECK(haplotype_bounded.contains(get_reference_allele(a_reference_part2, human)));
    BOOST_CHECK(haplotype_bounded.contains(get_reference_allele(a_reference_part3, human)));
    BOOST_CHECK(haplotype_bounded.contains(get_reference_allele(reference_end_bit, human)));
}

//BOOST_AUTO_TEST_CASE(if_a_haplotype_contains_an_insertion_then_it_should_also_contain_the_parsimonious_version)
//{
//    auto human = make_reference(human_reference_fasta);
//    
//    auto region = parse_region("3:1000000-1000020", human);
//    
//    Allele a1 {parse_region("3:1000005-1000005", human), "ACGT"};
//    Allele a2 {parse_region("3:1000004-1000005", human), "CACGT"}; // parsimonious version
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
    auto human = make_reference(human_reference_fasta);
    
    auto region = parse_region("3:1000000-1000020", human);
    
    Allele an_allele {parse_region("3:1000010-1000012", human), "GG"};
    
    Allele a_sub_allele {parse_region("3:1000010-1000011", human), "G"};
    Allele another_sub_allele {parse_region("3:1000011-1000012", human), "G"};
    Allele not_a_sub_allele {parse_region("3:1000010-1000011", human), "C"};
    
    Haplotype hap {human, region};
    hap.push_back(an_allele);
    
    BOOST_CHECK(hap.contains(an_allele));
    BOOST_CHECK(hap.contains(a_sub_allele));
    BOOST_CHECK(hap.contains(another_sub_allele));
    BOOST_CHECK(!hap.contains(not_a_sub_allele));
}

BOOST_AUTO_TEST_CASE(deletions_decompose)
{
    auto human = make_reference(human_reference_fasta);
    
    auto region = parse_region("3:1000000-1000020", human);
    
    Allele an_allele {parse_region("3:1000010-1000012", human), ""};
    
    Allele a_sub_allele {parse_region("3:1000010-1000011", human), ""};
    Allele another_sub_allele {parse_region("3:1000011-1000012", human), ""};
    Allele not_a_sub_allele1 {parse_region("3:1000010-1000011", human), "C"};
    Allele not_a_sub_allele2 {parse_region("3:1000010-1000013", human), ""};
    
    Haplotype hap {human, region};
    hap.push_back(an_allele);
    
    BOOST_CHECK(hap.contains(an_allele));
    BOOST_CHECK(hap.contains(a_sub_allele));
    BOOST_CHECK(hap.contains(another_sub_allele));
    BOOST_CHECK(!hap.contains(not_a_sub_allele1));
    BOOST_CHECK(!hap.contains(not_a_sub_allele2));
}

BOOST_AUTO_TEST_CASE(insertions_decompose)
{
    auto human = make_reference(human_reference_fasta);
    
    auto region = parse_region("3:1000000-1000020", human);
    
    Allele an_allele {parse_region("3:1000010-1000010", human), "AT"};
    
    Allele a_sub_allele {parse_region("3:1000010-1000010", human), "A"};
    Allele another_sub_allele {parse_region("3:1000010-1000010", human), "T"};
    Allele not_a_sub_allele1 {parse_region("3:1000010-1000010", human), "C"};
    Allele not_a_sub_allele2 {parse_region("3:1000010-1000011", human), "A"};
    
    Haplotype hap {human, region};
    hap.push_back(an_allele);
    
    BOOST_CHECK(hap.contains(an_allele));
    BOOST_CHECK(hap.contains(a_sub_allele));
    BOOST_CHECK(hap.contains(another_sub_allele));
    BOOST_CHECK(!hap.contains(not_a_sub_allele1));
    BOOST_CHECK(!hap.contains(not_a_sub_allele2));
}

BOOST_AUTO_TEST_CASE(haplotype_equate_when_alleles_infer_same_sequence)
{
    auto human = make_reference(human_reference_fasta);
    
    auto region = parse_region("16:9300000-9300100", human);
    
    Allele allele1 {parse_region("16:9300037-9300037", human), "TG"};
    Allele allele2 {parse_region("16:9300039-9300051", human), ""};
    Allele allele3 {parse_region("16:9300041-9300051", human), ""};
    
    Haplotype hap1 {human, region};
    hap1.push_back(allele3);
    
    Haplotype hap2 {human, region};
    hap2.push_back(allele1);
    hap2.push_back(allele2);
    
    BOOST_CHECK(hap1.get_sequence() == hap2.get_sequence());
    BOOST_CHECK(hap1 == hap2);
    
    Allele allele4 {parse_region("16:9300037-9300038", human), "T"};
    Allele allele5 {parse_region("16:9300038-9300039", human), "C"};
    Allele allele6 {parse_region("16:9300037-9300039", human), "TC"};
    
    Haplotype hap3 {human, region};
    hap3.push_back(allele4);
    hap3.push_back(allele5);
    
    Haplotype hap4 {human, region};
    hap4.push_back(allele6);
    
    BOOST_CHECK(hap3.get_sequence() == hap4.get_sequence());
    BOOST_CHECK(hap3 == hap4);
}

BOOST_AUTO_TEST_CASE(haplotypes_can_be_compared_for_structural_complexity)
{
    auto human = make_reference(human_reference_fasta);
    
    auto region = parse_region("16:9300000-9300100", human);
    
    Allele allele1 {parse_region("16:9300037-9300037", human), "TG"};
    Allele allele2 {parse_region("16:9300039-9300051", human), ""};
    Allele allele3 {parse_region("16:9300041-9300051", human), ""};
    
    Haplotype hap1 {human, region};
    hap1.push_back(allele3);
    
    Haplotype hap2 {human, region};
    hap2.push_back(allele1);
    hap2.push_back(allele2);
    
    BOOST_CHECK(IsLessComplex()(hap1, hap2));
}

BOOST_AUTO_TEST_CASE(haplotypes_behave_at_boundries)
{
    auto human = make_reference(human_reference_fasta);
    
    auto region = parse_region("16:9299940-9300100", human);
    
    Allele allele1 {parse_region("16:9299945-9299946", human), "T"};
    Allele allele2 {parse_region("16:9299946-9299957", human), "CGCATTACAAC"};
    Allele allele3 {parse_region("16:9299957-9299958", human), "C"};
    Allele allele4 {get_reference_allele(parse_region("16:9299958-9300037", human), human)};
    Allele allele5 {parse_region("16:9300037-9300037", human), ""};
    Allele allele6 {parse_region("16:9300037-9300039", human), "TG"};
    Allele allele7 {parse_region("16:9300039-9300051", human), "TGTGTGTGCGTT"};
    Allele allele8 {parse_region("16:9300051-9300061", human), "TGTGTGTGTG"};
    Allele allele9 {parse_region("16:9300061-9300062", human), "G"};
    Allele allele10 {parse_region("16:9300062-9300072", human), "GTGTGTGTGT"};
    Allele allele11 {parse_region("16:9300072-9300073", human), "G"};
    Allele allele12 {parse_region("16:9300073-9300074", human), "G"};
    Allele allele13 {parse_region("16:9300074-9300075", human), "G"};
    
    Haplotype haplotype {human, region};
    
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
    
    Allele test_allele1 {parse_region("16:9300037-9300050", human), ""};
    Allele test_allele2 {parse_region("16:9300037-9300051", human), ""};
    Allele test_allele3 {parse_region("16:9300037-9300052", human), ""};
    
    BOOST_CHECK(!haplotype.contains(test_allele1));
    BOOST_CHECK(!haplotype.contains(test_allele2));
    BOOST_CHECK(!haplotype.contains(test_allele3));
}

BOOST_AUTO_TEST_CASE(haplotypes_can_be_copied_and_moved)
{
    auto human = make_reference(human_reference_fasta);
    
    auto region = parse_region("16:9299940-9300100", human);
    
    Allele allele1 {parse_region("16:9300037-9300037", human), "TG"};
    Allele allele2 {parse_region("16:9300039-9300051", human), ""};
    
    Haplotype hap {human, region};
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
    auto human = make_reference(human_reference_fasta);
    
    auto region = parse_region("16:9299940-9300100", human);
    
    Allele allele1 {parse_region("16:9299945-9299946", human), "T"};
    Allele allele2 {parse_region("16:9299946-9299957", human), "CGCATTACAAC"};
    Allele allele3 {parse_region("16:9299957-9299958", human), "C"};
    Allele allele4 {get_reference_allele(parse_region("16:9299958-9300037", human), human)};
    Allele allele5 {parse_region("16:9300037-9300037", human), ""};
    Allele allele6 {parse_region("16:9300037-9300039", human), "TG"};
    Allele allele7 {parse_region("16:9300039-9300051", human), "TGTGTGTGCGTT"};
    Allele allele8 {parse_region("16:9300051-9300061", human), "TGTGTGTGTG"};
    Allele allele9 {parse_region("16:9300061-9300062", human), "G"};
    Allele allele10 {parse_region("16:9300062-9300072", human), "GTGTGTGTGT"};
    Allele allele11 {parse_region("16:9300072-9300073", human), "G"};
    Allele allele12 {parse_region("16:9300073-9300074", human), "G"};
    Allele allele13 {parse_region("16:9300074-9300075", human), "G"};
    
    Haplotype haplotype {human, region};
    
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
    
    auto splice_region = parse_region("16:9299964-9300083", human);
    
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
//    auto region = parse_region("16:9300000-9300100", human);
//    
//    Allele allele1 {parse_region("16:9300037-9300037", human), "TG"};
//    Allele allele2 {parse_region("16:9300037-9300051", human), ""};
//    Allele allele3 {parse_region("16:9300051-9300051", human), ""};
//    Allele allele4 {parse_region("16:9300051-9300061", human), "TGTGTGTGTG"};
//    Allele allele5 {parse_region("16:9300061-9300062", human), "T"};
//    Allele allele6 {parse_region("16:9300062-9300072", human), "GTGTGTGTGT"};
//    Allele allele7 {parse_region("16:9300072-9300073", human), "T"};
//    Allele allele8 {parse_region("16:9300073-9300074", human), "G"};
//    Allele allele9 {parse_region("16:9300074-9300075", human), "G"};
//    
//    Allele allele10 {parse_region("16:9300037-9300039", human), "TG"};
//    Allele allele11 {parse_region("16:9300039-9300051", human), ""};
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
