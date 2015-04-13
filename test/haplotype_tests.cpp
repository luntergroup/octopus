//
//  haplotype_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 02/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>
#include <cstddef>
#include <set>

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "read_manager.h"
#include "variant.h"
#include "variant_utils.h"
#include "variant_candidate_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"

TEST_CASE("test_make_haplotypes", "[haplotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    auto a_region = parse_region("3:1000000-1000010", human);
    
    auto the_reference_sequence = human.get_sequence(a_region); // CCAACAAGCA
    
    Haplotype reference_haplotype {human, a_region};
    
    REQUIRE(reference_haplotype.get_sequence() == the_reference_sequence);
    
    Variant variant_1 {"3", 1000004, "C", "A", 0, 0};
    
    Haplotype haplotype2 {human};
    add_to_back(variant_1, haplotype2);
    
    REQUIRE(haplotype2.get_sequence(a_region) == "CCAAAAAGCA");
    
    Haplotype haplotype3 {human, a_region};
    add_to_back(variant_1, haplotype3);
    
    REQUIRE(haplotype2.get_sequence(a_region) == haplotype3.get_sequence());
    
    Variant variant_2 {"3", 1000004, "CA", "", 0, 0};
    Variant variant_3 {"3", 1000008, "", "C", 0, 0};
    
    Haplotype haplotype4 {human, a_region};
    add_to_back(variant_2, haplotype4);
    add_to_back(variant_3, haplotype4);
    
    REQUIRE(haplotype4.get_sequence() == "CCAAAGCCA");
    
    Haplotype haplotype5 {human};
    add_to_back(variant_2, haplotype5);
    add_to_back(variant_3, haplotype5);
    
    REQUIRE(haplotype5.get_sequence() == "AGC");
    
    Variant variant_4 {"3", 1000004, "CA", "GG", 0, 0};
    
    Haplotype haplotype6 {human, a_region};
    add_to_back(variant_4, haplotype6);
    
    REQUIRE(haplotype6.get_sequence() == "CCAAGGAGCA");
    
    Variant variant_5 {"3", 1000004, "C", "G", 0, 0};
    Variant variant_6 {"3", 1000005, "A", "G", 0, 0};
    
    Haplotype haplotype7 {human, a_region};
    add_to_back(variant_6, haplotype7);
    add_to_front(variant_5, haplotype7);
    
    REQUIRE(haplotype7.get_sequence() == haplotype6.get_sequence());
}

TEST_CASE("test_make_haplotype_from_candidates", "[haplotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome ecoli(a_factory.make(ecoli_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {ecoli_bam});
    
    VariantCandidateGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(ecoli, 0));
    
    auto a_region = parse_region("R00000042:99640-99745", ecoli);
    
    auto reference_sequence = ecoli.get_sequence(a_region);
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(some_reads.begin(), some_reads.end());
    
    auto variants = candidate_generator.get_candidates(a_region);
    
    REQUIRE(variants.size() == 12);
    
    Haplotype haplotype1 {ecoli, a_region};
    for (const auto& variant : variants) {
        if (is_snp(variant)) {
            add_to_back(variant, haplotype1);
        }
    }
    
    REQUIRE(haplotype1.get_sequence() == "AGCGTGGGTAAACAAAGCCATGCTATCAGCACCGCCAGCGGCGTTGGCGAACA"
            "TTTTGCTGATAAAACTGCGTTAATTACGCGTCTTAAATTACTGATTGCTGAG");
}

TEST_CASE("haplotype_reference_haplotype_test", "[haplotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    GenomicRegion a_region {"7", 1000000, 1000100};
    
    Haplotype a_reference_haplotype {human, a_region};
    
    REQUIRE(a_reference_haplotype.contains(a_region, human.get_sequence(a_region)));
    
    GenomicRegion a_sub_region {"7", 1000010, 1000090};
    
    REQUIRE(a_reference_haplotype.contains(a_sub_region, human.get_sequence(a_sub_region)));
    
    GenomicRegion a_left_overlapping_region {"7", 999999, 1000090};
    
    REQUIRE(!a_reference_haplotype.contains(a_left_overlapping_region, human.get_sequence(a_left_overlapping_region)));
    
    GenomicRegion a_right_overlapping_region {"7", 1000090, 1000101};
    
    REQUIRE(!a_reference_haplotype.contains(a_right_overlapping_region, human.get_sequence(a_right_overlapping_region)));
}

// This is absolutely crucial to make the allele posteriors correct
TEST_CASE("haplotype_variant_containment_test", "[haplotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    GenomicRegion a_region {"3", 1000000, 1000020}; // CCAACAAGCATTGGTGTGGC
    
    // Variants the haplotype will contain
    Variant variant_1 {"3", 1000002, "A", "T", 0, 0};
    Variant variant_2 {"3", 1000004, "CA", "", 0, 0};
    Variant variant_3 {"3", 1000008, "", "C", 0, 0};
    Variant variant_4 {"3", 1000010, "TT", "GG", 0, 0};
    Variant variant_5 {"3", 1000014, "T", "C", 0, 0};
    Variant variant_6 {"3", 1000018, "G", "A", 0, 0};
    
    // Parts of the haplotype which remain reference
    GenomicRegion a_reference_part1 {"3", 1000003, 1000004};
    GenomicRegion a_reference_part2 {"3", 1000012, 1000014};
    GenomicRegion a_reference_part3 {"3", 1000015, 1000018};
    
    // Variants which the haplotype does not contain
    Variant false_variant_1 {"3", 1000002, "A", "C", 0, 0};
    Variant false_variant_2 {"3", 1000008, "", "T", 0, 0};
    Variant false_variant_3 {"3", 1000010, "TT", "AC", 0, 0};
    Variant false_variant_4 {"3", 1000014, "T", "A", 0, 0};
    
    Haplotype haplotype_unbounded {human};
    add_to_back(variant_1, haplotype_unbounded);
    add_to_back(variant_2, haplotype_unbounded);
    add_to_back(variant_3, haplotype_unbounded);
    add_to_back(variant_4, haplotype_unbounded);
    add_to_back(variant_5, haplotype_unbounded);
    add_to_back(variant_6, haplotype_unbounded);
    
    REQUIRE(haplotype_unbounded.get_sequence() == "TAAGCCAGGGGCGTGA");
    
    REQUIRE(contains(haplotype_unbounded, variant_1));
    REQUIRE(contains(haplotype_unbounded, variant_2));
    REQUIRE(contains(haplotype_unbounded, variant_3));
    REQUIRE(contains(haplotype_unbounded, variant_4));
    REQUIRE(contains(haplotype_unbounded, variant_5));
    REQUIRE(contains(haplotype_unbounded, variant_6));
    
    REQUIRE(!haplotype_unbounded.contains(variant_1.get_region(), variant_1.get_reference_allele_sequence()));
    REQUIRE(!haplotype_unbounded.contains(variant_2.get_region(), variant_2.get_reference_allele_sequence()));
    REQUIRE(!haplotype_unbounded.contains(variant_3.get_region(), variant_3.get_reference_allele_sequence()));
    REQUIRE(!haplotype_unbounded.contains(variant_4.get_region(), variant_4.get_reference_allele_sequence()));
    REQUIRE(!haplotype_unbounded.contains(variant_5.get_region(), variant_5.get_reference_allele_sequence()));
    REQUIRE(!haplotype_unbounded.contains(variant_6.get_region(), variant_6.get_reference_allele_sequence()));
    
    REQUIRE(!contains(haplotype_unbounded, false_variant_1));
    REQUIRE(!contains(haplotype_unbounded, false_variant_2));
    REQUIRE(!contains(haplotype_unbounded, false_variant_3));
    REQUIRE(!contains(haplotype_unbounded, false_variant_4));
    
    REQUIRE(haplotype_unbounded.contains(a_reference_part1, human.get_sequence(a_reference_part1)));
    REQUIRE(haplotype_unbounded.contains(a_reference_part2, human.get_sequence(a_reference_part2)));
    REQUIRE(haplotype_unbounded.contains(a_reference_part3, human.get_sequence(a_reference_part3)));
    
    Haplotype haplotype_bounded {human, a_region};
    
    REQUIRE(haplotype_bounded.get_sequence() == human.get_sequence(a_region));
    
    add_to_back(variant_1, haplotype_bounded);
    add_to_back(variant_2, haplotype_bounded);
    add_to_back(variant_3, haplotype_bounded);
    add_to_back(variant_4, haplotype_bounded);
    add_to_back(variant_5, haplotype_bounded);
    add_to_back(variant_6, haplotype_bounded);
    
    REQUIRE(haplotype_bounded.get_sequence() == "CCTAAGCCAGGGGCGTGAC");
    
    REQUIRE(contains(haplotype_bounded, variant_1));
    REQUIRE(contains(haplotype_bounded, variant_2));
    REQUIRE(contains(haplotype_bounded, variant_3));
    REQUIRE(contains(haplotype_bounded, variant_4));
    REQUIRE(contains(haplotype_bounded, variant_5));
    REQUIRE(contains(haplotype_bounded, variant_6));
    
    REQUIRE(!haplotype_bounded.contains(variant_1.get_region(), variant_1.get_reference_allele_sequence()));
    REQUIRE(!haplotype_bounded.contains(variant_2.get_region(), variant_2.get_reference_allele_sequence()));
    REQUIRE(!haplotype_bounded.contains(variant_3.get_region(), variant_3.get_reference_allele_sequence()));
    REQUIRE(!haplotype_bounded.contains(variant_4.get_region(), variant_4.get_reference_allele_sequence()));
    REQUIRE(!haplotype_bounded.contains(variant_5.get_region(), variant_5.get_reference_allele_sequence()));
    REQUIRE(!haplotype_bounded.contains(variant_6.get_region(), variant_6.get_reference_allele_sequence()));
    
    REQUIRE(!contains(haplotype_bounded, false_variant_1));
    REQUIRE(!contains(haplotype_bounded, false_variant_2));
    REQUIRE(!contains(haplotype_bounded, false_variant_3));
    REQUIRE(!contains(haplotype_bounded, false_variant_4));
    
    GenomicRegion reference_begin_bit {"3", 1000000, 1000002};
    GenomicRegion reference_end_bit {"3", 1000019, 1000020};
    
    REQUIRE(haplotype_bounded.contains(reference_begin_bit, human.get_sequence(reference_begin_bit)));
    REQUIRE(haplotype_bounded.contains(a_reference_part1, human.get_sequence(a_reference_part1)));
    REQUIRE(haplotype_bounded.contains(a_reference_part2, human.get_sequence(a_reference_part2)));
    REQUIRE(haplotype_bounded.contains(a_reference_part3, human.get_sequence(a_reference_part3)));
    REQUIRE(haplotype_bounded.contains(reference_end_bit, human.get_sequence(reference_end_bit)));
}

TEST_CASE("haplotype_equality_test", "[haplotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    auto a_region = parse_region("16:9300000-9300100", human);
    
    Allele allele1 {parse_region("16:9300037-9300037", human), "TG"};
    Allele allele2 {parse_region("16:9300037-9300051", human), ""};
    Allele allele3 {parse_region("16:9300039-9300051", human), ""};
    
    Haplotype hap1 {human, a_region};
    hap1.push_back(allele3);
    
    Haplotype hap2 {human, a_region};
    hap2.push_back(allele1);
    hap2.push_back(allele2);
    
    std::cout << hap1 << std::endl;
    std::cout << hap2 << std::endl;
    
    //std::cout << hap2.contains(allele3) << std::endl; // CAUSES EXCEPTION!
    
    //REQUIRE(hap1 == hap2);
}
