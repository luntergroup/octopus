//
//  variant_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include <iostream>
#include <vector>
#include <algorithm>

#include "catch.hpp"

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "genomic_region.h"
#include "region_utils.h"
#include "allele.h"
#include "variant.h"
#include "variant_utils.h"
#include "mock_objects.h"

TEST_CASE("< is consistent with ==", "[variant]")
{
    Variant snp1 {"chr1", 100, "C", "A", 0, 0};
    Variant snp2 {"chr1", 99, "C", "A", 0, 0};
    Variant snp3 {"chr1", 100, "C", "T", 0, 0};
    Variant snp4 {"chr1", 100, "C", "A", 0, 0};
    
    bool r1 = !(snp1 < snp2) && !(snp2 < snp1);
    bool r2 = (snp1 == snp2);
    REQUIRE(!r2);
    REQUIRE(r1 == r2);
    
    r1 = !(snp1 < snp3) && !(snp3 < snp1);
    r2 = snp1 == snp3;
    REQUIRE(!r2);
    REQUIRE(r1 == r2);
    
    r1 = !(snp1 < snp4) && !(snp4 < snp1);
    r2 = snp1 == snp4;
    REQUIRE(r2);
    REQUIRE(r1 == r2);
    
    Variant del1 {"chr1", 100, "C", "", 0, 0};
    Variant del2 {"chr1", 100, "CA", "", 0, 0};
    Variant del3 {"chr1", 100, "C", "", 0, 0};
    
    r1 = !(del1 < del2) && !(del2 < del1);
    r2 = del1 == del2;
    REQUIRE(!r2);
    REQUIRE(r1 == r2);
    
    r1 = !(del1 < del3) && !(del3 < del1);
    r2 = del1 == del3;
    REQUIRE(r2);
    REQUIRE(r1 == r2);
    
    Variant ins1 {"chr1", 100, "", "C", 0, 0};
    Variant ins2 {"chr1", 100, "", "CA", 0, 0};
    Variant ins3 {"chr1", 100, "", "C", 0, 0};
    
    r1 = !(ins1 < ins2) && !(ins2 < ins1);
    r2 = ins1 == ins2;
    REQUIRE(!r2);
    REQUIRE(r1 == r2);
    
    r1 = !(ins1 < ins3) && !(ins3 < ins1);
    r2 = ins1 == ins3;
    REQUIRE(r2);
    REQUIRE(r1 == r2);
}

TEST_CASE("can binary search for ordered variants with < operator", "[variant]")
{
    Variant snp1 {"chr1", 100, "C", "A", 0, 0};
    Variant snp2 {"chr1", 105, "G", "T", 0, 0};
    Variant del1 {"chr1", 103, "C", "", 0, 0};
    Variant del2 {"chr1", 115, "T", "", 0, 0};
    Variant ins1 {"chr1", 107, "", "G", 0, 0};
    Variant ins2 {"chr1", 110, "", "CA", 0, 0};
    
    std::vector<Variant> variants {snp1, del1, snp2, ins1, ins2, del2};
    
    REQUIRE(std::binary_search(variants.cbegin(), variants.cend(), snp1));
    REQUIRE(std::binary_search(variants.cbegin(), variants.cend(), snp2));
    REQUIRE(std::binary_search(variants.cbegin(), variants.cend(), del1));
    REQUIRE(std::binary_search(variants.cbegin(), variants.cend(), del2));
    REQUIRE(std::binary_search(variants.cbegin(), variants.cend(), ins1));
    REQUIRE(std::binary_search(variants.cbegin(), variants.cend(), ins2));
    
    Variant non_snp1 {"chr1", 111, "C", "A", 0, 0};
    Variant non_snp2 {"chr1", 105, "G", "A", 0, 0};
    Variant non_del1 {"chr1", 104, "C", "", 0, 0};
    Variant non_del2 {"chr1", 115, "TA", "", 0, 0};
    Variant non_ins1 {"chr1", 108, "", "G", 0, 0};
    Variant non_ins2 {"chr1", 110, "", "TA", 0, 0};
    
    REQUIRE(!std::binary_search(variants.cbegin(), variants.cend(), non_snp1));
    REQUIRE(!std::binary_search(variants.cbegin(), variants.cend(), non_snp2));
    REQUIRE(!std::binary_search(variants.cbegin(), variants.cend(), non_del1));
    REQUIRE(!std::binary_search(variants.cbegin(), variants.cend(), non_del2));
    REQUIRE(!std::binary_search(variants.cbegin(), variants.cend(), non_ins1));
    REQUIRE(!std::binary_search(variants.cbegin(), variants.cend(), non_ins2));
    
    Variant del3 {"chr1", 101, "C", "", 0, 0};
    Variant snp3 {"chr1", 101, "C", "", 0, 0};
    
    std::vector<Variant> variants2 {snp1, del3, snp3};
    
    REQUIRE(std::binary_search(variants2.cbegin(), variants2.cend(), snp1));
    REQUIRE(std::binary_search(variants2.cbegin(), variants2.cend(), del3));
    REQUIRE(std::binary_search(variants2.cbegin(), variants2.cend(), snp3));
}

TEST_CASE("snps do not overlap adjacent snps", "[snps]")
{
    Variant snp1 {"chr1", 100, "C", "A", 0, 0};
    Variant snp2 {"chr1", 99, "C", "A", 0, 0};
    Variant snp3 {"chr1", 101, "C", "A", 0, 0};
    Variant snp4 {"chr1", 100, "C", "T", 0, 0};
    Variant snp5 {"chr1", 99, "C", "T", 0, 0};
    Variant snp6 {"chr1", 101, "C", "T", 0, 0};
    
    REQUIRE(overlaps(snp1, snp1));
    REQUIRE(overlaps(snp1, snp4));
    REQUIRE(!overlaps(snp1, snp2));
    REQUIRE(!overlaps(snp1, snp3));
    REQUIRE(!overlaps(snp1, snp5));
    REQUIRE(!overlaps(snp1, snp6));
    
    REQUIRE(are_adjacent(snp1, snp2));
    REQUIRE(!are_adjacent(snp2, snp3));
}

TEST_CASE("mnps overlap correctly", "[mnp]")
{
    Variant mnp1 {"chr1", 100, "CAT", "TAC", 0, 0};
    Variant mnp2 {"chr1", 99, "CAT", "TAC", 0, 0};
    Variant mnp3 {"chr1", 101, "CAT", "TAC", 0, 0};
    Variant mnp4 {"chr1", 100, "CAT", "TAG", 0, 0};
    Variant mnp5 {"chr1", 99, "CAT", "TAG", 0, 0};
    Variant mnp6 {"chr1", 101, "CAT", "TAG", 0, 0};
    
    // edge cases
    Variant mnp7  {"chr1", 98, "CAT", "TAC", 0, 0};
    Variant mnp8  {"chr1", 102, "CAT", "TAC", 0, 0};
    Variant mnp9  {"chr1", 97, "CAT", "TAC", 0, 0};
    Variant mnp10 {"chr1", 103, "CAT", "TAC", 0, 0};
    
    REQUIRE(overlaps(mnp1, mnp1));
    REQUIRE(overlaps(mnp1, mnp2));
    REQUIRE(overlaps(mnp1, mnp3));
    REQUIRE(overlaps(mnp1, mnp4));
    REQUIRE(overlaps(mnp1, mnp5));
    REQUIRE(overlaps(mnp1, mnp6));
    
    REQUIRE(overlaps(mnp1, mnp7));
    REQUIRE(overlaps(mnp1, mnp8));
    REQUIRE(!overlaps(mnp1, mnp9));
    REQUIRE(!overlaps(mnp1, mnp10));
}

TEST_CASE("insertions overlap other insertions with same region", "[insertion]")
{
    Variant insert1 {"chr1", 100, "", "TAG", 0, 0};
    Variant insert2 {"chr1", 99, "", "TAG", 0, 0};
    Variant insert3 {"chr1", 101, "", "TAG", 0, 0};
    
    REQUIRE(overlaps(insert1, insert1));
    REQUIRE(!overlaps(insert1, insert2));
    REQUIRE(!overlaps(insert1, insert3));
}

TEST_CASE("insertions overlap with other variants when contained by their region", "[insertion]")
{
    Variant insert {"chr1", 100, "", "TAG", 0, 0};
    Variant del {"chr1", 99, "TAG", "", 0, 0};
    Variant mnp {"chr1", 100, "TAG", "CAT", 0, 0};
    
    REQUIRE(overlaps(insert, del));
    REQUIRE(overlaps(insert, mnp));
}

TEST_CASE("deletions overlap in the same way as mnps", "[deletion]")
{
    Variant del1 {"chr1", 100, "TAG", "", 0, 0};
    Variant del2 {"chr1", 99, "TAG", "", 0, 0};
    Variant del3 {"chr1", 101, "TAG", "", 0, 0};
    
    REQUIRE(overlaps(del1, del1));
    REQUIRE(overlaps(del1, del2));
    REQUIRE(overlaps(del1, del3));
}

TEST_CASE("variants are ordered by region and lexicographically by sequence", "[variant]")
{
    Variant snp1 {"chr1", 100, "T", "A", 0, 0};
    Variant snp2 {"chr1", 100, "T", "C", 0, 0};
    Variant snp3 {"chr1", 100, "T", "G", 0, 0};
    Variant ins1 {"chr1", 100, "", "AG", 0, 0};
    Variant ins2 {"chr1", 100, "", "CC", 0, 0};
    Variant ins3 {"chr1", 100, "", "CCA", 0, 0};
    Variant del1 {"chr1", 100, "TA", "", 0, 0};
    Variant del2 {"chr1", 100, "TAG", "", 0, 0};
    
    std::vector<Variant> variants1 {snp1, ins1, snp2};
    
    std::sort(variants1.begin(), variants1.end());
    
    std::vector<Variant> variants1_required_sort {ins1, snp1, snp2};
    
    bool is_required_sort1 = std::equal(variants1.cbegin(), variants1.cend(), variants1_required_sort.cbegin());
    
    REQUIRE(is_required_sort1);
    
    std::vector<Variant> variants2 {snp1, del1, snp2, ins1, snp3};
    
    std::sort(variants2.begin(), variants2.end());
    
    std::vector<Variant> variants2_required_sort {ins1, snp1, snp2, snp3, del1};
    
    auto is_required_sort2 = std::equal(variants2.cbegin(), variants2.cend(), variants2_required_sort.cbegin());
    
    REQUIRE(is_required_sort2);
    
    std::vector<Variant> variants3 {del2, snp1, del1, snp2, ins1, snp3, ins2};
    
    std::sort(variants3.begin(), variants3.end());
    
    std::vector<Variant> variants3_required_sort {ins1, ins2, snp1, snp2, snp3, del1, del2};
    
    auto is_required_sort3 = std::equal(variants3.cbegin(), variants3.cend(), variants3_required_sort.cbegin());
    
    REQUIRE(is_required_sort3);
    
    std::vector<Variant> variants4 {ins2, del2, snp1, del1, snp2, ins1, snp3, ins3};
    
    std::sort(variants4.begin(), variants4.end());
    
    std::vector<Variant> variants4_required_sort {ins1, ins2, ins3, snp1, snp2, snp3, del1, del2};
    
    auto is_required_sort4 = std::equal(variants4.cbegin(), variants4.cend(), variants4_required_sort.cbegin());
    
    REQUIRE(is_required_sort4);
}

TEST_CASE("overlap_range includes insertions on boundry", "[variant]")
{
    Variant snp1 {"chr1", 100, "T", "A", 0, 0};
    Variant snp2 {"chr1", 110, "T", "C", 0, 0};
    Variant ins1 {"chr1", 105, "", "A", 0, 0};
    
    std::vector<Variant> variants {snp1, ins1, snp2};
    
    GenomicRegion region1 {"chr1", 104, 106};
    GenomicRegion region2 {"chr1", 104, 105};
    GenomicRegion region3 {"chr1", 105, 106};
    GenomicRegion region4 {"chr1", 105, 105};
    
    auto overlapped1 = overlap_range(variants.cbegin(), variants.cend(), region1);
    auto overlapped2 = overlap_range(variants.cbegin(), variants.cend(), region2);
    auto overlapped3 = overlap_range(variants.cbegin(), variants.cend(), region3);
    auto overlapped4 = overlap_range(variants.cbegin(), variants.cend(), region4);
    
    REQUIRE(std::distance(overlapped1.first, overlapped1.second) == 1);
    REQUIRE(std::distance(overlapped2.first, overlapped2.second) == 0);
    REQUIRE(std::distance(overlapped3.first, overlapped3.second) == 1);
    REQUIRE(std::distance(overlapped4.first, overlapped4.second) == 1);
}

TEST_CASE("indels can be left aligned", "[left_alignment]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    // Huntingtin region CCAGCAGCAGCAGCAG...
    auto a_region = parse_region("4:3076657-3076660", human);
    
    auto the_sequence = human.get_sequence(a_region);
    
    REQUIRE(the_sequence == "CAG");
    
    Variant a_deletion {a_region, the_sequence, "", 0, 0};
    
    Variant left_aligned_deletion = left_align(a_deletion, human);
    
    REQUIRE(left_aligned_deletion.get_region() ==
            parse_region("4:3076603-3076606", human));
    REQUIRE(left_aligned_deletion.get_reference_allele_sequence() == "CAG");
    REQUIRE(left_aligned_deletion.get_alternative_allele_sequence() == "");
    
    Variant an_insertion {parse_region("4:3076660-3076660", human),
                                               "", the_sequence, 0, 0};
    
    auto left_aligned_insertion = left_align(an_insertion, human);
    
    REQUIRE(left_aligned_insertion.get_region() ==
            parse_region("4:3076603-3076603", human));
    REQUIRE(left_aligned_insertion.get_reference_allele_sequence() == "");
    REQUIRE(left_aligned_insertion.get_alternative_allele_sequence() == "CAG");
    
    // Region is CCAACAACAACAACAC (94594947-94594962)
    a_region = parse_region("5:94594956-94594959", human);
    
    the_sequence = human.get_sequence(a_region);
    
    REQUIRE(the_sequence == "CAA");
    
    a_deletion = Variant {a_region, the_sequence, "", 0, 0};
    
    left_aligned_deletion = left_align(a_deletion, human);
    
    REQUIRE(left_aligned_deletion.get_region() ==
            parse_region("5:94594949-94594952", human));
    REQUIRE(left_aligned_deletion.get_reference_allele_sequence() == "ACA");
    REQUIRE(left_aligned_deletion.get_alternative_allele_sequence() == "");
    
    an_insertion = Variant {parse_region("5:94594959-94594959", human), "", the_sequence, 0, 0};
    
    left_aligned_insertion = left_align(an_insertion, human);
    
    REQUIRE(left_aligned_insertion.get_region() ==
            parse_region("5:94594949-94594949", human));
    REQUIRE(left_aligned_insertion.get_reference_allele_sequence() == "");
    REQUIRE(left_aligned_insertion.get_alternative_allele_sequence() == "ACA");
}

TEST_CASE("can make variants parsimonious", "[parsimonious]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    Variant a_snp {parse_region("12:10001330-10001331", human), std::string {"G"}, std::string {"C"}, 0, 0};
    
    REQUIRE(is_parsimonious(a_snp));
    REQUIRE(make_parsimonious(a_snp, human) == a_snp);
    
    Variant an_unparsimonious_snp {parse_region("12:10001330-10001332", human), std::string {"GT"}, std::string {"CT"}, 0, 0};
    
    REQUIRE(!is_parsimonious(an_unparsimonious_snp));
    
    auto parsimonised_snp = make_parsimonious(an_unparsimonious_snp, human);
    
    REQUIRE(is_parsimonious(parsimonised_snp));
    REQUIRE(parsimonised_snp == a_snp);
    
    Variant another_unparsimonious_snp {parse_region("12:10001329-10001332", human), std::string {"TGT"}, std::string {"TCT"}, 0, 0};
    
    REQUIRE(!is_parsimonious(another_unparsimonious_snp));
    
    auto another_parsimonised_snp = make_parsimonious(another_unparsimonious_snp, human);
    
    REQUIRE(is_parsimonious(another_parsimonised_snp));
    REQUIRE(another_parsimonised_snp == a_snp);
    
    auto a_region = parse_region("12:10001330-10001335", human);
    
    auto the_sequence = human.get_sequence(a_region);
    
    REQUIRE(the_sequence == "GTGGA");
    
    Variant a_deletion {a_region, the_sequence, "", 0, 0};
    
    auto parsimonious_deletion = make_parsimonious(a_deletion, human);
    
    REQUIRE(parsimonious_deletion.get_region() ==
            parse_region("12:10001329-10001335", human));
    REQUIRE(parsimonious_deletion.get_reference_allele_sequence() == "CGTGGA");
    REQUIRE(parsimonious_deletion.get_alternative_allele_sequence() == "C");
    
    Variant an_insertion {parse_region("12:10001330-10001330", human), "", the_sequence, 0, 0};
    
    auto parsimonious_insertion = make_parsimonious(an_insertion, human);
    
    REQUIRE(parsimonious_insertion.get_region() ==
            parse_region("12:10001329-10001330", human));
    REQUIRE(parsimonious_insertion.get_reference_allele_sequence() == "C");
    REQUIRE(parsimonious_insertion.get_alternative_allele_sequence() == "CGTGGA");
    
    Variant an_unparsimonious_deletion {parse_region("12:10001328-10001335", human), "TCGTGGA", "TC", 0, 0};
    
    REQUIRE(!is_parsimonious(an_unparsimonious_deletion));
    
    auto parsimonised_deletion = make_parsimonious(an_unparsimonious_deletion, human);
    
    REQUIRE(is_parsimonious(parsimonised_deletion));
    
    Variant an_unparsimonious_insertion {parse_region("12:10001329-10001331", human), "CG", "CGTGGA", 0, 0};
    
    REQUIRE(!is_parsimonious(an_unparsimonious_insertion));
    
    auto parsimonised_insertion = make_parsimonious(an_unparsimonious_insertion, human);
    
    REQUIRE(is_parsimonious(parsimonised_insertion));
}

TEST_CASE("can normalise variants", "[normalisation]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    // Huntingtin region CCAGCAGCAGCAGCAG...
    
    Variant a_snp {parse_region("4:3076657-3076658", human), std::string {"G"}, std::string {"C"}, 0, 0};
    
    REQUIRE(is_parsimonious(a_snp));
    
    auto a_normalised_snp = normalise(a_snp, human);
    
    REQUIRE(a_normalised_snp == a_snp);
    
    auto a_region = parse_region("4:3076657-3076660", human);
    
    auto the_sequence = human.get_sequence(a_region);
    
    REQUIRE(the_sequence == "CAG");
    
    Variant a_mnp {a_region, the_sequence, std::string {"GAC"}, 0, 0};
    
    REQUIRE(is_parsimonious(a_mnp));
    
    auto a_normalised_mnp = normalise(a_mnp, human);
    
    REQUIRE(a_normalised_mnp == a_mnp);
    
    Variant a_deletion {a_region, the_sequence, "", 0, 0};
    
    REQUIRE(!is_parsimonious(a_deletion));
    
    auto left_aligned_unparsimonious_deletion = left_align(a_deletion, human);
    
    REQUIRE(!is_parsimonious(left_aligned_unparsimonious_deletion));
    
    auto normilised_deletion = normalise(a_deletion, human);
    
    REQUIRE(is_parsimonious(normilised_deletion));
    REQUIRE(normilised_deletion.get_region() ==
            parse_region("4:3076602-3076606", human));
    REQUIRE(normilised_deletion.get_reference_allele_sequence() == "CCAG");
    REQUIRE(normilised_deletion.get_alternative_allele_sequence() == "C");
    
    Variant an_insertion {parse_region("4:3076660-3076660", human), "", the_sequence, 0, 0};
    
    REQUIRE(!is_parsimonious(an_insertion));
    
    auto left_aligned_unparsimonious_insertion = left_align(an_insertion, human);
    
    REQUIRE(!is_parsimonious(left_aligned_unparsimonious_insertion));
    
    auto normilised_insertion = normalise(an_insertion, human);
    
    REQUIRE(is_parsimonious(normilised_insertion));
    REQUIRE(normilised_insertion.get_region() ==
            parse_region("4:3076602-3076603", human));
    REQUIRE(normilised_insertion.get_reference_allele_sequence() == "C");
    REQUIRE(normilised_insertion.get_alternative_allele_sequence() == "CCAG");
    
    // Some hard ones
    
    Variant an_unormilised_snp {parse_region("4:3076656-3076659", human), std::string {"AGC"}, std::string {"ACC"}, 0, 0};
    
    REQUIRE(!is_parsimonious(an_unormilised_snp));
    
    a_normalised_snp = normalise(an_unormilised_snp, human);
    
    REQUIRE(is_parsimonious(a_normalised_mnp));
    REQUIRE(a_normalised_snp.get_region() == parse_region("4:3076657-3076658", human));
    REQUIRE(a_normalised_snp.get_reference_allele_sequence() == "G");
    REQUIRE(a_normalised_snp.get_alternative_allele_sequence() == "C");
    
    Variant an_unormilised_mnp {parse_region("4:3076656-3076661", human), std::string {"GCAGC"}, std::string {"GGACC"}, 0, 0};
    
    REQUIRE(!is_parsimonious(an_unormilised_mnp));
    
    a_normalised_mnp = normalise(an_unormilised_mnp, human);
    
    REQUIRE(is_parsimonious(a_normalised_mnp));
    REQUIRE(a_normalised_mnp.get_region() == parse_region("4:3076657-3076660", human));
    REQUIRE(a_normalised_mnp.get_reference_allele_sequence() == "CAG");
    REQUIRE(a_normalised_mnp.get_alternative_allele_sequence() == "GAC");
    
    Variant an_unnormilised_deletion {parse_region("4:3076655-3076660", human), std::string {"AGCAG"}, std::string {"AG"}, 0, 0};
    
    REQUIRE(!is_parsimonious(an_unnormilised_deletion));
    
    auto a_normalised_deletion = normalise(an_unnormilised_deletion, human);
    
    REQUIRE(is_parsimonious(a_normalised_deletion));
    REQUIRE(a_normalised_deletion.get_region() ==
            parse_region("4:3076602-3076606", human));
    REQUIRE(a_normalised_deletion.get_reference_allele_sequence() == "CCAG");
    REQUIRE(a_normalised_deletion.get_alternative_allele_sequence() == "C");
    
    Variant an_unnormilised_insertion {parse_region("4:3076655-3076657", human),
                                    std::string {"AG"}, std::string {"AGCAG"}, 0, 0};
    
    REQUIRE(!is_parsimonious(an_unnormilised_insertion));
    
    auto a_normalised_insertion = normalise(an_unnormilised_insertion, human);
    
    REQUIRE(is_parsimonious(a_normalised_insertion));
    REQUIRE(a_normalised_insertion.get_region() ==
            parse_region("4:3076602-3076603", human));
    REQUIRE(a_normalised_insertion.get_reference_allele_sequence() == "C");
    REQUIRE(a_normalised_insertion.get_alternative_allele_sequence() == "CCAG");
}
