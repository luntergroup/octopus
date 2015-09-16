//
//  variant_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>
#include <algorithm>

#include "test_common.h"
#include "reference_genome.h"
#include "genomic_region.h"
#include "mappable_algorithms.h"
#include "allele.h"
#include "variant.h"
#include "variant_utils.h"
#include "mock_objects.h"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(operator_are_consistent)
{
    Variant snp1 {"chr1", 100, "C", "A"};
    Variant snp2 {"chr1", 99, "C", "A"};
    Variant snp3 {"chr1", 100, "C", "T"};
    Variant snp4 {"chr1", 100, "C", "A"};
    
    bool r1 = !(snp1 < snp2) && !(snp2 < snp1);
    bool r2 = (snp1 == snp2);
    BOOST_CHECK(!r2);
    BOOST_CHECK(r1 == r2);
    
    r1 = !(snp1 < snp3) && !(snp3 < snp1);
    r2 = snp1 == snp3;
    BOOST_CHECK(!r2);
    BOOST_CHECK(r1 == r2);
    
    r1 = !(snp1 < snp4) && !(snp4 < snp1);
    r2 = snp1 == snp4;
    BOOST_CHECK(r2);
    BOOST_CHECK(r1 == r2);
    
    Variant del1 {"chr1", 100, "C", ""};
    Variant del2 {"chr1", 100, "CA", ""};
    Variant del3 {"chr1", 100, "C", ""};
    
    r1 = !(del1 < del2) && !(del2 < del1);
    r2 = del1 == del2;
    BOOST_CHECK(!r2);
    BOOST_CHECK(r1 == r2);
    
    r1 = !(del1 < del3) && !(del3 < del1);
    r2 = del1 == del3;
    BOOST_CHECK(r2);
    BOOST_CHECK(r1 == r2);
    
    Variant ins1 {"chr1", 100, "", "C"};
    Variant ins2 {"chr1", 100, "", "CA"};
    Variant ins3 {"chr1", 100, "", "C"};
    
    r1 = !(ins1 < ins2) && !(ins2 < ins1);
    r2 = ins1 == ins2;
    BOOST_CHECK(!r2);
    BOOST_CHECK(r1 == r2);
    
    r1 = !(ins1 < ins3) && !(ins3 < ins1);
    r2 = ins1 == ins3;
    BOOST_CHECK(r2);
    BOOST_CHECK(r1 == r2);
}

BOOST_AUTO_TEST_CASE(can_binary_search_for_ordered_variants)
{
    Variant snp1 {"chr1", 100, "C", "A"};
    Variant snp2 {"chr1", 105, "G", "T"};
    Variant del1 {"chr1", 103, "C", ""};
    Variant del2 {"chr1", 115, "T", ""};
    Variant ins1 {"chr1", 107, "", "G"};
    Variant ins2 {"chr1", 110, "", "CA"};
    
    std::vector<Variant> variants {snp1, del1, snp2, ins1, ins2, del2};
    
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), snp1));
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), snp2));
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), del1));
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), del2));
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), ins1));
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), ins2));
    
    Variant non_snp1 {"chr1", 111, "C", "A"};
    Variant non_snp2 {"chr1", 105, "G", "A"};
    Variant non_del1 {"chr1", 104, "C", ""};
    Variant non_del2 {"chr1", 115, "TA", ""};
    Variant non_ins1 {"chr1", 108, "", "G"};
    Variant non_ins2 {"chr1", 110, "", "TA"};
    
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_snp1));
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_snp2));
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_del1));
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_del2));
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_ins1));
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_ins2));
    
    Variant del3 {"chr1", 101, "C", ""};
    Variant snp3 {"chr1", 101, "C", ""};
    
    std::vector<Variant> variants2 {snp1, del3, snp3};
    
    BOOST_CHECK(std::binary_search(variants2.cbegin(), variants2.cend(), snp1));
    BOOST_CHECK(std::binary_search(variants2.cbegin(), variants2.cend(), del3));
    BOOST_CHECK(std::binary_search(variants2.cbegin(), variants2.cend(), snp3));
}

BOOST_AUTO_TEST_CASE(snps_do_not_overlap_adjacent_snps)
{
    Variant snp1 {"chr1", 100, "C", "A"};
    Variant snp2 {"chr1", 99, "C", "A"};
    Variant snp3 {"chr1", 101, "C", "A"};
    Variant snp4 {"chr1", 100, "C", "T"};
    Variant snp5 {"chr1", 99, "C", "T"};
    Variant snp6 {"chr1", 101, "C", "T"};
    
    BOOST_CHECK(overlaps(snp1, snp1));
    BOOST_CHECK(overlaps(snp1, snp4));
    BOOST_CHECK(!overlaps(snp1, snp2));
    BOOST_CHECK(!overlaps(snp1, snp3));
    BOOST_CHECK(!overlaps(snp1, snp5));
    BOOST_CHECK(!overlaps(snp1, snp6));
    
    BOOST_CHECK(are_adjacent(snp1, snp2));
    BOOST_CHECK(!are_adjacent(snp2, snp3));
}

BOOST_AUTO_TEST_CASE(mnps_overlap_correctly)
{
    Variant mnp1 {"chr1", 100, "CAT", "TAC"};
    Variant mnp2 {"chr1", 99, "CAT", "TAC"};
    Variant mnp3 {"chr1", 101, "CAT", "TAC"};
    Variant mnp4 {"chr1", 100, "CAT", "TAG"};
    Variant mnp5 {"chr1", 99, "CAT", "TAG"};
    Variant mnp6 {"chr1", 101, "CAT", "TAG"};
    
    // edge cases
    Variant mnp7  {"chr1", 98, "CAT", "TAC"};
    Variant mnp8  {"chr1", 102, "CAT", "TAC"};
    Variant mnp9  {"chr1", 97, "CAT", "TAC"};
    Variant mnp10 {"chr1", 103, "CAT", "TAC"};
    
    BOOST_CHECK(overlaps(mnp1, mnp1));
    BOOST_CHECK(overlaps(mnp1, mnp2));
    BOOST_CHECK(overlaps(mnp1, mnp3));
    BOOST_CHECK(overlaps(mnp1, mnp4));
    BOOST_CHECK(overlaps(mnp1, mnp5));
    BOOST_CHECK(overlaps(mnp1, mnp6));
    
    BOOST_CHECK(overlaps(mnp1, mnp7));
    BOOST_CHECK(overlaps(mnp1, mnp8));
    BOOST_CHECK(!overlaps(mnp1, mnp9));
    BOOST_CHECK(!overlaps(mnp1, mnp10));
}

BOOST_AUTO_TEST_CASE(insertions_overlap_other_insertions_with_same_region)
{
    Variant insert1 {"chr1", 100, "", "TAG"};
    Variant insert2 {"chr1", 99, "", "TAG"};
    Variant insert3 {"chr1", 101, "", "TAG"};
    
    BOOST_CHECK(overlaps(insert1, insert1));
    BOOST_CHECK(!overlaps(insert1, insert2));
    BOOST_CHECK(!overlaps(insert1, insert3));
}

BOOST_AUTO_TEST_CASE(insertions_overlap_with_other_variants_when_contained_by_their_region)
{
    Variant insert {"chr1", 100, "", "TAG"};
    Variant del {"chr1", 99, "TAG", ""};
    Variant mnp {"chr1", 100, "TAG", "CAT"};
    
    BOOST_CHECK(overlaps(insert, del));
    BOOST_CHECK(overlaps(insert, mnp));
}

BOOST_AUTO_TEST_CASE(deletions_overlap_in_the_same_way_as_mnps)
{
    Variant del1 {"chr1", 100, "TAG", ""};
    Variant del2 {"chr1", 99, "TAG", ""};
    Variant del3 {"chr1", 101, "TAG", ""};
    
    BOOST_CHECK(overlaps(del1, del1));
    BOOST_CHECK(overlaps(del1, del2));
    BOOST_CHECK(overlaps(del1, del3));
}

BOOST_AUTO_TEST_CASE(variants_are_ordered_by_region_and_lexicographically_by_sequence)
{
    Variant snp1 {"chr1", 100, "T", "A"};
    Variant snp2 {"chr1", 100, "T", "C"};
    Variant snp3 {"chr1", 100, "T", "G"};
    Variant ins1 {"chr1", 100, "", "AG"};
    Variant ins2 {"chr1", 100, "", "CC"};
    Variant ins3 {"chr1", 100, "", "CCA"};
    Variant del1 {"chr1", 100, "TA", ""};
    Variant del2 {"chr1", 100, "TAG", ""};
    
    std::vector<Variant> variants1 {snp1, ins1, snp2};
    
    std::sort(variants1.begin(), variants1.end());
    
    std::vector<Variant> variants1_required_sort {ins1, snp1, snp2};
    
    bool is_required_sort1 = std::equal(variants1.cbegin(), variants1.cend(), variants1_required_sort.cbegin());
    
    BOOST_CHECK(is_required_sort1);
    
    std::vector<Variant> variants2 {snp1, del1, snp2, ins1, snp3};
    
    std::sort(variants2.begin(), variants2.end());
    
    std::vector<Variant> variants2_required_sort {ins1, snp1, snp2, snp3, del1};
    
    auto is_required_sort2 = std::equal(variants2.cbegin(), variants2.cend(), variants2_required_sort.cbegin());
    
    BOOST_CHECK(is_required_sort2);
    
    std::vector<Variant> variants3 {del2, snp1, del1, snp2, ins1, snp3, ins2};
    
    std::sort(variants3.begin(), variants3.end());
    
    std::vector<Variant> variants3_required_sort {ins1, ins2, snp1, snp2, snp3, del1, del2};
    
    auto is_required_sort3 = std::equal(variants3.cbegin(), variants3.cend(), variants3_required_sort.cbegin());
    
    BOOST_CHECK(is_required_sort3);
    
    std::vector<Variant> variants4 {ins2, del2, snp1, del1, snp2, ins1, snp3, ins3};
    
    std::sort(variants4.begin(), variants4.end());
    
    std::vector<Variant> variants4_required_sort {ins1, ins2, ins3, snp1, snp2, snp3, del1, del2};
    
    auto is_required_sort4 = std::equal(variants4.cbegin(), variants4.cend(), variants4_required_sort.cbegin());
    
    BOOST_CHECK(is_required_sort4);
}

BOOST_AUTO_TEST_CASE(overlap_range_includes_insertions_on_left_boundry_but_not_the_right)
{
    Variant snp1 {"chr1", 100, "T", "A"};
    Variant snp2 {"chr1", 110, "T", "C"};
    Variant ins1 {"chr1", 105, "", "A"};
    
    std::vector<Variant> variants {snp1, ins1, snp2};
    
    GenomicRegion region1 {"chr1", 104, 106};
    GenomicRegion region2 {"chr1", 104, 105};
    GenomicRegion region3 {"chr1", 105, 106};
    GenomicRegion region4 {"chr1", 105, 105};
    
    auto overlapped1 = overlap_range(variants.cbegin(), variants.cend(), region1);
    auto overlapped2 = overlap_range(variants.cbegin(), variants.cend(), region2);
    auto overlapped3 = overlap_range(variants.cbegin(), variants.cend(), region3);
    auto overlapped4 = overlap_range(variants.cbegin(), variants.cend(), region4);
    
    BOOST_CHECK(size(overlapped1) == 1);
    BOOST_CHECK(size(overlapped2) == 0);
    BOOST_CHECK(size(overlapped3) == 1);
    BOOST_CHECK(size(overlapped4) == 1);
}

BOOST_AUTO_TEST_CASE(inner_distance_respects_insertion_lhs_ordering_rule)
{
    Variant snp1 {"chr1", 99, "T", "A"};
    Variant snp2 {"chr1", 101, "T", "C"};
    Variant ins {"chr1", 100, "", "A"};
    
    BOOST_CHECK(inner_distance(snp1, ins) == 0);
    BOOST_CHECK(inner_distance(ins, snp1) == 0);
    BOOST_CHECK(inner_distance(snp2, ins) == -1);
    BOOST_CHECK(inner_distance(ins, snp2) == 1);
}

BOOST_AUTO_TEST_CASE(indels_can_be_left_aligned)
{
    auto human = make_reference(human_reference_fasta);
    
    // Huntingtin region CCAGCAGCAGCAGCAG...
    auto a_region = parse_region("4:3076657-3076660", human);
    
    auto the_sequence = human.get_sequence(a_region);
    
    BOOST_CHECK(the_sequence == "CAG");
    
    Variant a_deletion {a_region, the_sequence, ""};
    
    Variant left_aligned_deletion = left_align(a_deletion, human);
    
    BOOST_CHECK(left_aligned_deletion.get_region() ==
            parse_region("4:3076603-3076606", human));
    BOOST_CHECK(left_aligned_deletion.get_reference_allele_sequence() == "CAG");
    BOOST_CHECK(left_aligned_deletion.get_alternative_allele_sequence() == "");
    
    Variant an_insertion {parse_region("4:3076660-3076660", human),
                                               "", the_sequence};
    
    auto left_aligned_insertion = left_align(an_insertion, human);
    
    BOOST_CHECK(left_aligned_insertion.get_region() ==
            parse_region("4:3076603-3076603", human));
    BOOST_CHECK(left_aligned_insertion.get_reference_allele_sequence() == "");
    BOOST_CHECK(left_aligned_insertion.get_alternative_allele_sequence() == "CAG");
    
    // Region is CCAACAACAACAACAC (94594947-94594962)
    a_region = parse_region("5:94594956-94594959", human);
    
    the_sequence = human.get_sequence(a_region);
    
    BOOST_CHECK(the_sequence == "CAA");
    
    a_deletion = Variant {a_region, the_sequence, ""};
    
    left_aligned_deletion = left_align(a_deletion, human);
    
    BOOST_CHECK(left_aligned_deletion.get_region() ==
            parse_region("5:94594949-94594952", human));
    BOOST_CHECK(left_aligned_deletion.get_reference_allele_sequence() == "ACA");
    BOOST_CHECK(left_aligned_deletion.get_alternative_allele_sequence() == "");
    
    an_insertion = Variant {parse_region("5:94594959-94594959", human), "", the_sequence};
    
    left_aligned_insertion = left_align(an_insertion, human);
    
    BOOST_CHECK(left_aligned_insertion.get_region() ==
            parse_region("5:94594949-94594949", human));
    BOOST_CHECK(left_aligned_insertion.get_reference_allele_sequence() == "");
    BOOST_CHECK(left_aligned_insertion.get_alternative_allele_sequence() == "ACA");
}

BOOST_AUTO_TEST_CASE(can_make_variants_parsimonious)
{
    auto human = make_reference(human_reference_fasta);
    
    Variant a_snp {parse_region("12:10001330-10001331", human), std::string {"G"}, std::string {"C"}};
    
    BOOST_CHECK(is_parsimonious(a_snp));
    BOOST_CHECK(make_parsimonious(a_snp, human) == a_snp);
    
    Variant an_unparsimonious_snp {parse_region("12:10001330-10001332", human), std::string {"GT"}, std::string {"CT"}};
    
    BOOST_CHECK(!is_parsimonious(an_unparsimonious_snp));
    
    auto parsimonised_snp = make_parsimonious(an_unparsimonious_snp, human);
    
    BOOST_CHECK(is_parsimonious(parsimonised_snp));
    BOOST_CHECK(parsimonised_snp == a_snp);
    
    Variant another_unparsimonious_snp {parse_region("12:10001329-10001332", human), std::string {"TGT"}, std::string {"TCT"}};
    
    BOOST_CHECK(!is_parsimonious(another_unparsimonious_snp));
    
    auto another_parsimonised_snp = make_parsimonious(another_unparsimonious_snp, human);
    
    BOOST_CHECK(is_parsimonious(another_parsimonised_snp));
    BOOST_CHECK(another_parsimonised_snp == a_snp);
    
    auto a_region = parse_region("12:10001330-10001335", human);
    
    auto the_sequence = human.get_sequence(a_region);
    
    BOOST_CHECK(the_sequence == "GTGGA");
    
    Variant a_deletion {a_region, the_sequence, ""};
    
    auto parsimonious_deletion = make_parsimonious(a_deletion, human);
    
    BOOST_CHECK(parsimonious_deletion.get_region() ==
            parse_region("12:10001329-10001335", human));
    BOOST_CHECK(parsimonious_deletion.get_reference_allele_sequence() == "CGTGGA");
    BOOST_CHECK(parsimonious_deletion.get_alternative_allele_sequence() == "C");
    
    Variant an_insertion {parse_region("12:10001330-10001330", human), "", the_sequence};
    
    auto parsimonious_insertion = make_parsimonious(an_insertion, human);
    
    BOOST_CHECK(parsimonious_insertion.get_region() ==
            parse_region("12:10001329-10001330", human));
    BOOST_CHECK(parsimonious_insertion.get_reference_allele_sequence() == "C");
    BOOST_CHECK(parsimonious_insertion.get_alternative_allele_sequence() == "CGTGGA");
    
    Variant an_unparsimonious_deletion {parse_region("12:10001328-10001335", human), "TCGTGGA", "TC"};
    
    BOOST_CHECK(!is_parsimonious(an_unparsimonious_deletion));
    
    auto parsimonised_deletion = make_parsimonious(an_unparsimonious_deletion, human);
    
    BOOST_CHECK(is_parsimonious(parsimonised_deletion));
    
    Variant an_unparsimonious_insertion {parse_region("12:10001329-10001331", human), "CG", "CGTGGA"};
    
    BOOST_CHECK(!is_parsimonious(an_unparsimonious_insertion));
    
    auto parsimonised_insertion = make_parsimonious(an_unparsimonious_insertion, human);
    
    BOOST_CHECK(is_parsimonious(parsimonised_insertion));
}

BOOST_AUTO_TEST_CASE(can_normalise_variants)
{
    auto human = make_reference(human_reference_fasta);
    
    // Huntingtin region CCAGCAGCAGCAGCAG...
    
    Variant a_snp {parse_region("4:3076657-3076658", human), std::string {"G"}, std::string {"C"}};
    
    BOOST_CHECK(is_parsimonious(a_snp));
    
    auto a_normalised_snp = normalise(a_snp, human);
    
    BOOST_CHECK(a_normalised_snp == a_snp);
    
    auto a_region = parse_region("4:3076657-3076660", human);
    
    auto the_sequence = human.get_sequence(a_region);
    
    BOOST_CHECK(the_sequence == "CAG");
    
    Variant a_mnp {a_region, the_sequence, std::string {"GAC"}};
    
    BOOST_CHECK(is_parsimonious(a_mnp));
    
    auto a_normalised_mnp = normalise(a_mnp, human);
    
    BOOST_CHECK(a_normalised_mnp == a_mnp);
    
    Variant a_deletion {a_region, the_sequence, ""};
    
    BOOST_CHECK(!is_parsimonious(a_deletion));
    
    auto left_aligned_unparsimonious_deletion = left_align(a_deletion, human);
    
    BOOST_CHECK(!is_parsimonious(left_aligned_unparsimonious_deletion));
    
    auto normilised_deletion = normalise(a_deletion, human);
    
    BOOST_CHECK(is_parsimonious(normilised_deletion));
    BOOST_CHECK(normilised_deletion.get_region() ==
            parse_region("4:3076602-3076606", human));
    BOOST_CHECK(normilised_deletion.get_reference_allele_sequence() == "CCAG");
    BOOST_CHECK(normilised_deletion.get_alternative_allele_sequence() == "C");
    
    Variant an_insertion {parse_region("4:3076660-3076660", human), "", the_sequence};
    
    BOOST_CHECK(!is_parsimonious(an_insertion));
    
    auto left_aligned_unparsimonious_insertion = left_align(an_insertion, human);
    
    BOOST_CHECK(!is_parsimonious(left_aligned_unparsimonious_insertion));
    
    auto normilised_insertion = normalise(an_insertion, human);
    
    BOOST_CHECK(is_parsimonious(normilised_insertion));
    BOOST_CHECK(normilised_insertion.get_region() ==
            parse_region("4:3076602-3076603", human));
    BOOST_CHECK(normilised_insertion.get_reference_allele_sequence() == "C");
    BOOST_CHECK(normilised_insertion.get_alternative_allele_sequence() == "CCAG");
    
    // Some hard ones
    
    Variant an_unormilised_snp {parse_region("4:3076656-3076659", human), std::string {"AGC"}, std::string {"ACC"}};
    
    BOOST_CHECK(!is_parsimonious(an_unormilised_snp));
    
    a_normalised_snp = normalise(an_unormilised_snp, human);
    
    BOOST_CHECK(is_parsimonious(a_normalised_mnp));
    BOOST_CHECK(a_normalised_snp.get_region() == parse_region("4:3076657-3076658", human));
    BOOST_CHECK(a_normalised_snp.get_reference_allele_sequence() == "G");
    BOOST_CHECK(a_normalised_snp.get_alternative_allele_sequence() == "C");
    
    Variant an_unormilised_mnp {parse_region("4:3076656-3076661", human), std::string {"GCAGC"}, std::string {"GGACC"}};
    
    BOOST_CHECK(!is_parsimonious(an_unormilised_mnp));
    
    a_normalised_mnp = normalise(an_unormilised_mnp, human);
    
    BOOST_CHECK(is_parsimonious(a_normalised_mnp));
    BOOST_CHECK(a_normalised_mnp.get_region() == parse_region("4:3076657-3076660", human));
    BOOST_CHECK(a_normalised_mnp.get_reference_allele_sequence() == "CAG");
    BOOST_CHECK(a_normalised_mnp.get_alternative_allele_sequence() == "GAC");
    
    Variant an_unnormilised_deletion {parse_region("4:3076655-3076660", human), std::string {"AGCAG"}, std::string {"AG"}};
    
    BOOST_CHECK(!is_parsimonious(an_unnormilised_deletion));
    
    auto a_normalised_deletion = normalise(an_unnormilised_deletion, human);
    
    BOOST_CHECK(is_parsimonious(a_normalised_deletion));
    BOOST_CHECK(a_normalised_deletion.get_region() ==
            parse_region("4:3076602-3076606", human));
    BOOST_CHECK(a_normalised_deletion.get_reference_allele_sequence() == "CCAG");
    BOOST_CHECK(a_normalised_deletion.get_alternative_allele_sequence() == "C");
    
    Variant an_unnormilised_insertion {parse_region("4:3076655-3076657", human),
                                    std::string {"AG"}, std::string {"AGCAG"}};
    
    BOOST_CHECK(!is_parsimonious(an_unnormilised_insertion));
    
    auto a_normalised_insertion = normalise(an_unnormilised_insertion, human);
    
    BOOST_CHECK(is_parsimonious(a_normalised_insertion));
    BOOST_CHECK(a_normalised_insertion.get_region() ==
            parse_region("4:3076602-3076603", human));
    BOOST_CHECK(a_normalised_insertion.get_reference_allele_sequence() == "C");
    BOOST_CHECK(a_normalised_insertion.get_alternative_allele_sequence() == "CCAG");
}

//BOOST_AUTO_TEST_CASE(decompose_returns_all_alleles_from_given_variants)
//{
//    std::vector<Variant> variants {};
//}

BOOST_AUTO_TEST_SUITE_END()
