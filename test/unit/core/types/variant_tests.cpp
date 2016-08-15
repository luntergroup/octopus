// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <vector>
#include <algorithm>
#include <iterator>

#include <basics/genomic_region.hpp>
#include <core/types/allele.hpp>
#include <core/types/variant.hpp>
#include <utils/mappable_algorithms.hpp>

#include "mock/mock_reference.hpp"

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(core)
BOOST_AUTO_TEST_SUITE(variant)

static void sort(std::vector<Variant>& variants)
{
    std::sort(std::begin(variants), std::end(variants));
}

BOOST_AUTO_TEST_CASE(operators_less_and_equal_are_consistent)
{
    const Variant snp1 {"test", 100, "C", "A"};
    const Variant snp2 {"test", 99, "C", "A"};
    const Variant snp3 {"test", 100, "C", "T"};
    const Variant snp4 {"test", 100, "C", "A"};
    
    bool r1 = !(snp1 < snp2) && !(snp2 < snp1);
    bool r2 = (snp1 == snp2);
    BOOST_REQUIRE(!r2);
    BOOST_CHECK_EQUAL(r1, r2);
    
    r1 = !(snp1 < snp3) && !(snp3 < snp1);
    r2 = snp1 == snp3;
    BOOST_REQUIRE(!r2);
    BOOST_CHECK_EQUAL(r1, r2);
    
    r1 = !(snp1 < snp4) && !(snp4 < snp1);
    r2 = snp1 == snp4;
    BOOST_REQUIRE(r2);
    BOOST_CHECK_EQUAL(r1, r2);
    
    const Variant del1 {"test", 100, "C", ""};
    const Variant del2 {"test", 100, "CA", ""};
    const Variant del3 {"test", 100, "C", ""};
    
    r1 = !(del1 < del2) && !(del2 < del1);
    r2 = del1 == del2;
    BOOST_REQUIRE(!r2);
    BOOST_CHECK_EQUAL(r1, r2);
    
    r1 = !(del1 < del3) && !(del3 < del1);
    r2 = del1 == del3;
    BOOST_REQUIRE(r2);
    BOOST_CHECK_EQUAL(r1, r2);
    
    const Variant ins1 {"test", 100, "", "C"};
    const Variant ins2 {"test", 100, "", "CA"};
    const Variant ins3 {"test", 100, "", "C"};
    
    r1 = !(ins1 < ins2) && !(ins2 < ins1);
    r2 = ins1 == ins2;
    BOOST_REQUIRE(!r2);
    BOOST_CHECK_EQUAL(r1, r2);
    
    r1 = !(ins1 < ins3) && !(ins3 < ins1);
    r2 = ins1 == ins3;
    BOOST_REQUIRE(r2);
    BOOST_CHECK_EQUAL(r1, r2);
}

BOOST_AUTO_TEST_CASE(can_binary_search_for_ordered_variants)
{
    const Variant snp1 {"test", 100, "C", "A"};
    const Variant snp2 {"test", 105, "G", "T"};
    const Variant del1 {"test", 103, "C", ""};
    const Variant del2 {"test", 115, "T", ""};
    const Variant ins1 {"test", 107, "", "G"};
    const Variant ins2 {"test", 110, "", "CA"};
    
    std::vector<Variant> variants {snp1, del1, snp2, ins1, ins2, del2};
    
    sort(variants);
    
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), snp1));
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), snp2));
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), del1));
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), del2));
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), ins1));
    BOOST_CHECK(std::binary_search(variants.cbegin(), variants.cend(), ins2));
    
    const Variant non_snp1 {"test", 111, "C", "A"};
    const Variant non_snp2 {"test", 105, "G", "A"};
    const Variant non_del1 {"test", 104, "C", ""};
    const Variant non_del2 {"test", 115, "TA", ""};
    const Variant non_ins1 {"test", 108, "", "G"};
    const Variant non_ins2 {"test", 110, "", "TA"};
    
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_snp1));
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_snp2));
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_del1));
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_del2));
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_ins1));
    BOOST_CHECK(!std::binary_search(variants.cbegin(), variants.cend(), non_ins2));
    
    const Variant del3 {"test", 101, "C", ""};
    const Variant snp3 {"test", 101, "C", ""};
    
    std::vector<Variant> variants2 {snp1, del3, snp3};
    
    sort(variants2);
    
    BOOST_CHECK(std::binary_search(variants2.cbegin(), variants2.cend(), snp1));
    BOOST_CHECK(std::binary_search(variants2.cbegin(), variants2.cend(), del3));
    BOOST_CHECK(std::binary_search(variants2.cbegin(), variants2.cend(), snp3));
}

BOOST_AUTO_TEST_CASE(snps_do_not_overlap_adjacent_snps)
{
    const Variant snp1 {"test", 100, "C", "A"};
    const Variant snp2 {"test", 99, "C", "A"};
    const Variant snp3 {"test", 101, "C", "A"};
    const Variant snp4 {"test", 100, "C", "T"};
    const Variant snp5 {"test", 99, "C", "T"};
    const Variant snp6 {"test", 101, "C", "T"};
    
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
    const Variant mnp1 {"test", 100, "CAT", "TAC"};
    const Variant mnp2 {"test", 99, "CAT", "TAC"};
    const Variant mnp3 {"test", 101, "CAT", "TAC"};
    const Variant mnp4 {"test", 100, "CAT", "TAG"};
    const Variant mnp5 {"test", 99, "CAT", "TAG"};
    const Variant mnp6 {"test", 101, "CAT", "TAG"};
    
    // edge cases
    const Variant mnp7  {"test", 98, "CAT", "TAC"};
    const Variant mnp8  {"test", 102, "CAT", "TAC"};
    const Variant mnp9  {"test", 97, "CAT", "TAC"};
    const Variant mnp10 {"test", 103, "CAT", "TAC"};
    
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
    const Variant insert1 {"test", 1, "", "TAG"};
    const Variant insert2 {"test", 0, "", "TAG"};
    const Variant insert3 {"test", 2, "", "TAG"};
    
    BOOST_CHECK(overlaps(insert1, insert1));
    BOOST_CHECK(!overlaps(insert1, insert2));
    BOOST_CHECK(!overlaps(insert1, insert3));
}

BOOST_AUTO_TEST_CASE(insertions_overlap_with_other_variants_when_contained_by_their_region)
{
    const Variant ins {"test", 1, "",    "TAG"};
    const Variant del {"test", 0, "TAG", ""};
    const Variant mnp {"test", 0, "TAG", "CAT"};
    
    BOOST_CHECK(overlaps(ins, del));
    BOOST_CHECK(overlaps(ins, mnp));
}

BOOST_AUTO_TEST_CASE(variants_are_ordered_by_region_and_then_lexicographically_by_sequence)
{
    const Variant snp1 {"test", 0, "T", "A"};
    const Variant snp2 {"test", 0, "T", "C"};
    const Variant snp3 {"test", 0, "T", "G"};
    const Variant ins1 {"test", 0, "", "AG"};
    const Variant ins2 {"test", 0, "", "CC"};
    const Variant ins3 {"test", 0, "", "CCA"};
    const Variant del1 {"test", 0, "TA", ""};
    const Variant del2 {"test", 0, "TAG", ""};
    
    std::vector<Variant> variants1 {snp1, ins1, snp2};
    
    sort(variants1);
    
    std::vector<Variant> variants1_required_sort {ins1, snp1, snp2};
    
    bool is_required_sort1 = std::equal(variants1.cbegin(), variants1.cend(), variants1_required_sort.cbegin());
    
    BOOST_CHECK(is_required_sort1);
    
    std::vector<Variant> variants2 {snp1, del1, snp2, ins1, snp3};
    
    sort(variants2);
    
    std::vector<Variant> variants2_required_sort {ins1, snp1, snp2, snp3, del1};
    
    auto is_required_sort2 = std::equal(variants2.cbegin(), variants2.cend(), variants2_required_sort.cbegin());
    
    BOOST_CHECK(is_required_sort2);
    
    std::vector<Variant> variants3 {del2, snp1, del1, snp2, ins1, snp3, ins2};
    
    sort(variants3);
    
    std::vector<Variant> variants3_required_sort {ins1, ins2, snp1, snp2, snp3, del1, del2};
    
    auto is_required_sort3 = std::equal(variants3.cbegin(), variants3.cend(), variants3_required_sort.cbegin());
    
    BOOST_CHECK(is_required_sort3);
    
    std::vector<Variant> variants4 {ins2, del2, snp1, del1, snp2, ins1, snp3, ins3};
    
    sort(variants4);
    
    std::vector<Variant> variants4_required_sort {ins1, ins2, ins3, snp1, snp2, snp3, del1, del2};
    
    auto is_required_sort4 = std::equal(variants4.cbegin(), variants4.cend(), variants4_required_sort.cbegin());
    
    BOOST_CHECK(is_required_sort4);
}

BOOST_AUTO_TEST_CASE(insertions_are_overlapped_on_the_left_and_right_boundry)
{
    const Variant snp1 {"test", 0, "T", "A"};
    const Variant ins1 {"test", 5, "", "A"};
    const Variant snp2 {"test", 0, "T", "C"};
    
    std::vector<Variant> variants {snp1, ins1, snp2};
    
    sort(variants);
    
    const GenomicRegion region1 {"test", 4, 4};
    const GenomicRegion region2 {"test", 4, 5};
    const GenomicRegion region3 {"test", 4, 6};
    const GenomicRegion region4 {"test", 5, 5};
    const GenomicRegion region5 {"test", 5, 6};
    const GenomicRegion region6 {"test", 6, 6};
    
    BOOST_CHECK_EQUAL(count_overlapped(variants, region1), 0);
    BOOST_CHECK_EQUAL(count_overlapped(variants, region2), 1);
    BOOST_CHECK_EQUAL(count_overlapped(variants, region3), 1);
    BOOST_CHECK_EQUAL(count_overlapped(variants, region4), 1);
    BOOST_CHECK_EQUAL(count_overlapped(variants, region5), 1);
    BOOST_CHECK_EQUAL(count_overlapped(variants, region6), 0);
    
    const Variant del1 {"test", 0, "TTTTTTTTTT", ""};
    
    variants.push_back(del1);
    
    sort(variants);
    
    BOOST_CHECK_EQUAL(count_overlapped(variants, region1), 1);
    BOOST_CHECK_EQUAL(count_overlapped(variants, region2), 2);
    BOOST_CHECK_EQUAL(count_overlapped(variants, region3), 2);
    BOOST_CHECK_EQUAL(count_overlapped(variants, region4), 2);
    BOOST_CHECK_EQUAL(count_overlapped(variants, region5), 2);
    BOOST_CHECK_EQUAL(count_overlapped(variants, region6), 1);
}

BOOST_AUTO_TEST_CASE(inner_distance_respects_insertion_lhs_ordering_rule)
{
    const Variant snp1 {"test", 0, "T", "A"};
    const Variant ins  {"test", 1, "",  "A"};
    const Variant snp2 {"test", 2, "T", "C"};
    
    BOOST_CHECK_EQUAL(inner_distance(snp1, ins), 0);
    BOOST_CHECK_EQUAL(inner_distance(ins, snp1), 0);
    BOOST_CHECK_EQUAL(inner_distance(snp2, ins), -1);
    BOOST_CHECK_EQUAL(inner_distance(ins, snp2), 1);
}

BOOST_AUTO_TEST_CASE(indels_can_be_left_aligned)
{
    const auto reference = make_mock_reference();
    
    auto region = parse_region("mock4:657-660", reference); // sequence around region is CCAGCAGCAGCAGCAG...
    
    auto sequence = reference.fetch_sequence(region);
    
    BOOST_REQUIRE_EQUAL(sequence, "CAG");
    
    Variant raw_deletion {region, sequence, ""};
    
    Variant left_aligned_deletion = left_align(raw_deletion, reference);
    
    BOOST_CHECK_EQUAL(left_aligned_deletion.mapped_region(), parse_region("mock4:603-606", reference));
    BOOST_CHECK_EQUAL(ref_sequence(left_aligned_deletion), "CAG");
    BOOST_CHECK_EQUAL(alt_sequence(left_aligned_deletion), "");
    
    Variant raw_insertion {parse_region("mock4:660-660", reference), "", sequence};
    
    auto left_aligned_insertion = left_align(raw_insertion, reference);
    
    BOOST_CHECK_EQUAL(left_aligned_insertion.mapped_region(), parse_region("mock4:603-603", reference));
    BOOST_CHECK_EQUAL(ref_sequence(left_aligned_insertion), "");
    BOOST_CHECK_EQUAL(alt_sequence(left_aligned_insertion), "CAG");
    
    region = parse_region("mock5:956-959", reference); // sequence around region is CCAACAACAACAACAC
    
    sequence = reference.fetch_sequence(region);
    
    BOOST_REQUIRE_EQUAL(sequence, "CAA");
    
    raw_deletion = Variant {region, sequence, ""};
    
    left_aligned_deletion = left_align(raw_deletion, reference);
    
    BOOST_CHECK_EQUAL(left_aligned_deletion.mapped_region(), parse_region("mock5:949-952", reference));
    BOOST_CHECK_EQUAL(ref_sequence(left_aligned_deletion), "ACA");
    BOOST_CHECK_EQUAL(alt_sequence(left_aligned_deletion), "");
    
    raw_insertion = Variant {parse_region("mock5:959-959", reference), "", sequence};
    
    left_aligned_insertion = left_align(raw_insertion, reference);
    
    BOOST_CHECK_EQUAL(left_aligned_insertion.mapped_region(), parse_region("mock5:949-949", reference));
    BOOST_CHECK_EQUAL(ref_sequence(left_aligned_insertion), "");
    BOOST_CHECK_EQUAL(alt_sequence(left_aligned_insertion), "ACA");
}

BOOST_AUTO_TEST_CASE(can_make_variants_parsimonious)
{
    const auto reference = make_mock_reference();
    
    Variant raw_snp {parse_region("mock6:330-331", reference), std::string {"G"}, std::string {"C"}};
    
    BOOST_CHECK(is_parsimonious(raw_snp));
    BOOST_CHECK_EQUAL(make_parsimonious(raw_snp, reference), raw_snp);
    
    Variant unparsimonious_snp {parse_region("mock6:330-332", reference), std::string {"GT"}, std::string {"CT"}};
    
    BOOST_REQUIRE(!is_parsimonious(unparsimonious_snp));
    
    auto parsimonised_snp = make_parsimonious(unparsimonious_snp, reference);
    
    BOOST_CHECK(is_parsimonious(parsimonised_snp));
    BOOST_CHECK_EQUAL(parsimonised_snp, raw_snp);
    
    Variant another_unparsimonious_snp {parse_region("mock6:329-332", reference), std::string {"TGT"}, std::string {"TCT"}};
    
    BOOST_REQUIRE(!is_parsimonious(another_unparsimonious_snp));
    
    auto another_parsimonised_snp = make_parsimonious(another_unparsimonious_snp, reference);
    
    BOOST_CHECK(is_parsimonious(another_parsimonised_snp));
    BOOST_CHECK_EQUAL(another_parsimonised_snp, raw_snp);
    
    auto region = parse_region("mock6:330-335", reference);
    
    auto sequence = reference.fetch_sequence(region);
    
    BOOST_REQUIRE_EQUAL(sequence, "GTGGA");
    
    Variant raw_deletion {region, sequence, ""};
    
    auto parsimonious_deletion = make_parsimonious(raw_deletion, reference);
    
    BOOST_CHECK_EQUAL(parsimonious_deletion.mapped_region(), parse_region("mock6:329-335", reference));
    BOOST_CHECK_EQUAL(ref_sequence(parsimonious_deletion), "CGTGGA");
    BOOST_CHECK_EQUAL(alt_sequence(parsimonious_deletion), "C");
    
    Variant raw_insertion {parse_region("mock6:330-330", reference), "", sequence};
    
    auto parsimonious_insertion = make_parsimonious(raw_insertion, reference);
    
    BOOST_CHECK_EQUAL(parsimonious_insertion.mapped_region(), parse_region("mock6:329-330", reference));
    BOOST_CHECK_EQUAL(ref_sequence(parsimonious_insertion), "C");
    BOOST_CHECK_EQUAL(alt_sequence(parsimonious_insertion), "CGTGGA");
    
    Variant unparsimonious_deletion {parse_region("mock6:328-335", reference), "TCGTGGA", "TC"};
    
    BOOST_REQUIRE(!is_parsimonious(unparsimonious_deletion));
    
    auto parsimonised_deletion = make_parsimonious(unparsimonious_deletion, reference);
    
    BOOST_CHECK(is_parsimonious(parsimonised_deletion));
    
    Variant unparsimonious_insertion {parse_region("mock6:329-331", reference), "CG", "CGTGGA"};
    
    BOOST_REQUIRE(!is_parsimonious(unparsimonious_insertion));
    
    auto parsimonised_insertion = make_parsimonious(unparsimonious_insertion, reference);
    
    BOOST_CHECK(is_parsimonious(parsimonised_insertion));
}

BOOST_AUTO_TEST_CASE(can_normalise_variants)
{
   const auto reference = make_mock_reference();
    
    Variant raw_snp {parse_region("mock4:657-658", reference), std::string {"G"}, std::string {"C"}};
    
    BOOST_REQUIRE(is_parsimonious(raw_snp));
    
    auto normalised_snp = normalise(raw_snp, reference);
    
    BOOST_CHECK_EQUAL(normalised_snp, raw_snp);
    
    auto region = parse_region("mock4:657-660", reference);
    
    auto sequence = reference.fetch_sequence(region);
    
    BOOST_CHECK_EQUAL(sequence, "CAG");
    
    Variant raw_mnp {region, sequence, std::string {"GAC"}};
    
    BOOST_REQUIRE(is_parsimonious(raw_mnp));
    
    auto normalised_mnp = normalise(raw_mnp, reference);
    
    BOOST_CHECK_EQUAL(normalised_mnp, raw_mnp);
    
    Variant raw_deletion {region, sequence, ""};
    
    BOOST_REQUIRE(!is_parsimonious(raw_deletion));
    
    auto left_aligned_unparsimonious_deletion = left_align(raw_deletion, reference);
    
    BOOST_REQUIRE(!is_parsimonious(left_aligned_unparsimonious_deletion));
    
    auto normilised_deletion = normalise(raw_deletion, reference);
    
    BOOST_CHECK(is_parsimonious(normilised_deletion));
    BOOST_CHECK_EQUAL(normilised_deletion.mapped_region(), parse_region("mock4:602-606", reference));
    BOOST_CHECK_EQUAL(ref_sequence(normilised_deletion), "CCAG");
    BOOST_CHECK_EQUAL(alt_sequence(normilised_deletion), "C");
    
    Variant raw_insertion {parse_region("mock4:660-660", reference), "", sequence};
    
    BOOST_REQUIRE(!is_parsimonious(raw_insertion));
    
    auto left_aligned_unparsimonious_insertion = left_align(raw_insertion, reference);
    
    BOOST_REQUIRE(!is_parsimonious(left_aligned_unparsimonious_insertion));
    
    auto normilised_insertion = normalise(raw_insertion, reference);
    
    BOOST_CHECK(is_parsimonious(normilised_insertion));
    BOOST_CHECK_EQUAL(normilised_insertion.mapped_region(), parse_region("mock4:602-603", reference));
    BOOST_CHECK_EQUAL(ref_sequence(normilised_insertion), "C");
    BOOST_CHECK_EQUAL(alt_sequence(normilised_insertion), "CCAG");
    
    // Some hard ones
    
    Variant unormilised_snp {parse_region("mock4:656-659", reference), std::string {"AGC"}, std::string {"ACC"}};
    
    BOOST_REQUIRE(!is_parsimonious(unormilised_snp));
    
    normalised_snp = normalise(unormilised_snp, reference);
    
    BOOST_CHECK(is_parsimonious(normalised_mnp));
    BOOST_CHECK_EQUAL(normalised_snp.mapped_region(), parse_region("mock4:657-658", reference));
    BOOST_CHECK_EQUAL(ref_sequence(normalised_snp), "G");
    BOOST_CHECK_EQUAL(alt_sequence(normalised_snp), "C");
    
    Variant unormilised_mnp {parse_region("mock4:656-661", reference), std::string {"GCAGC"}, std::string {"GGACC"}};
    
    BOOST_REQUIRE(!is_parsimonious(unormilised_mnp));
    
    normalised_mnp = normalise(unormilised_mnp, reference);
    
    BOOST_CHECK(is_parsimonious(normalised_mnp));
    BOOST_CHECK_EQUAL(normalised_mnp.mapped_region(), parse_region("mock4:657-660", reference));
    BOOST_CHECK_EQUAL(ref_sequence(normalised_mnp), "CAG");
    BOOST_CHECK_EQUAL(alt_sequence(normalised_mnp), "GAC");
    
    Variant unnormilised_deletion {parse_region("mock4:655-660", reference), std::string {"AGCAG"}, std::string {"AG"}};
    
    BOOST_REQUIRE(!is_parsimonious(unnormilised_deletion));
    
    auto normalised_deletion = normalise(unnormilised_deletion, reference);
    
    BOOST_CHECK(is_parsimonious(normalised_deletion));
    BOOST_CHECK_EQUAL(normalised_deletion.mapped_region(), parse_region("mock4:602-606", reference));
    BOOST_CHECK_EQUAL(ref_sequence(normalised_deletion), "CCAG");
    BOOST_CHECK_EQUAL(alt_sequence(normalised_deletion), "C");
    
    Variant unnormilised_insertion {parse_region("mock4:655-657", reference),
        std::string {"AG"}, std::string {"AGCAG"}};
    
    BOOST_REQUIRE(!is_parsimonious(unnormilised_insertion));
    
    auto normalised_insertion = normalise(unnormilised_insertion, reference);
    
    BOOST_CHECK(is_parsimonious(normalised_insertion));
    BOOST_CHECK_EQUAL(normalised_insertion.mapped_region(), parse_region("mock4:602-603", reference));
    BOOST_CHECK_EQUAL(ref_sequence(normalised_insertion), "C");
    BOOST_CHECK_EQUAL(alt_sequence(normalised_insertion), "CCAG");
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
