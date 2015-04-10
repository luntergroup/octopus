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
#include "variant.h"
#include "variant_factory.h"
#include "variant_utils.h"
#include "mock_objects.h"

TEST_CASE("comparison_consistency_test", "[variant]")
{
    VariantFactory a_variant_factory {};
    
    auto snp1 = a_variant_factory.make("chr1", 100, "C", "A");
    auto snp2 = a_variant_factory.make("chr1", 99, "C", "A");
    auto snp3 = a_variant_factory.make("chr1", 100, "C", "T");
    auto snp4 = a_variant_factory.make("chr1", 100, "C", "A");
    
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
}

TEST_CASE("snp_overlap_test", "[snps]")
{
    VariantFactory a_variant_factory {};
    
    auto snp1 = a_variant_factory.make("chr1", 100, "C", "A");
    auto snp2 = a_variant_factory.make("chr1", 99, "C", "A");
    auto snp3 = a_variant_factory.make("chr1", 101, "C", "A");
    auto snp4 = a_variant_factory.make("chr1", 100, "C", "T");
    auto snp5 = a_variant_factory.make("chr1", 99, "C", "T");
    auto snp6 = a_variant_factory.make("chr1", 101, "C", "T");
    
    REQUIRE(overlaps(snp1, snp1));
    REQUIRE(overlaps(snp1, snp4));
    REQUIRE(!overlaps(snp1, snp2));
    REQUIRE(!overlaps(snp1, snp3));
    REQUIRE(!overlaps(snp1, snp5));
    REQUIRE(!overlaps(snp1, snp6));
}

TEST_CASE("mnp_overlap_test", "[mnp]")
{
    VariantFactory a_variant_factory {};
    
    auto mnp1 = a_variant_factory.make("chr1", 100, "CAT", "TAC");
    auto mnp2 = a_variant_factory.make("chr1", 99, "CAT", "TAC");
    auto mnp3 = a_variant_factory.make("chr1", 101, "CAT", "TAC");
    auto mnp4 = a_variant_factory.make("chr1", 100, "CAT", "TAG");
    auto mnp5 = a_variant_factory.make("chr1", 99, "CAT", "TAG");
    auto mnp6 = a_variant_factory.make("chr1", 101, "CAT", "TAG");
    
    // edge cases
    auto mnp7  = a_variant_factory.make("chr1", 98, "CAT", "TAC");
    auto mnp8  = a_variant_factory.make("chr1", 102, "CAT", "TAC");
    auto mnp9  = a_variant_factory.make("chr1", 97, "CAT", "TAC");
    auto mnp10 = a_variant_factory.make("chr1", 103, "CAT", "TAC");
    
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

TEST_CASE("insertion_overlap_test", "[insertion]")
{
    VariantFactory a_variant_factory {};
    
    auto insert1 = a_variant_factory.make("chr1", 100, "", "TAG");
    auto insert2 = a_variant_factory.make("chr1", 99, "", "TAG");
    auto insert3 = a_variant_factory.make("chr1", 101, "", "TAG");
    
    // Insertions never overlap. IS THIS RIGHT!?!
    REQUIRE(!overlaps(insert1, insert1));
    REQUIRE(!overlaps(insert1, insert2));
    REQUIRE(!overlaps(insert1, insert3));
}

TEST_CASE("deletion_overlap_test", "[deletion]")
{
    VariantFactory a_variant_factory {};
    
    auto del1 = a_variant_factory.make("chr1", 100, "TAG", "");
    auto del2 = a_variant_factory.make("chr1", 99, "TAG", "");
    auto del3 = a_variant_factory.make("chr1", 101, "TAG", "");
    
    REQUIRE(overlaps(del1, del1));
    REQUIRE(overlaps(del1, del2));
    REQUIRE(overlaps(del1, del3));
}

TEST_CASE("variant_ordering_test", "[variant]")
{
    VariantFactory a_variant_factory {};
    
    auto snp1 = a_variant_factory.make("chr1", 100, "T", "A");
    auto snp2 = a_variant_factory.make("chr1", 100, "T", "C");
    auto snp3 = a_variant_factory.make("chr1", 100, "T", "G");
    auto ins1 = a_variant_factory.make("chr1", 100, "", "AG");
    auto ins2 = a_variant_factory.make("chr1", 100, "", "CC");
    auto ins3 = a_variant_factory.make("chr1", 100, "", "CCA");
    auto del1 = a_variant_factory.make("chr1", 100, "TA", "");
    auto del2 = a_variant_factory.make("chr1", 100, "TAG", "");
    
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

TEST_CASE("left_alignment", "[left_alignment]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    VariantFactory a_variant_factory {};
    
    // Huntingtin region CCAGCAGCAGCAGCAG...
    auto a_region = parse_region("4:3076657-3076660", human);
    
    auto the_sequence = human.get_sequence(a_region);
    
    REQUIRE(the_sequence == "CAG");
    
    auto a_deletion = a_variant_factory.make(a_region, the_sequence, "");
    
    auto left_aligned_deletion = left_align(a_deletion, human, a_variant_factory);
    
    REQUIRE(left_aligned_deletion.get_region() ==
            parse_region("4:3076603-3076606", human));
    REQUIRE(left_aligned_deletion.get_reference_allele() == "CAG");
    REQUIRE(left_aligned_deletion.get_alternative_allele() == "");
    
    auto an_insertion = a_variant_factory.make(parse_region("4:3076660-3076660", human),
                                               "", the_sequence);
    
    auto left_aligned_insertion = left_align(an_insertion, human, a_variant_factory);
    
    REQUIRE(left_aligned_insertion.get_region() ==
            parse_region("4:3076603-3076603", human));
    REQUIRE(left_aligned_insertion.get_reference_allele() == "");
    REQUIRE(left_aligned_insertion.get_alternative_allele() == "CAG");
    
    // Region is CCAACAACAACAACAC (94594947-94594962)
    a_region = parse_region("5:94594956-94594959", human);
    
    the_sequence = human.get_sequence(a_region);
    
    REQUIRE(the_sequence == "CAA");
    
    a_deletion = a_variant_factory.make(a_region, the_sequence, "");
    
    left_aligned_deletion = left_align(a_deletion, human, a_variant_factory);
    
    REQUIRE(left_aligned_deletion.get_region() ==
            parse_region("5:94594949-94594952", human));
    REQUIRE(left_aligned_deletion.get_reference_allele() == "ACA");
    REQUIRE(left_aligned_deletion.get_alternative_allele() == "");
    
    an_insertion = a_variant_factory.make(parse_region("5:94594959-94594959", human),
                                               "", the_sequence);
    
    left_aligned_insertion = left_align(an_insertion, human, a_variant_factory);
    
    REQUIRE(left_aligned_insertion.get_region() ==
            parse_region("5:94594949-94594949", human));
    REQUIRE(left_aligned_insertion.get_reference_allele() == "");
    REQUIRE(left_aligned_insertion.get_alternative_allele() == "ACA");
}

TEST_CASE("parsimonious_test", "[parsimonious]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    VariantFactory a_variant_factory {};
    
    auto a_snp = a_variant_factory.make(parse_region("12:10001330-10001331", human),
                                        std::string {"G"}, std::string {"C"});
    
    REQUIRE(is_parsimonious(a_snp));
    REQUIRE(make_parsimonious(a_snp, human, a_variant_factory) == a_snp);
    
    auto an_unparsimonious_snp = a_variant_factory.make(parse_region("12:10001330-10001332", human),
                                                        std::string {"GT"}, std::string {"CT"});
    
    REQUIRE(!is_parsimonious(an_unparsimonious_snp));
    
    auto parsimonised_snp = make_parsimonious(an_unparsimonious_snp, human, a_variant_factory);
    
    REQUIRE(is_parsimonious(parsimonised_snp));
    REQUIRE(parsimonised_snp == a_snp);
    
    auto another_unparsimonious_snp = a_variant_factory.make(parse_region("12:10001329-10001332", human),
                                                             std::string {"TGT"}, std::string {"TCT"});
    
    REQUIRE(!is_parsimonious(another_unparsimonious_snp));
    
    auto another_parsimonised_snp = make_parsimonious(another_unparsimonious_snp, human, a_variant_factory);
    
    REQUIRE(is_parsimonious(another_parsimonised_snp));
    REQUIRE(another_parsimonised_snp == a_snp);
    
    auto a_region = parse_region("12:10001330-10001335", human);
    
    auto the_sequence = human.get_sequence(a_region);
    
    REQUIRE(the_sequence == "GTGGA");
    
    auto a_deletion = a_variant_factory.make(a_region, the_sequence, "");
    
    auto parsimonious_deletion = make_parsimonious(a_deletion, human, a_variant_factory);
    
    REQUIRE(parsimonious_deletion.get_region() ==
            parse_region("12:10001329-10001335", human));
    REQUIRE(parsimonious_deletion.get_reference_allele() == "CGTGGA");
    REQUIRE(parsimonious_deletion.get_alternative_allele() == "C");
    
    auto an_insertion = a_variant_factory.make(parse_region("12:10001330-10001330", human),
                                               "", the_sequence);
    
    auto parsimonious_insertion = make_parsimonious(an_insertion, human, a_variant_factory);
    
    REQUIRE(parsimonious_insertion.get_region() ==
            parse_region("12:10001329-10001330", human));
    REQUIRE(parsimonious_insertion.get_reference_allele() == "C");
    REQUIRE(parsimonious_insertion.get_alternative_allele() == "CGTGGA");
    
    auto an_unparsimonious_deletion = a_variant_factory.make(parse_region("12:10001328-10001335", human),
                                                             "TCGTGGA", "TC");
    
    REQUIRE(!is_parsimonious(an_unparsimonious_deletion));
    
    auto parsimonised_deletion = make_parsimonious(an_unparsimonious_deletion, human, a_variant_factory);
    
    REQUIRE(is_parsimonious(parsimonised_deletion));
    
    auto an_unparsimonious_insertion = a_variant_factory.make(parse_region("12:10001329-10001331", human),
                                                              "CG", "CGTGGA");
    
    REQUIRE(!is_parsimonious(an_unparsimonious_insertion));
    
    auto parsimonised_insertion = make_parsimonious(an_unparsimonious_insertion, human, a_variant_factory);
    
    REQUIRE(is_parsimonious(parsimonised_insertion));
}

TEST_CASE("normalisation_test", "[normalisation]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    VariantFactory a_variant_factory {};
    
    // Huntingtin region CCAGCAGCAGCAGCAG...
    
    auto a_snp = a_variant_factory.make(parse_region("4:3076657-3076658", human),
                                        std::string {"G"}, std::string {"C"});
    
    REQUIRE(is_parsimonious(a_snp));
    
    auto a_normalised_snp = normalise(a_snp, human, a_variant_factory);
    
    REQUIRE(a_normalised_snp == a_snp);
    
    auto a_region = parse_region("4:3076657-3076660", human);
    
    auto the_sequence = human.get_sequence(a_region);
    
    REQUIRE(the_sequence == "CAG");
    
    auto a_mnp = a_variant_factory.make(a_region, the_sequence, std::string {"GAC"});
    
    REQUIRE(is_parsimonious(a_mnp));
    
    auto a_normalised_mnp = normalise(a_mnp, human, a_variant_factory);
    
    REQUIRE(a_normalised_mnp == a_mnp);
    
    auto a_deletion = a_variant_factory.make(a_region, the_sequence, "");
    
    REQUIRE(!is_parsimonious(a_deletion));
    
    auto left_aligned_unparsimonious_deletion = left_align(a_deletion, human, a_variant_factory);
    
    REQUIRE(!is_parsimonious(left_aligned_unparsimonious_deletion));
    
    auto normilised_deletion = normalise(a_deletion, human, a_variant_factory);
    
    REQUIRE(is_parsimonious(normilised_deletion));
    REQUIRE(normilised_deletion.get_region() ==
            parse_region("4:3076602-3076606", human));
    REQUIRE(normilised_deletion.get_reference_allele() == "CCAG");
    REQUIRE(normilised_deletion.get_alternative_allele() == "C");
    
    auto an_insertion = a_variant_factory.make(parse_region("4:3076660-3076660", human),
                                               "", the_sequence);
    
    REQUIRE(!is_parsimonious(an_insertion));
    
    auto left_aligned_unparsimonious_insertion = left_align(an_insertion, human, a_variant_factory);
    
    REQUIRE(!is_parsimonious(left_aligned_unparsimonious_insertion));
    
    auto normilised_insertion = normalise(an_insertion, human, a_variant_factory);
    
    REQUIRE(is_parsimonious(normilised_insertion));
    REQUIRE(normilised_insertion.get_region() ==
            parse_region("4:3076602-3076603", human));
    REQUIRE(normilised_insertion.get_reference_allele() == "C");
    REQUIRE(normilised_insertion.get_alternative_allele() == "CCAG");
    
    // Some hard ones
    
    auto an_unormilised_snp = a_variant_factory.make(parse_region("4:3076656-3076659", human),
                                                     std::string {"AGC"}, std::string {"ACC"});
    
    REQUIRE(!is_parsimonious(an_unormilised_snp));
    
    a_normalised_snp = normalise(an_unormilised_snp, human, a_variant_factory);
    
    REQUIRE(is_parsimonious(a_normalised_mnp));
    REQUIRE(a_normalised_snp.get_region() == parse_region("4:3076657-3076658", human));
    REQUIRE(a_normalised_snp.get_reference_allele() == "G");
    REQUIRE(a_normalised_snp.get_alternative_allele() == "C");
    
    auto an_unormilised_mnp = a_variant_factory.make(parse_region("4:3076656-3076661", human),
                                                     std::string {"GCAGC"}, std::string {"GGACC"});
    
    REQUIRE(!is_parsimonious(an_unormilised_mnp));
    
    a_normalised_mnp = normalise(an_unormilised_mnp, human, a_variant_factory);
    
    REQUIRE(is_parsimonious(a_normalised_mnp));
    REQUIRE(a_normalised_mnp.get_region() == parse_region("4:3076657-3076660", human));
    REQUIRE(a_normalised_mnp.get_reference_allele() == "CAG");
    REQUIRE(a_normalised_mnp.get_alternative_allele() == "GAC");
    
    auto an_unnormilised_deletion = a_variant_factory.make(parse_region("4:3076655-3076660", human),
                                                          std::string {"AGCAG"}, std::string {"AG"});
    
    REQUIRE(!is_parsimonious(an_unnormilised_deletion));
    
    auto a_normalised_deletion = normalise(an_unnormilised_deletion, human, a_variant_factory);
    
    REQUIRE(is_parsimonious(a_normalised_deletion));
    REQUIRE(a_normalised_deletion.get_region() ==
            parse_region("4:3076602-3076606", human));
    REQUIRE(a_normalised_deletion.get_reference_allele() == "CCAG");
    REQUIRE(a_normalised_deletion.get_alternative_allele() == "C");
    
    auto an_unnormilised_insertion = a_variant_factory.make(parse_region("4:3076655-3076657", human),
                                                          std::string {"AG"}, std::string {"AGCAG"});
    
    REQUIRE(!is_parsimonious(an_unnormilised_insertion));
    
    auto a_normalised_insertion = normalise(an_unnormilised_insertion, human, a_variant_factory);
    
    REQUIRE(is_parsimonious(a_normalised_insertion));
    REQUIRE(a_normalised_insertion.get_region() ==
            parse_region("4:3076602-3076603", human));
    REQUIRE(a_normalised_insertion.get_reference_allele() == "C");
    REQUIRE(a_normalised_insertion.get_alternative_allele() == "CCAG");
}
