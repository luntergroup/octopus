//
//  genotype_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <cstddef>
#include <set>
#include <vector>

#include "test_common.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "variant.hpp"
#include "variant_utils.hpp"
#include "candidate_variant_generator.hpp"
#include "alignment_candidate_variant_generator.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(can_iterate_Genotype_Haplotypes_with_range_based_for)
{
    auto human = make_reference(human_reference_fasta);
    
    Haplotype hap1 {human};
    hap1.push_back(Allele {parse_region("1:1000000-1000001", human), "A"});
    
    Haplotype hap2 {human};
    hap2.push_back(Allele {parse_region("1:1000000-1000001", human), "C"});
    
    Genotype<Haplotype> genotype {hap1, hap2};
    
    std::vector<Haplotype> r {};
    
    for (const auto& haplotype : genotype) {
        r.push_back(haplotype);
    }
    
    BOOST_CHECK(r.front() == hap1);
    BOOST_CHECK(r.back() == hap2);
}

BOOST_AUTO_TEST_CASE(Genotype_can_be_tested_for_haplotype_occurence)
{
    auto human = make_reference(human_reference_fasta);
    
    Haplotype hap1 {human};
    hap1.push_back(Allele {parse_region("1:1000000-1000001", human), "A"});
    
    Haplotype hap2 {human};
    hap2.push_back(Allele {parse_region("1:1000000-1000001", human), "C"});
    
    Haplotype hap3 {human};
    hap3.push_back(Allele {parse_region("1:1000000-1000001", human), "G"});
    
    Haplotype hap4 {human};
    hap4.push_back(Allele {parse_region("1:1000000-1000001", human), "T"});
    
    Genotype<Haplotype> g1 {hap1, hap2, hap3};
    
    BOOST_CHECK(g1.contains(hap1));
    BOOST_CHECK(g1.contains(hap2));
    BOOST_CHECK(g1.contains(hap3));
    BOOST_CHECK(!g1.contains(hap4));
    
    BOOST_CHECK(g1.num_occurences(hap1) == 1);
    BOOST_CHECK(g1.num_occurences(hap2) == 1);
    BOOST_CHECK(g1.num_occurences(hap3) == 1);
    BOOST_CHECK(g1.num_occurences(hap4) == 0);
    
    Genotype<Haplotype> g2 {hap1, hap1, hap2};
    
    BOOST_CHECK(g2.contains(hap1));
    BOOST_CHECK(g2.contains(hap2));
    BOOST_CHECK(!g2.contains(hap3));
    BOOST_CHECK(!g2.contains(hap4));
    
    BOOST_CHECK(g2.num_occurences(hap1) == 2);
    BOOST_CHECK(g2.num_occurences(hap2) == 1);
    BOOST_CHECK(g2.num_occurences(hap3) == 0);
    BOOST_CHECK(g2.num_occurences(hap4) == 0);
    
    Genotype<Haplotype> g3 {hap1, hap3, hap4};
    
    BOOST_CHECK(g3.contains(hap1));
    BOOST_CHECK(!g3.contains(hap2));
    BOOST_CHECK(g3.contains(hap3));
    BOOST_CHECK(g3.contains(hap4));
    
    BOOST_CHECK(g3.num_occurences(hap1) == 1);
    BOOST_CHECK(g3.num_occurences(hap2) == 0);
    BOOST_CHECK(g3.num_occurences(hap3) == 1);
    BOOST_CHECK(g3.num_occurences(hap4) == 1);
    
    Genotype<Haplotype> g4 {hap4, hap4, hap4};
    
    BOOST_CHECK(!g4.contains(hap1));
    BOOST_CHECK(!g4.contains(hap2));
    BOOST_CHECK(!g4.contains(hap3));
    BOOST_CHECK(g4.contains(hap4));
    
    BOOST_CHECK(g4.num_occurences(hap1) == 0);
    BOOST_CHECK(g4.num_occurences(hap2) == 0);
    BOOST_CHECK(g4.num_occurences(hap3) == 0);
    BOOST_CHECK(g4.num_occurences(hap4) == 3);
}

BOOST_AUTO_TEST_CASE(Genotypes_are_equal_if_they_contain_the_same_haplotypes)
{
    auto human = make_reference(human_reference_fasta);
    
    Haplotype hap1 {human};
    hap1.push_back(Allele {parse_region("1:1000000-1000001", human), "A"});
    
    Haplotype hap2 {human};
    hap2.push_back(Allele {parse_region("1:1000000-1000001", human), "C"});
    
    Haplotype hap3 {human};
    hap3.push_back(Allele {parse_region("1:1000000-1000001", human), "G"});
    
    Genotype<Haplotype> g1 {};
    g1.emplace(hap1);
    g1.emplace(hap2);
    g1.emplace(hap3);
    
    Genotype<Haplotype> g2 {};
    g2.emplace(hap1);
    g2.emplace(hap2);
    g2.emplace(hap2);
    
    Genotype<Haplotype> g3 {};
    g3.emplace(hap1);
    g3.emplace(hap2);
    g3.emplace(hap3);
    
    Genotype<Haplotype> g4 {};
    g4.emplace(hap1);
    g4.emplace(hap3);
    g4.emplace(hap3);
    
    Genotype<Haplotype> g5 {};
    g5.emplace(hap1);
    g5.emplace(hap2);
    g5.emplace(hap2);
    
    BOOST_CHECK(g1 == g1);
    BOOST_CHECK(g1 != g2);
    BOOST_CHECK(g1 == g3);
    BOOST_CHECK(g1 != g4);
    BOOST_CHECK(g1 != g5);
    
    BOOST_CHECK(g2 == g2);
    BOOST_CHECK(g2 != g3);
    BOOST_CHECK(g2 != g4);
    BOOST_CHECK(g2 == g5);
    
    BOOST_CHECK(g3 == g3);
    BOOST_CHECK(g3 != g4);
    BOOST_CHECK(g3 != g5);
    
    BOOST_CHECK(g4 == g4);
    BOOST_CHECK(g4 != g5);
    
    BOOST_CHECK(g5 == g5);
}

BOOST_AUTO_TEST_CASE(Genotypes_are_not_influenced_by_haplotype_entry_order)
{
    auto human = make_reference(human_reference_fasta);
    
    Haplotype hap1 {human};
    hap1.push_back(Allele {parse_region("1:1000000-1000001", human), "A"});
    
    Haplotype hap2 {human};
    hap2.push_back(Allele {parse_region("1:1000000-1000001", human), "T"});
    
    Genotype<Haplotype> g1 {};
    g1.emplace(hap1);
    g1.emplace(hap2);
    g1.emplace(hap2);
    
    Genotype<Haplotype> g2 {};
    g2.emplace(hap2);
    g2.emplace(hap1);
    g2.emplace(hap2);
    
    BOOST_CHECK(g1.num_occurences(hap1) == g2.num_occurences(hap1));
    BOOST_CHECK(g1.num_occurences(hap2) == g2.num_occurences(hap2));
    
    BOOST_CHECK(g1 == g2);
    BOOST_CHECK(std::hash<Genotype<Haplotype>>()(g1) == std::hash<Genotype<Haplotype>>()(g2));
}

BOOST_AUTO_TEST_CASE(generate_all_genotypes_works_when_num_elements_is_less_than_ploidy)
{
    auto human = make_reference(human_reference_fasta);
    
    Haplotype hap1 {human};
    hap1.push_back(Allele {parse_region("1:1000000-1000001", human), "A"});
    
    Haplotype hap2 {human};
    hap2.push_back(Allele {parse_region("1:1000000-1000001", human), "T"});
    
    std::vector<Haplotype> haplotypes {hap1};
    
    auto genotypes = generate_all_genotypes(haplotypes, 2);
    
    BOOST_CHECK(genotypes.size() == 1);
    
    genotypes = generate_all_genotypes(haplotypes, 3);
    
    BOOST_CHECK(genotypes.size() == 1);
    
    genotypes = generate_all_genotypes(haplotypes, 4);
    
    BOOST_CHECK(genotypes.size() == 1);
    
    genotypes = generate_all_genotypes(haplotypes, 5);
    
    BOOST_CHECK(genotypes.size() == 1);
    
    haplotypes.push_back(hap2);
    
    genotypes = generate_all_genotypes(haplotypes, 3);
    
    BOOST_CHECK(genotypes.size() == 4);
    
    genotypes = generate_all_genotypes(haplotypes, 4);
    
    BOOST_CHECK(genotypes.size() == 5);
}

BOOST_AUTO_TEST_CASE(generate_all_genotypes_returns_all_possible_unique_genotypes)
{
    auto human = make_reference(human_reference_fasta);
    
    Haplotype hap1 {human};
    hap1.push_back(Allele {parse_region("1:1000000-1000001", human), "A"});
    
    Haplotype hap2 {human};
    hap2.push_back(Allele {parse_region("1:1000000-1000001", human), "C"});
    
    Haplotype hap3 {human};
    hap3.push_back(Allele {parse_region("1:1000000-1000001", human), "G"});
    
    Haplotype hap4 {human};
    hap4.push_back(Allele {parse_region("1:1000000-1000001", human), "T"});
    
    std::vector<Haplotype> haplotypes {hap1, hap2, hap3, hap4};
    
    unsigned num_haplotypes {4};
    
    auto genotypes_1 = generate_all_genotypes(haplotypes, 1);
    
    BOOST_CHECK(genotypes_1.size() == num_genotypes(num_haplotypes, 1));
    
    std::unordered_set<Genotype<Haplotype>> unique_1 {genotypes_1.cbegin(), genotypes_1.cend()};
    
    BOOST_CHECK(genotypes_1.size() == unique_1.size());
    
    auto genotypes_2 = generate_all_genotypes(haplotypes, 2);
    
    BOOST_CHECK(genotypes_2.size() == num_genotypes(num_haplotypes, 2));
    
    std::unordered_set<Genotype<Haplotype>> unique_2 {genotypes_2.cbegin(), genotypes_2.cend()};
    
    BOOST_CHECK(genotypes_2.size() == unique_2.size());
    
    auto genotypes_3 = generate_all_genotypes(haplotypes, 3);
    
    BOOST_CHECK(genotypes_3.size() == num_genotypes(num_haplotypes, 3));
    
    std::unordered_set<Genotype<Haplotype>> unique_3 {genotypes_3.cbegin(), genotypes_3.cend()};
    
    BOOST_CHECK(genotypes_3.size() == unique_3.size());
    
    auto genotypes_4 = generate_all_genotypes(haplotypes, 4);
    
    BOOST_CHECK(genotypes_4.size() == num_genotypes(num_haplotypes, 4));
    
    std::unordered_set<Genotype<Haplotype>> unique_4 {genotypes_4.cbegin(), genotypes_4.cend()};
    
    BOOST_CHECK(genotypes_4.size() == unique_4.size());
}

BOOST_AUTO_TEST_CASE(generate_all_genotypes_results_in_correct_ploidy)
{
    auto human = make_reference(human_reference_fasta);
    
    GenomicRegion region1 {"3", 1000000, 1000001};
    GenomicRegion region2 {"3", 1000010, 1000011};
    
    Haplotype haplotype1 {human};
    haplotype1.push_back(Allele{region1, "A"});
    haplotype1.push_back(Allele{region2, "A"});
    Haplotype haplotype2 {human};
    haplotype2.push_back(Allele{region1, "C"});
    haplotype2.push_back(Allele{region2, "C"});
    Haplotype haplotype3 {human};
    haplotype3.push_back(Allele{region1, "G"});
    haplotype3.push_back(Allele{region2, "G"});
    Haplotype haplotype4 {human};
    haplotype4.push_back(Allele{region1, "A"});
    haplotype4.push_back(Allele{region2, "C"});
    Haplotype haplotype5 {human};
    haplotype5.push_back(Allele{region1, "C"});
    haplotype5.push_back(Allele{region2, "G"});
    Haplotype haplotype6 {human};
    haplotype6.push_back(Allele{region1, "G"});
    haplotype6.push_back(Allele{region2, "C"});
    
    std::vector<Haplotype> haplotypes {haplotype1, haplotype2, haplotype3, haplotype4, haplotype5, haplotype6};
    
    auto genotypes1 = generate_all_genotypes(haplotypes, 1);
    
    for (const auto& genotype : genotypes1) {
        BOOST_CHECK(genotype.ploidy() == 1);
    }
    
    auto genotypes2 = generate_all_genotypes(haplotypes, 2);
    
    for (const auto& genotype : genotypes2) {
        BOOST_CHECK(genotype.ploidy() == 2);
    }
    
    auto genotypes3 = generate_all_genotypes(haplotypes, 3);
    
    for (const auto& genotype : genotypes3) {
        BOOST_CHECK(genotype.ploidy() == 3);
    }
}

BOOST_AUTO_TEST_CASE(get_unique_returns_all_the_unique_haplotypes_in_a_Genotype)
{
    auto human = make_reference(human_reference_fasta);
    
    Haplotype hap1 {human};
    hap1.push_back(Allele {parse_region("1:1000000-1000001", human), "A"});
    
    Haplotype hap2 {human};
    hap2.push_back(Allele {parse_region("1:1000000-1000001", human), "C"});
    
    Haplotype hap3 {human};
    hap3.push_back(Allele {parse_region("1:1000000-1000001", human), "G"});
    
    Genotype<Haplotype> g1 {};
    g1.emplace(hap1);
    g1.emplace(hap2);
    g1.emplace(hap3);
    
    auto g1_unique = g1.get_unique();
    
    BOOST_CHECK(std::count(g1_unique.cbegin(), g1_unique.cend(), hap1) == 1);
    BOOST_CHECK(std::count(g1_unique.cbegin(), g1_unique.cend(), hap2) == 1);
    BOOST_CHECK(std::count(g1_unique.cbegin(), g1_unique.cend(), hap3) == 1);
    
    Genotype<Haplotype> g2 {};
    g2.emplace(hap1);
    g2.emplace(hap3);
    g2.emplace(hap3);
    
    auto g2_unique = g2.get_unique();
    
    BOOST_CHECK(std::count(g2_unique.cbegin(), g2_unique.cend(), hap1) == 1);
    BOOST_CHECK(std::count(g2_unique.cbegin(), g2_unique.cend(), hap2) == 0);
    BOOST_CHECK(std::count(g2_unique.cbegin(), g2_unique.cend(), hap3) == 1);
    
    Genotype<Haplotype> g3 {};
    g3.emplace(hap3);
    g3.emplace(hap3);
    g3.emplace(hap3);
    
    auto g3_unique = g3.get_unique();
    
    BOOST_CHECK(std::count(g3_unique.cbegin(), g3_unique.cend(), hap1) == 0);
    BOOST_CHECK(std::count(g3_unique.cbegin(), g3_unique.cend(), hap2) == 0);
    BOOST_CHECK(std::count(g3_unique.cbegin(), g3_unique.cend(), hap3) == 1);
}

BOOST_AUTO_TEST_SUITE_END()
