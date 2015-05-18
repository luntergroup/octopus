//
//  genotype_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>
#include <cstddef>
#include <set>
#include <vector>

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "read_manager.h"
#include "variant.h"
#include "variant_utils.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"

using std::cout;
using std::endl;

TEST_CASE("Genotype can be tested for haplotype occurence", "[genotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    Haplotype hap1 {human};
    hap1.push_back(Allele {parse_region("1:1000000-1000001", human), "A"});
    
    Haplotype hap2 {human};
    hap2.push_back(Allele {parse_region("1:1000000-1000001", human), "C"});
    
    Haplotype hap3 {human};
    hap3.push_back(Allele {parse_region("1:1000000-1000001", human), "G"});
    
    Haplotype hap4 {human};
    hap4.push_back(Allele {parse_region("1:1000000-1000001", human), "T"});
    
    Genotype<Haplotype> g1 {};
    g1.emplace(hap1);
    g1.emplace(hap2);
    g1.emplace(hap3);
    
    REQUIRE(g1.contains(hap1));
    REQUIRE(g1.contains(hap2));
    REQUIRE(g1.contains(hap3));
    REQUIRE(!g1.contains(hap4));
    
    REQUIRE(g1.num_occurences(hap1) == 1);
    REQUIRE(g1.num_occurences(hap2) == 1);
    REQUIRE(g1.num_occurences(hap3) == 1);
    REQUIRE(g1.num_occurences(hap4) == 0);
    
    Genotype<Haplotype> g2 {};
    g2.emplace(hap1);
    g2.emplace(hap1);
    g2.emplace(hap2);
    
    REQUIRE(g2.contains(hap1));
    REQUIRE(g2.contains(hap2));
    REQUIRE(!g2.contains(hap3));
    REQUIRE(!g2.contains(hap4));
    
    REQUIRE(g2.num_occurences(hap1) == 2);
    REQUIRE(g2.num_occurences(hap2) == 1);
    REQUIRE(g2.num_occurences(hap3) == 0);
    REQUIRE(g2.num_occurences(hap4) == 0);
    
    Genotype<Haplotype> g3 {};
    g3.emplace(hap1);
    g3.emplace(hap3);
    g3.emplace(hap4);
    
    REQUIRE(g3.contains(hap1));
    REQUIRE(!g3.contains(hap2));
    REQUIRE(g3.contains(hap3));
    REQUIRE(g3.contains(hap4));
    
    REQUIRE(g3.num_occurences(hap1) == 1);
    REQUIRE(g3.num_occurences(hap2) == 0);
    REQUIRE(g3.num_occurences(hap3) == 1);
    REQUIRE(g3.num_occurences(hap4) == 1);
    
    Genotype<Haplotype> g4 {};
    g4.emplace(hap4);
    g4.emplace(hap4);
    g4.emplace(hap4);
    
    REQUIRE(!g4.contains(hap1));
    REQUIRE(!g4.contains(hap2));
    REQUIRE(!g4.contains(hap3));
    REQUIRE(g4.contains(hap4));
    
    REQUIRE(g4.num_occurences(hap1) == 0);
    REQUIRE(g4.num_occurences(hap2) == 0);
    REQUIRE(g4.num_occurences(hap3) == 0);
    REQUIRE(g4.num_occurences(hap4) == 3);
}

TEST_CASE("Genotypes are equal if they contain the same haplotypes", "[genotype")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
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
    
    REQUIRE(g1 == g1);
    REQUIRE(g1 != g2);
    REQUIRE(g1 == g3);
    REQUIRE(g1 != g4);
    REQUIRE(g1 != g5);
    
    REQUIRE(g2 == g2);
    REQUIRE(g2 != g3);
    REQUIRE(g2 != g4);
    REQUIRE(g2 == g5);
    
    REQUIRE(g3 == g3);
    REQUIRE(g3 != g4);
    REQUIRE(g3 != g5);
    
    REQUIRE(g4 == g4);
    REQUIRE(g4 != g5);
    
    REQUIRE(g5 == g5);
}

TEST_CASE("Genotypes are not influenced by haplotype entry order", "[genotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
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
    
    REQUIRE(g1.num_occurences(hap1) == g2.num_occurences(hap1));
    REQUIRE(g1.num_occurences(hap2) == g2.num_occurences(hap2));
    
    REQUIRE(g1 == g2);
    REQUIRE(std::hash<Genotype<Haplotype>>()(g1) == std::hash<Genotype<Haplotype>>()(g2));
}

TEST_CASE("get_all_genotypes returns all possible unique genotypes", "[genotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
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
    
    auto genotypes_1 = get_all_genotypes(haplotypes, 1);
    
    REQUIRE(genotypes_1.size() == num_genotypes(num_haplotypes, 1));
    
    std::unordered_set<Genotype<Haplotype>> unique_1 {genotypes_1.cbegin(), genotypes_1.cend()};
    
    REQUIRE(genotypes_1.size() == unique_1.size());
    
    auto genotypes_2 = get_all_genotypes(haplotypes, 2);
    
    REQUIRE(genotypes_2.size() == num_genotypes(num_haplotypes, 2));
    
    std::unordered_set<Genotype<Haplotype>> unique_2 {genotypes_2.cbegin(), genotypes_2.cend()};
    
    REQUIRE(genotypes_2.size() == unique_2.size());
    
    auto genotypes_3 = get_all_genotypes(haplotypes, 3);
    
    REQUIRE(genotypes_3.size() == num_genotypes(num_haplotypes, 3));
    
    std::unordered_set<Genotype<Haplotype>> unique_3 {genotypes_3.cbegin(), genotypes_3.cend()};
    
    REQUIRE(genotypes_3.size() == unique_3.size());
    
    auto genotypes_4 = get_all_genotypes(haplotypes, 4);
    
    REQUIRE(genotypes_4.size() == num_genotypes(num_haplotypes, 4));
    
    std::unordered_set<Genotype<Haplotype>> unique_4 {genotypes_4.cbegin(), genotypes_4.cend()};
    
    REQUIRE(genotypes_4.size() == unique_4.size());
}

TEST_CASE("get_all_genotypes results in correct ploidy", "[genotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
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
    
    auto genotypes1 = get_all_genotypes(haplotypes, 1);
    
    for (const auto& genotype : genotypes1) {
        REQUIRE(genotype.ploidy() == 1);
    }
    
    auto genotypes2 = get_all_genotypes(haplotypes, 2);
    
    for (const auto& genotype : genotypes2) {
        REQUIRE(genotype.ploidy() == 2);
    }
    
    auto genotypes3 = get_all_genotypes(haplotypes, 3);
    
    for (const auto& genotype : genotypes3) {
        REQUIRE(genotype.ploidy() == 3);
    }
}

TEST_CASE("Genotype::get_unique returns all the unique haplotypes in a Genotype", "[genotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
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
    
    REQUIRE(std::count(g1_unique.cbegin(), g1_unique.cend(), hap1) == 1);
    REQUIRE(std::count(g1_unique.cbegin(), g1_unique.cend(), hap2) == 1);
    REQUIRE(std::count(g1_unique.cbegin(), g1_unique.cend(), hap3) == 1);
    
    Genotype<Haplotype> g2 {};
    g2.emplace(hap1);
    g2.emplace(hap3);
    g2.emplace(hap3);
    
    auto g2_unique = g2.get_unique();
    
    REQUIRE(std::count(g2_unique.cbegin(), g2_unique.cend(), hap1) == 1);
    REQUIRE(std::count(g2_unique.cbegin(), g2_unique.cend(), hap2) == 0);
    REQUIRE(std::count(g2_unique.cbegin(), g2_unique.cend(), hap3) == 1);
    
    Genotype<Haplotype> g3 {};
    g3.emplace(hap3);
    g3.emplace(hap3);
    g3.emplace(hap3);
    
    auto g3_unique = g3.get_unique();
    
    REQUIRE(std::count(g3_unique.cbegin(), g3_unique.cend(), hap1) == 0);
    REQUIRE(std::count(g3_unique.cbegin(), g3_unique.cend(), hap2) == 0);
    REQUIRE(std::count(g3_unique.cbegin(), g3_unique.cend(), hap3) == 1);
}
