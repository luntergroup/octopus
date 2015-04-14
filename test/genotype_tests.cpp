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
#include "variant_candidate_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"

TEST_CASE("get_all_genotypes gives in all possible unique genotypes", "[genotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    Haplotype hap1 {human};
    hap1.push_back(Allele {parse_region("1", human), ""});
    
    Haplotype hap2 {human};
    hap2.push_back(Allele {parse_region("2", human), ""});
    
    Haplotype hap3 {human};
    hap3.push_back(Allele {parse_region("3", human), ""});
    
    Haplotype hap4 {human};
    hap4.push_back(Allele {parse_region("4", human), ""});
    
    std::vector<Haplotype> haplotypes {hap1, hap2, hap3, hap4};
    
    unsigned num_haplotypes {4};
    
    auto genotypes_1 = get_all_genotypes(haplotypes, 1);
    
    REQUIRE(genotypes_1.size() == num_genotypes(num_haplotypes, 1));
    
    std::unordered_set<Genotype> unique_1 {genotypes_1.cbegin(), genotypes_1.cend()};
    
    REQUIRE(genotypes_1.size() == unique_1.size());
    
    auto genotypes_2 = get_all_genotypes(haplotypes, 2);
    
    REQUIRE(genotypes_2.size() == num_genotypes(num_haplotypes, 2));
    
    std::unordered_set<Genotype> unique_2 {genotypes_2.cbegin(), genotypes_2.cend()};
    
    REQUIRE(genotypes_2.size() == unique_2.size());
    
    auto genotypes_3 = get_all_genotypes(haplotypes, 3);
    
    REQUIRE(genotypes_3.size() == num_genotypes(num_haplotypes, 3));
    
    std::unordered_set<Genotype> unique_3 {genotypes_3.cbegin(), genotypes_3.cend()};
    
    REQUIRE(genotypes_3.size() == unique_3.size());
    
    auto genotypes_4 = get_all_genotypes(haplotypes, 4);
    
    REQUIRE(genotypes_4.size() == num_genotypes(num_haplotypes, 4));
    
    std::unordered_set<Genotype> unique_4 {genotypes_4.cbegin(), genotypes_4.cend()};
    
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
