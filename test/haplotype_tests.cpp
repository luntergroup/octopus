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
    
}

TEST_CASE("test_extend_haplotypes", "[haplotype]")
{
    
}

TEST_CASE("test_make_genotypes", "[genotype]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    Haplotype hap1 {human};
    hap1.emplace_back(parse_region("1", human), "");
    
    Haplotype hap2 {human};
    hap2.emplace_back(parse_region("2", human), "");
    
    Haplotype hap3 {human};
    hap3.emplace_back(parse_region("3", human), "");
    
    Haplotype hap4 {human};
    hap4.emplace_back(parse_region("4", human), "");
    
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
