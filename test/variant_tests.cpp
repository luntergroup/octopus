//
//  variant_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include <iostream>

#include "catch.hpp"

#include "genomic_region.h"
#include "variant.h"
#include "variant_factory.h"

TEST_CASE("snp_overlap_test", "[snps]")
{
    auto snp1 = VariantFactory::get_instance().make("chr1", 100, "C", "A");
    auto snp2 = VariantFactory::get_instance().make("chr1", 99, "C", "A");
    auto snp3 = VariantFactory::get_instance().make("chr1", 101, "C", "A");
    auto snp4 = VariantFactory::get_instance().make("chr1", 100, "C", "T");
    auto snp5 = VariantFactory::get_instance().make("chr1", 99, "C", "T");
    auto snp6 = VariantFactory::get_instance().make("chr1", 101, "C", "T");
    
    REQUIRE(overlaps(snp1, snp1));
    REQUIRE(overlaps(snp1, snp4));
    REQUIRE(!overlaps(snp1, snp2));
    REQUIRE(!overlaps(snp1, snp3));
    REQUIRE(!overlaps(snp1, snp5));
    REQUIRE(!overlaps(snp1, snp6));
}

TEST_CASE("mnp_overlap_test", "[mnp]")
{
    auto mnp1 = VariantFactory::get_instance().make("chr1", 100, "CAT", "TAC");
    auto mnp2 = VariantFactory::get_instance().make("chr1", 99, "CAT", "TAC");
    auto mnp3 = VariantFactory::get_instance().make("chr1", 101, "CAT", "TAC");
    auto mnp4 = VariantFactory::get_instance().make("chr1", 100, "CAT", "TAG");
    auto mnp5 = VariantFactory::get_instance().make("chr1", 99, "CAT", "TAG");
    auto mnp6 = VariantFactory::get_instance().make("chr1", 101, "CAT", "TAG");
    
    // edge cases
    auto mnp7  = VariantFactory::get_instance().make("chr1", 98, "CAT", "TAC");
    auto mnp8  = VariantFactory::get_instance().make("chr1", 102, "CAT", "TAC");
    auto mnp9  = VariantFactory::get_instance().make("chr1", 97, "CAT", "TAC");
    auto mnp10 = VariantFactory::get_instance().make("chr1", 103, "CAT", "TAC");
    
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
    auto insert1 = VariantFactory::get_instance().make("chr1", 100, "", "TAG");
    auto insert2 = VariantFactory::get_instance().make("chr1", 99, "", "TAG");
    auto insert3 = VariantFactory::get_instance().make("chr1", 101, "", "TAG");
    
    // Insertions never overlap. IS THIS RIGHT!?!
    REQUIRE(!overlaps(insert1, insert1));
    REQUIRE(!overlaps(insert1, insert2));
    REQUIRE(!overlaps(insert1, insert3));
}

TEST_CASE("deletion_overlap_test", "[deletion]")
{
    auto del1 = VariantFactory::get_instance().make("chr1", 100, "TAG", "");
    auto del2 = VariantFactory::get_instance().make("chr1", 99, "TAG", "");
    auto del3 = VariantFactory::get_instance().make("chr1", 101, "TAG", "");
    
    REQUIRE(overlaps(del1, del1));
    REQUIRE(overlaps(del1, del2));
    REQUIRE(overlaps(del1, del3));
}
