//
//  reference_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"

TEST_CASE("ReferemceGenome handles basic queries", "[reference]")
{
    ReferenceGenomeFactory a_factory {};
    
    ReferenceGenome ecoli {a_factory.make(ecoli_reference_fasta)};
    
    REQUIRE(ecoli.get_name() == "R00000042");
    REQUIRE(ecoli.contains_region(GenomicRegion("R00000042", 10000, 2000000)));
    REQUIRE(ecoli.get_contig_size("R00000042") == 5231428);
    REQUIRE(!ecoli.has_contig("X"));
    REQUIRE(ecoli.get_contig_region("R00000042") == GenomicRegion("R00000042", 0, 5231428));
    REQUIRE(ecoli.get_sequence(GenomicRegion("R00000042", 0, 10)) == "AGCTTTTCAT"); // first line
    REQUIRE(ecoli.get_sequence(GenomicRegion("R00000042", 69, 80)) == "CTTCTGAACTG"); // accross lines
    
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    REQUIRE(human.get_name() == "human_g1k_v37");
    REQUIRE(human.contains_region(GenomicRegion("1", 100, 10000)));
    REQUIRE(!human.contains_region(GenomicRegion("1", 100, 1e10))); // too big
    REQUIRE(human.get_contig_size("20") == 63025520);
    REQUIRE(human.has_contig("X"));
    REQUIRE(!human.has_contig("y"));
    REQUIRE(human.get_contig_region("X") == GenomicRegion("X", 0, 155270560));
    REQUIRE(human.get_sequence(GenomicRegion("15", 51265690, 51265700)) == "ACAATGTTGT");
    REQUIRE(human.get_sequence(GenomicRegion("5", 100000, 100010)) == "AGGAAGTTTC");
}

TEST_CASE("ReferemceGenome handles edge cases", "[reference]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    REQUIRE(human.get_sequence(GenomicRegion {"1", 100, 100}) == "");
}

TEST_CASE("parse_region works with correctly formatted region input", "[region, reference]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    auto r1 = parse_region("3", human);
    REQUIRE(r1.get_contig_name() == "3");
    REQUIRE(r1.get_begin() == 0);
    REQUIRE(r1.get_end() == human.get_contig_size("3"));
    
    auto r2 = parse_region("10:100-200", human);
    REQUIRE(r2.get_contig_name() == "10");
    REQUIRE(r2.get_begin() == 100);
    REQUIRE(r2.get_end() == 200);
    
    auto r3 = parse_region("18:102029", human);
    REQUIRE(r3.get_contig_name() == "18");
    REQUIRE(r3.get_begin() == 102029);
    REQUIRE(r3.get_end() == 102029);
    
    auto r4 = parse_region("MT:100-", human);
    REQUIRE(r4.get_contig_name() == "MT");
    REQUIRE(r4.get_begin() == 100);
    REQUIRE(r4.get_end() == human.get_contig_size("MT"));
}

TEST_CASE("parse_region throws when region is not formatted correctly", "[region, reference]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    bool throwed {};
    
    try {
        auto r1 = parse_region("-", human);
        throwed = false;
    } catch (...) {
        throwed = true;
    }
    
    REQUIRE(throwed);
    
    try {
        auto r2 = parse_region("5:100-99", human);
        throwed = false;
    } catch (...) {
        throwed = true;
    }
    
    REQUIRE(throwed);
    
//    try {
//        auto r3 = parse_region("2::0-100", human);
//        throwed = false;
//    } catch (...) {
//        throwed = true;
//    }
//    
//    REQUIRE(throwed);
}
