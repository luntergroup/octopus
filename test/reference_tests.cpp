//
//  reference_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(ReferemceGenome_handles_basic_queries)
{
    ReferenceGenomeFactory a_factory {};
    
    ReferenceGenome ecoli {a_factory.make(ecoli_reference_fasta)};
    
    BOOST_CHECK(ecoli.get_name() == "R00000042");
    BOOST_CHECK(ecoli.contains_region(GenomicRegion("R00000042", 10000, 2000000)));
    BOOST_CHECK(ecoli.get_contig_size("R00000042") == 5231428);
    BOOST_CHECK(!ecoli.has_contig("X"));
    BOOST_CHECK(ecoli.get_contig_region("R00000042") == GenomicRegion("R00000042", 0, 5231428));
    BOOST_CHECK(ecoli.get_sequence(GenomicRegion("R00000042", 0, 10)) == "AGCTTTTCAT"); // first line
    BOOST_CHECK(ecoli.get_sequence(GenomicRegion("R00000042", 69, 80)) == "CTTCTGAACTG"); // accross lines
    
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    BOOST_CHECK(human.get_name() == "human_g1k_v37");
    BOOST_CHECK(human.contains_region(GenomicRegion("1", 100, 10000)));
    BOOST_CHECK(!human.contains_region(GenomicRegion("1", 100, 1e10))); // too big
    BOOST_CHECK(human.get_contig_size("20") == 63025520);
    BOOST_CHECK(human.has_contig("X"));
    BOOST_CHECK(!human.has_contig("y"));
    BOOST_CHECK(human.get_contig_region("X") == GenomicRegion("X", 0, 155270560));
    BOOST_CHECK(human.get_sequence(GenomicRegion("15", 51265690, 51265700)) == "ACAATGTTGT");
    BOOST_CHECK(human.get_sequence(GenomicRegion("5", 100000, 100010)) == "AGGAAGTTTC");
}

BOOST_AUTO_TEST_CASE(ReferemceGenome_handles_edge_cases)
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    BOOST_CHECK(human.get_sequence(GenomicRegion {"1", 100, 100}) == "");
}

BOOST_AUTO_TEST_CASE(parse_region_works_with_correctly_formatted_region_input)
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    auto r1 = parse_region("3", human);
    BOOST_CHECK(r1.get_contig_name() == "3");
    BOOST_CHECK(r1.get_begin() == 0);
    BOOST_CHECK(r1.get_end() == human.get_contig_size("3"));
    
    auto r2 = parse_region("10:100-200", human);
    BOOST_CHECK(r2.get_contig_name() == "10");
    BOOST_CHECK(r2.get_begin() == 100);
    BOOST_CHECK(r2.get_end() == 200);
    
    auto r3 = parse_region("18:102029", human);
    BOOST_CHECK(r3.get_contig_name() == "18");
    BOOST_CHECK(r3.get_begin() == 102029);
    BOOST_CHECK(r3.get_end() == 102029);
    
    auto r4 = parse_region("MT:100-", human);
    BOOST_CHECK(r4.get_contig_name() == "MT");
    BOOST_CHECK(r4.get_begin() == 100);
    BOOST_CHECK(r4.get_end() == human.get_contig_size("MT"));
}

BOOST_AUTO_TEST_CASE(parse_region_throws_when_region_is_not_formatted_correctly)
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
    
    BOOST_CHECK(throwed);
    
    try {
        auto r2 = parse_region("5:100-99", human);
        throwed = false;
    } catch (...) {
        throwed = true;
    }
    
    BOOST_CHECK(throwed);
    
//    try {
//        auto r3 = parse_region("2::0-100", human);
//        throwed = false;
//    } catch (...) {
//        throwed = true;
//    }
//    
//    BOOST_CHECK(throwed);
}

BOOST_AUTO_TEST_SUITE_END()
