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

#include "reference_genome.h"
#include "reference_genome_implementor_factory.h"

TEST_CASE("initialisation_test", "[reference]")
{
    ReferenceGenomeImplementorFactory a_factory;
    
    std::string homedir {getenv("HOME")};
    
    // test for a small single genome contig
    std::string ecoli_fasta_path {homedir + "/Genomics/References/R00000042.fasta"};
    ReferenceGenome ecoli(a_factory.make(ecoli_fasta_path));
    
    REQUIRE(ecoli.get_name() == "R00000042");
    REQUIRE(ecoli.contains_region(GenomicRegion("R00000042", 10000, 2000000)));
    REQUIRE(ecoli.get_contig_size("R00000042") == 5231428);
    REQUIRE(!ecoli.has_contig("X"));
    REQUIRE(ecoli.get_contig_region("R00000042") == GenomicRegion("R00000042", 0, 5231428));
    REQUIRE(ecoli.get_sequence(GenomicRegion("R00000042", 0, 10)) == "AGCTTTTCAT"); // first line
    REQUIRE(ecoli.get_sequence(GenomicRegion("R00000042", 69, 80)) == "CTTCTGAACTG"); // accross lines
    
    // test for a large multi-contig genome
    std::string human_fasta_path {homedir + "/Genomics/References/human_g1k_v37.fasta"};
    ReferenceGenome human(a_factory.make(human_fasta_path));
    
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
