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
    std::string fasta_path {homedir + "/Genomics/References/human_g1k_v37.fasta"};
    
    ReferenceGenome reference(a_factory.make(fasta_path));
    
    REQUIRE(reference.get_name() == "NOT IMPLEMENTED!");
    REQUIRE(reference.contains_region(GenomicRegion("1", 100, 10000)));
    REQUIRE(reference.get_contig_size("20") == 63025520);
    REQUIRE(reference.has_contig("X"));
    REQUIRE(reference.get_contig_region("X") == GenomicRegion("X", 0, 155270560));
    REQUIRE(reference.get_sequence(GenomicRegion("15", 51265690, 51265700)) == "ATACAATGTT");
}
