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
    
    auto f = a_factory.make(fasta_path);
    
    auto f2 = std::move(f);
    
    std::cout << f2->get_reference_name() << std::endl;
    std::cout << f2->get_contig_size("1") << std::endl;
    
    ReferenceGenome reference(std::move(f));
    
    REQUIRE(reference.get_name() == "NOT IMPLEMENTED!");
    REQUIRE(reference.contains_region(GenomicRegion("1", 100, 10000)));
}
