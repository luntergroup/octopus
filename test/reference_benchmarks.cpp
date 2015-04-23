//
//  reference_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include "benchmark_utils.h"

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "bioio.h"

TEST_CASE("reference_benchmark", "[reference,benchmark]")
{
    std::ifstream fasta(human_reference_fasta);
    auto index = bioio::read_fasta_index(human_reference_fasta_index)["5"];
    auto f_bioio = [&fasta, &index] () {
        bioio::read_fasta_contig(fasta, index, 100000, 100);
    };
    
    auto without_vptr = benchmark<std::chrono::nanoseconds>(f_bioio, 10).count();
    
    fasta.close();
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome reference(a_factory.make(human_reference_fasta));
    
    auto f_ref = [&reference] () {
        reference.get_sequence(GenomicRegion("5", 100000, 100100));
    };
    
    auto with_vptr = benchmark<std::chrono::nanoseconds>(f_ref, 10).count();
    
    //std::cout << "Without vptr: " << without_vptr << "ns" << std::endl;
    //std::cout << "With vptr: " << with_vptr << "ns" << std::endl;
    
    REQUIRE(true);
}
