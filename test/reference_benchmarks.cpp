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

#include "reference_genome.h"
#include "reference_genome_implementor_factory.h"
#include "bioio.h"

TEST_CASE("reference_benchmark", "[benchmark]")
{
    std::string homedir {getenv("HOME")};
    std::string human_fasta_path {homedir + "/Genomics/References/human_g1k_v37.fasta"};
    std::string human_fasta_index_path {homedir + "/Genomics/References/human_g1k_v37.fasta.fai"};
    
    std::ifstream fasta(human_fasta_path);
    auto index = bioio::read_fasta_index(human_fasta_index_path)["5"];
    auto f_bioio = [&fasta, &index] () {
        bioio::read_fasta_contig(fasta, index, 100000, 100);
    };
    
    auto bioio = benchmark<std::chrono::nanoseconds>(f_bioio, 1000).count();
    
    fasta.close();
    
    ReferenceGenomeImplementorFactory a_factory;
    ReferenceGenome reference(a_factory.make(human_fasta_path));
    
    auto f_ref = [&reference] () {
        reference.get_sequence(GenomicRegion("5", 100000, 100100));
    };
    
    auto ref = benchmark<std::chrono::nanoseconds>(f_ref, 1000).count();
    
    REQUIRE(ref <= bioio + 3000);
}
