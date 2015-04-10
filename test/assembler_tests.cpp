//
//  assembler_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>

#include "test_common.h"
#include "genomic_region.h"
#include "aligned_read.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "mock_objects.h"
//#include "assembler.h"
#include "kmer_graph.h"
#include "storage_policies.h"

//TEST_CASE("assembler_construct_test", "[assembler]")
//{
//    Assembler assembler {10};
//    
//    auto a_mock_region = get_mock_region();
//    
//    assembler.add_reference_contig(a_mock_region, "AAAAAAAAAACCCCCCCCCC");
//    
//    assembler.add_read(get_mock_aligned_read("AAAAAAAAAACCCCCCCCCC"));
//    
//    REQUIRE(assembler.get_num_verticies() == 10);
//    REQUIRE(assembler.get_num_edges() == 11);
//    
//    assembler.clear();
//    
//    REQUIRE(assembler.get_num_verticies() == 0);
//    REQUIRE(assembler.get_num_edges() == 0);
//    
//    assembler.add_read(get_mock_aligned_read("AAAAAAAAAACCCCCCCCCC"));
//    
//    REQUIRE(assembler.get_num_verticies() == 10);
//    REQUIRE(assembler.get_num_edges() == 11);
//    
//    assembler.add_read(get_mock_aligned_read("GGGGGGGGGGTTTTTTTTTT"));
//    
//    REQUIRE(assembler.get_num_verticies() == 20);
//    REQUIRE(assembler.get_num_edges() == 22);
//    
//    assembler.add_reference_contig(a_mock_region, "CCCCCCCCCCTTTTTTTTTT");
//    
//    REQUIRE(assembler.get_num_verticies() == 28);
//    REQUIRE(assembler.get_num_edges() == 31);
//    
//    assembler.add_reference_contig(a_mock_region, "GGGGGGGGGGAAAAAAAAAA");
//    
//    REQUIRE(assembler.get_num_verticies() == 36);
//    REQUIRE(assembler.get_num_edges() == 40);
//    
//    assembler.add_read(get_mock_aligned_read("CCCCCCCCCCGGGGGGGGGG"));
//    
//    REQUIRE(assembler.get_num_verticies() == 44);
//    REQUIRE(assembler.get_num_edges() == 51);
//    
//    assembler.add_reference_contig(a_mock_region, "CCCCCCCCCCGGGGGGGGGG");
//    
//    REQUIRE(assembler.get_num_verticies() == 44);
//    REQUIRE(assembler.get_num_edges() == 51);
//}

TEST_CASE("assembler_path_test", "[assembler]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome lambda(a_factory.make(lambda_reference_fasta));
    auto contig_name = lambda.get_contig_names()[0];
    auto contig_size = lambda.get_contig_size(contig_name);
    auto contig = lambda.get_sequence(GenomicRegion {contig_name, 0, contig_size});
    
    //Assembler assembler {15};
    
    KmerGraph<int, policies::StoreStringReference> kmer_assembler {10};
    kmer_assembler.add_sequence("AAAAAAAAAACCCCCCCCCC", 0, 1);
    kmer_assembler.add_sequence("CCCCCCCCCCGGGGGGGGGG", 10, 2);
    kmer_assembler.add_sequence("GGGGGGGGGGAAAAAAAAAA", 20, 2);
    //kmer_assembler.print_kmers();
    //kmer_assembler.add_sequence("GGGGGGGGGGTTTTTTTTTT", 1);
//    kmer_assembler.add_sequence(contig, 1);
    auto paths = kmer_assembler.get_contigs(1);
    //std::cout << paths[0] << std::endl;
    //std::cout << kmer_assembler.is_acyclic() << std::endl;
    
    //assembler.add_reference_contig(contig);
    
    //std::cout << assembler.get_contigs().at(0) << std::endl;
}

//TEST_CASE("assembler_cycle_test", "[assembler]")
//{
//    Assembler assembler {5};
//    
//    assembler.add_read("AAAAACCCCC");
//    
//    REQUIRE(!assembler.has_cycle());
//    
//    assembler.add_read("CCCCCGGGGG");
//    
//    //REQUIRE(!assembler.has_cycle());
//    
//    assembler.add_read("GGGGGAAAAA");
//    
//    //REQUIRE(assembler.has_cycle());
//}
