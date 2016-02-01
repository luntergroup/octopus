//
//  assembler_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>

#include "test_common.hpp"
#include "genomic_region.hpp"
#include "aligned_read.hpp"
#include "reference_genome.hpp"
#include "mock_objects.hpp"
//#include "assembler.hpp"
#include "kmer_graph.hpp"
#include "storage_policies.hpp"

BOOST_AUTO_TEST_SUITE(Components)

//BOOST_AUTO_TEST_CASE(assembler_construct_test)
//{
//    Assembler assembler {10};
//    
//    auto a_mock_region = get_mock_region();
//    
//    assembler.add_reference_contig(a_mock_region, "AAAAAAAAAACCCCCCCCCC");
//    
//    assembler.add_read(get_mock_aligned_read("AAAAAAAAAACCCCCCCCCC"));
//    
//    BOOST_CHECK(assembler.get_num_verticies() == 10);
//    BOOST_CHECK(assembler.get_num_edges() == 11);
//    
//    assembler.clear();
//    
//    BOOST_CHECK(assembler.get_num_verticies() == 0);
//    BOOST_CHECK(assembler.get_num_edges() == 0);
//    
//    assembler.add_read(get_mock_aligned_read("AAAAAAAAAACCCCCCCCCC"));
//    
//    BOOST_CHECK(assembler.get_num_verticies() == 10);
//    BOOST_CHECK(assembler.get_num_edges() == 11);
//    
//    assembler.add_read(get_mock_aligned_read("GGGGGGGGGGTTTTTTTTTT"));
//    
//    BOOST_CHECK(assembler.get_num_verticies() == 20);
//    BOOST_CHECK(assembler.get_num_edges() == 22);
//    
//    assembler.add_reference_contig(a_mock_region, "CCCCCCCCCCTTTTTTTTTT");
//    
//    BOOST_CHECK(assembler.get_num_verticies() == 28);
//    BOOST_CHECK(assembler.get_num_edges() == 31);
//    
//    assembler.add_reference_contig(a_mock_region, "GGGGGGGGGGAAAAAAAAAA");
//    
//    BOOST_CHECK(assembler.get_num_verticies() == 36);
//    BOOST_CHECK(assembler.get_num_edges() == 40);
//    
//    assembler.add_read(get_mock_aligned_read("CCCCCCCCCCGGGGGGGGGG"));
//    
//    BOOST_CHECK(assembler.get_num_verticies() == 44);
//    BOOST_CHECK(assembler.get_num_edges() == 51);
//    
//    assembler.add_reference_contig(a_mock_region, "CCCCCCCCCCGGGGGGGGGG");
//    
//    BOOST_CHECK(assembler.get_num_verticies() == 44);
//    BOOST_CHECK(assembler.get_num_edges() == 51);
//}

BOOST_AUTO_TEST_CASE(assembler_path_test)
{
    BOOST_REQUIRE(test_file_exists(lambda_reference_fasta));
    
    auto lambda = make_reference(lambda_reference_fasta);
    
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

//BOOST_AUTO_TEST_CASE(assembler_cycle_test)
//{
//    Assembler assembler {5};
//    
//    assembler.add_read("AAAAACCCCC");
//    
//    BOOST_CHECK(!assembler.has_cycle());
//    
//    assembler.add_read("CCCCCGGGGG");
//    
//    //BOOST_CHECK(!assembler.has_cycle());
//    
//    assembler.add_read("GGGGGAAAAA");
//    
//    //BOOST_CHECK(assembler.has_cycle());
//}

BOOST_AUTO_TEST_SUITE_END()
