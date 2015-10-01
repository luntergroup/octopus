//
//  read_model_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <cstddef>
#include <unordered_map>

#include "test_common.hpp"
#include "benchmark_utils.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "variant.hpp"
#include "variant_utils.hpp"
#include "candidate_variant_generator.hpp"
#include "alignment_candidate_variant_generator.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "read_model.hpp"

//BOOST_AUTO_TEST_CASE(read model benchmarks)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    auto a_region = parse_region("1:10000000-10001000", human);
//    
//    ReadManager a_read_manager {std::vector<std::string> {human_1000g_bam1}};
//    
//    auto sample_id = a_read_manager.get_sample_ids().at(0);
//    
//    auto reads = a_read_manager.fetch_reads(sample_id, a_region);
//    
//    unsigned ploidy {2};
//    ReadModel read_model {ploidy};
//    
//    Haplotype reference_haplotype {human, a_region};
//    
//    std::vector<Haplotype> haplotypes {reference_haplotype};
//    
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//    
//    auto f_genotype_prob = [&sample_id, &read_model, &reads, &genotypes] () {
//        read_model.log_probability(reads.cbegin(), reads.cend(), genotypes[0], sample_id);
//    };
//    
//    auto time = benchmark<std::chrono::microseconds>(f_genotype_prob, 1).count();
//    
//    std::cout << "time to align haplotype not in cache " << time << std::endl;
//    
//    time = benchmark<std::chrono::microseconds>(f_genotype_prob, 1000).count();
//    
//    std::cout << "time to align haplotype in cache " << time << std::endl;
//}
