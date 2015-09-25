//
//  haplotype_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 13/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>

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
#include "allele.hpp"
#include "haplotype_tree.hpp"

//BOOST_AUTO_TEST_CASE(haplotype hashing benchmark)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    auto a_region = parse_region("1:10000000-10001000", human);
//    
//    Haplotype haplotype {human, a_region};
//    
//    auto f_haplotype_hash = [&haplotype] () {
//        std::hash<Haplotype>()(haplotype);
//    };
//    
//    auto hash_time = benchmark<std::chrono::nanoseconds>(f_haplotype_hash, 10000).count();
//    
//    std::cout << "hash_time: " << hash_time << "ns" << std::endl;
//}

//BOOST_AUTO_TEST_CASE(haplotype containment benchmarks)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager {std::vector<std::string> {human_1000g_bam1}};
//    
//    auto a_region = parse_region("16:9300000-9300100", human);
//    
//    auto samples = a_read_manager.get_sample_ids();
//    
//    CandidateVariantGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//    
//    auto reads = a_read_manager.fetch_reads(samples[0], a_region);
//    
//    std::sort(reads.begin(), reads.end());
//    
//    candidate_generator.add_reads(reads.cbegin(), reads.cend());
//    
//    auto candidates = candidate_generator.get_candidates(a_region);
//    
//    HaplotypeTree haplotype_tree {human};
//    
//    for (const auto& candidate : candidates) {
//        haplotype_tree.extend(candidate.get_reference_allele());
//        haplotype_tree.extend(candidate.get_alternative_allele());
//    }
//    
//    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
//    
//    Allele allele1 {parse_region("16:9299945-9299946", human), "T"};
//    Allele allele2 {parse_region("16:9299946-9299957", human), "CGCATTACAAC"};
//    Allele allele3 {parse_region("16:9299957-9299958", human), "C"};
//    Allele allele4 {get_reference_allele(parse_region("16:9299958-9300037", human), human)};
//    Allele allele5 {parse_region("16:9300037-9300037", human), ""};
//    Allele allele6 {parse_region("16:9300037-9300039", human), "TG"};
//    Allele allele7 {parse_region("16:9300039-9300051", human), "TGTGTGTGCGTT"};
//    Allele allele8 {parse_region("16:9300051-9300061", human), "TGTGTGTGTG"};
//    Allele allele9 {parse_region("16:9300061-9300062", human), "G"};
//    Allele allele10 {parse_region("16:9300062-9300072", human), "GTGTGTGTGT"};
//    
//    std::vector<Allele> alleles {allele1, allele2, allele3, allele4, allele5, allele6, allele7, allele8, allele9, allele10};
//    
//    auto f_contains1 = [&haplotypes, &allele1] () {
//        for (const auto& haplotype : haplotypes) {
//            haplotype.contains(allele1);
//        }
//    };
//    auto f_contains2 = [&haplotypes, &allele2] () {
//        for (const auto& haplotype : haplotypes) {
//            haplotype.contains(allele2);
//        }
//    };
//    auto f_contains3 = [&haplotypes, &allele3] () {
//        for (const auto& haplotype : haplotypes) {
//            haplotype.contains(allele3);
//        }
//    };
//    auto f_contains4 = [&haplotypes, &allele4] () {
//        for (const auto& haplotype : haplotypes) {
//            haplotype.contains(allele4);
//        }
//    };
//    auto f_contains5 = [&haplotypes, &allele5] () {
//        for (const auto& haplotype : haplotypes) {
//            haplotype.contains(allele5);
//        }
//    };
//    auto f_contains6 = [&haplotypes, &allele6] () {
//        for (const auto& haplotype : haplotypes) {
//            haplotype.contains(allele6);
//        }
//    };
//    auto f_contains7 = [&haplotypes, &allele7] () {
//        for (const auto& haplotype : haplotypes) {
//            haplotype.contains(allele7);
//        }
//    };
//    auto f_contains8 = [&haplotypes, &allele8] () {
//        for (const auto& haplotype : haplotypes) {
//            haplotype.contains(allele8);
//        }
//    };
//    auto f_contains9 = [&haplotypes, &allele9] () {
//        for (const auto& haplotype : haplotypes) {
//            haplotype.contains(allele9);
//        }
//    };
//    auto f_contains10 = [&haplotypes, &allele10] () {
//        for (const auto& haplotype : haplotypes) {
//            haplotype.contains(allele10);
//        }
//    };
//    
//    auto time1  = benchmark<std::chrono::microseconds>(f_contains1, 10).count();
//    auto time2  = benchmark<std::chrono::microseconds>(f_contains2, 10).count();
//    auto time3  = benchmark<std::chrono::microseconds>(f_contains3, 10).count();
//    auto time4  = benchmark<std::chrono::microseconds>(f_contains4, 10).count();
//    auto time5  = benchmark<std::chrono::microseconds>(f_contains5, 10).count();
//    auto time6  = benchmark<std::chrono::microseconds>(f_contains6, 10).count();
//    auto time7  = benchmark<std::chrono::microseconds>(f_contains7, 10).count();
//    auto time8  = benchmark<std::chrono::microseconds>(f_contains8, 10).count();
//    auto time9  = benchmark<std::chrono::microseconds>(f_contains9, 10).count();
//    auto time10 = benchmark<std::chrono::microseconds>(f_contains10, 10).count();
//    
////    std::cout << time1 << std::endl;
////    std::cout << time2 << std::endl;
////    std::cout << time3 << std::endl;
////    std::cout << time4 << std::endl;
////    std::cout << time5 << std::endl;
////    std::cout << time6 << std::endl;
////    std::cout << time7 << std::endl;
////    std::cout << time8 << std::endl;
////    std::cout << time9 << std::endl;
////    std::cout << time10 << std::endl;
//}
