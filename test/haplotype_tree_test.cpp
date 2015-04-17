//
//  haplotype_tree_test.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <cstddef>
#include <algorithm>

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "read_manager.h"
#include "variant.h"
#include "variant_candidate_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "haplotype_tree.h"
#include "region_utils.h"
#include "read_filter.h"
#include "read_filters.h"

using std::cout;
using std::endl;

TEST_CASE("haplotype tree does not bifurcate on non-overlapping alleles", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000001-1000002", human), "C"};
    Allele allele3 {parse_region("4:1000002-1000002", human), "GC"};
    Allele allele4 {parse_region("4:1000005-1000007", human), ""};
    Allele allele5 {parse_region("4:1000007-1000008", human), "G"};
    
    HaplotypeTree haplotype_tree {human};
    haplotype_tree.extend_haplotypes(allele1);
    haplotype_tree.extend_haplotypes(allele2);
    haplotype_tree.extend_haplotypes(allele3);
    haplotype_tree.extend_haplotypes(allele4);
    haplotype_tree.extend_haplotypes(allele5);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 1);
}

TEST_CASE("haplotype tree ignores duplicate alleles coming from same allele", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000000-1000001", human), "C"};
    Allele allele3 {parse_region("4:1000000-1000001", human), "A"};
    
    HaplotypeTree haplotype_tree {human};
    haplotype_tree.extend_haplotypes(allele1);
    haplotype_tree.extend_haplotypes(allele2);
    haplotype_tree.extend_haplotypes(allele3);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
}

TEST_CASE("haplotype tree splits overlapping snps into different branches", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000000-1000001", human), "C"};
    Allele allele3 {parse_region("4:1000000-1000001", human), "G"};
    
    Allele allele4 {parse_region("4:1000001-1000002", human), "G"};
    Allele allele5 {parse_region("4:1000001-1000002", human), "C"};
    
    HaplotypeTree haplotype_tree {human};
    
    haplotype_tree.extend_haplotypes(allele1);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 1);
    
    haplotype_tree.extend_haplotypes(allele2);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    haplotype_tree.extend_haplotypes(allele3);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.extend_haplotypes(allele4);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.extend_haplotypes(allele5);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 6);
}

TEST_CASE("haplotype tree can selectively extend branches", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    HaplotypeTree haplotype_tree {human};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000000-1000003", human), ""};
    
    haplotype_tree.extend_haplotypes(allele1);
    haplotype_tree.extend_haplotypes(allele2);
    
    Allele allele3 {parse_region("4:1000001-1000002", human), "C"};
    Allele allele4 {parse_region("4:1000002-1000003", human), "G"};
    
    haplotype_tree.extend_haplotypes(allele2);
    haplotype_tree.extend_haplotypes(allele3);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    Allele allele5 {parse_region("4:1000003-1000004", human), "T"};
    
    haplotype_tree.extend_haplotypes(allele5);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    Allele allele6 {parse_region("4:1000003-1000004", human), "A"};
    
    haplotype_tree.extend_haplotypes(allele6);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 4);
}

TEST_CASE("haplotype tree can generate haplotypes in a region", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    HaplotypeTree haplotype_tree {human};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000002-1000003", human), "C"};
    Allele allele3 {parse_region("4:1000002-1000003", human), "G"};
    Allele allele4 {parse_region("4:1000004-1000005", human), "T"};
    
    haplotype_tree.extend_haplotypes(allele1);
    haplotype_tree.extend_haplotypes(allele2);
    haplotype_tree.extend_haplotypes(allele3);
    haplotype_tree.extend_haplotypes(allele4);
    
    auto a_region = parse_region("4:1000000-1000005", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    REQUIRE(haplotypes.size() == 2);
    std::sort(haplotypes.begin(), haplotypes.end());
    REQUIRE(haplotypes[0].get_sequence() == "ATCCT");
    REQUIRE(haplotypes[1].get_sequence() == "ATGCT");
}

TEST_CASE("haplotype tree can generate haplotypes ending in different regions", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    HaplotypeTree haplotype_tree {human};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000002-1000006", human), ""};
    Allele allele3 {parse_region("4:1000002-1000003", human), "G"};
    
    haplotype_tree.extend_haplotypes(allele1);
    haplotype_tree.extend_haplotypes(allele2);
    haplotype_tree.extend_haplotypes(allele3);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    auto a_region = parse_region("4:1000000-1000006", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    REQUIRE(haplotypes.size() == 2);
    
    a_region = parse_region("4:1000000-1000003", human);
    
    haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    REQUIRE(haplotypes.size() == 2);
    std::sort(haplotypes.begin(), haplotypes.end());
    REQUIRE(haplotypes[0].get_sequence() == "ATG");
    REQUIRE(haplotypes[1].get_sequence() == "AT");
}

TEST_CASE("leading haplotypes can be removed from the tree", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    HaplotypeTree haplotype_tree {human};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000002-1000003", human), "C"};
    Allele allele3 {parse_region("4:1000002-1000003", human), "G"};
    Allele allele4 {parse_region("4:1000004-1000005", human), "T"};
    Allele allele5 {parse_region("4:1000004-1000005", human), "C"};
    
    haplotype_tree.extend_haplotypes(allele1);
    haplotype_tree.extend_haplotypes(allele2);
    haplotype_tree.extend_haplotypes(allele3);
    haplotype_tree.extend_haplotypes(allele4);
    haplotype_tree.extend_haplotypes(allele5);
    
    auto a_region = parse_region("4:1000000-1000005", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    REQUIRE(haplotypes.size() == 4);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    haplotype_tree.prune_haplotype(haplotypes[0]); // ATCCC
    
    REQUIRE(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.prune_haplotype(haplotypes[1]); // ATCCT
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    REQUIRE(haplotypes[0].get_sequence() == "ATGCC");
    REQUIRE(haplotypes[1].get_sequence() == "ATGCT");
}

TEST_CASE("pruned branches can still be extended", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    HaplotypeTree haplotype_tree {human};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000002-1000003", human), "C"};
    Allele allele3 {parse_region("4:1000002-1000003", human), "G"};
    Allele allele4 {parse_region("4:1000004-1000005", human), "T"};
    Allele allele5 {parse_region("4:1000004-1000005", human), "C"};
    Allele allele6 {parse_region("4:1000006-1000007", human), "A"};
    
    haplotype_tree.extend_haplotypes(allele1);
    haplotype_tree.extend_haplotypes(allele2);
    haplotype_tree.extend_haplotypes(allele3);
    haplotype_tree.extend_haplotypes(allele4);
    haplotype_tree.extend_haplotypes(allele5);
    haplotype_tree.extend_haplotypes(allele6);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 4);
    
    Haplotype hap1 {human};
    hap1.push_back(allele3);
    hap1.push_back(allele4);
    hap1.push_back(allele6);
    
    Haplotype hap2 {human};
    hap2.push_back(allele5);
    hap2.push_back(allele6);
    
    haplotype_tree.prune_haplotype(hap1);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 3);
    
    auto a_region = parse_region("4:1000000-1000007", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    REQUIRE(haplotypes[0].get_sequence() == "ATCCCAA");
    REQUIRE(haplotypes[1].get_sequence() == "ATCCTAA");
    REQUIRE(haplotypes[2].get_sequence() == "ATGCCAA");
    
    haplotype_tree.prune_haplotype(hap2);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    REQUIRE(haplotypes[0].get_sequence() == "ATG");
    REQUIRE(haplotypes[1].get_sequence() == "ATCCTAA");
    
    haplotype_tree.extend_haplotypes(allele4);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    REQUIRE(haplotypes[0].get_sequence() == "ATGCT");
    REQUIRE(haplotypes[1].get_sequence() == "ATCCTAA");
}

//TEST_CASE("haplotype_tree_single_sample_test", "[haplotype_tree]")
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human(a_factory.make(human_reference_fasta));
//    
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
//    
//    VariantCandidateGenerator candidate_generator {};
//    
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//    
//    auto sample_ids = a_read_manager.get_sample_ids();
//    auto the_sample_id = sample_ids.at(0);
//    
//    //auto a_region = parse_region("12:0-10000000", human);
//    auto a_region = parse_region("16:9300000-9300100", human);
//    
//    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
//    
//    using ReadIterator = std::vector<AlignedRead>::const_iterator;
//    ReadFilter<ReadIterator> a_read_filter {};
//    a_read_filter.register_filter(is_mapped);
//    
//    std::vector<AlignedRead> good_reads {}, bad_reads {};
//    good_reads.reserve(reads.size());
//    bad_reads.reserve(reads.size());
//    a_read_filter.filter_reads(std::make_move_iterator(reads.begin()), std::make_move_iterator(reads.end()),
//                               ContextBackInserter(good_reads), ContextBackInserter(bad_reads));
//    
//    candidate_generator.add_reads(good_reads.cbegin(), good_reads.cend());
//    
//    auto candidates = candidate_generator.get_candidates(a_region);
//    
//    HaplotypeTree haplotype_tree {human};
//    
//    std::vector<Variant> aligned_candidates {};
//    std::transform(candidates.begin(), candidates.end(),
//                      std::back_inserter(aligned_candidates),
//                      [&human] (const auto& a_variant) {
//                          return left_align(a_variant, human);
//                      });
//    
//    for (const auto& candidate : aligned_candidates) {
//        //cout << candidate << endl;
//        //haplotype_tree.extend_haplotypes(candidate);
//    }
//}
