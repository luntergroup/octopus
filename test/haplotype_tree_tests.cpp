//
//  haplotype_tree_tests.cpp
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
#include <set>

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "read_manager.h"
#include "variant.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "haplotype_tree.h"
#include "region_algorithms.h"
#include "read_filter.h"
#include "read_filters.h"

using std::cout;
using std::endl;

TEST_CASE("haplotype tree does not bifurcate on alleles positioned past the leading alleles", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000001-1000002", human), "C"};
    Allele allele3 {parse_region("4:1000002-1000002", human), "GC"};
    Allele allele4 {parse_region("4:1000005-1000007", human), ""};
    Allele allele5 {parse_region("4:1000007-1000008", human), "G"};
    
    HaplotypeTree haplotype_tree {human};
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    haplotype_tree.extend(allele4);
    haplotype_tree.extend(allele5);
    
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
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
}

TEST_CASE("haplotype tree ignores insertions followed immediatly by deletions and vice versa", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    Allele allele1 {parse_region("16:9300037-9300037", human), "TG"};
    Allele allele2 {parse_region("16:9300037-9300051 ", human), ""};
    
    HaplotypeTree haplotype_tree {human};
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 1);
    
    auto a_region = parse_region("16:9300037-9300037", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    REQUIRE(haplotypes[0].contains(allele1));
    REQUIRE(!haplotypes[0].contains(allele2));
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
    
    haplotype_tree.extend(allele1);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 1);
    
    haplotype_tree.extend(allele2);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    haplotype_tree.extend(allele3);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.extend(allele4);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.extend(allele5);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 6);
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
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    haplotype_tree.extend(allele4);
    
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
    ReferenceGenome human {a_factory.make(human_reference_fasta)};;
    
    HaplotypeTree haplotype_tree {human};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000002-1000006", human), ""};
    Allele allele3 {parse_region("4:1000002-1000003", human), "G"};
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    
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
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    haplotype_tree.extend(allele4);
    haplotype_tree.extend(allele5);
    
    auto a_region = parse_region("4:1000000-1000005", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    REQUIRE(haplotypes.size() == 4);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    haplotype_tree.prune_all(haplotypes[0]); // ATCCC
    
    REQUIRE(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.prune_all(haplotypes[1]); // ATCCT
    
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
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    haplotype_tree.extend(allele4);
    haplotype_tree.extend(allele5);
    haplotype_tree.extend(allele6);
    
    auto a_region = parse_region("4:1000000-1000007", human);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 4);
    
    Haplotype hap1 {human};
    hap1.push_back(allele3);
    hap1.push_back(allele4);
    hap1.push_back(allele6);
    
    haplotype_tree.prune_all(hap1);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 3);
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    REQUIRE(haplotypes[0].get_sequence() == "ATCCCAA");
    REQUIRE(haplotypes[1].get_sequence() == "ATCCTAA");
    REQUIRE(haplotypes[2].get_sequence() == "ATGCCAA");
    
    Haplotype hap2 {human};
    hap2.push_back(allele5);
    hap2.push_back(allele6);
    
    haplotype_tree.prune_all(hap2);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    REQUIRE(haplotypes[0].get_sequence() == "ATCCTAA");
    REQUIRE(haplotypes[1].get_sequence() == "ATGCCAG");
    
    haplotype_tree.extend(allele4);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    REQUIRE(haplotypes[0].get_sequence() == "ATCCTAA");
    REQUIRE(haplotypes[1].get_sequence() == "ATGCTAG");
}

TEST_CASE("extending on mnps results in backtracked bifurification", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    Allele allele1 {parse_region("16:9300039-9300051", human), "TGTGTGTGCGTT"};
    Allele allele2 {parse_region("16:9300039-9300051", human), ""};
    
    Allele allele3 {parse_region("16:9300047-9300048", human), "C"};
    Allele allele4 {parse_region("16:9300047-9300048", human), "T"};
    Allele allele5 {parse_region("16:9300050-9300051", human), "T"};
    Allele allele6 {parse_region("16:9300050-9300051", human), "G"};
    
    HaplotypeTree haplotype_tree {human};
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    haplotype_tree.extend(allele3);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.extend(allele4);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 4);
    
    haplotype_tree.extend(allele5);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 5);
    
    haplotype_tree.extend(allele6);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 8);
    
    auto a_region = parse_region("16:9300039-9300051", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    REQUIRE(haplotypes.size() == 8);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    auto it = std::unique(haplotypes.begin(), haplotypes.end());
    
    haplotypes.erase(it, haplotypes.end());
    
    REQUIRE(haplotypes.size() == 5);
}

TEST_CASE("haplotype tree can selectively extend branches", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    HaplotypeTree haplotype_tree {human};
    
    Allele allele1 {parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {parse_region("4:1000000-1000003", human), ""};
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 2);
    
    Allele allele3 {parse_region("4:1000001-1000002", human), "C"};
    Allele allele4 {parse_region("4:1000002-1000003", human), "G"};
    
    haplotype_tree.extend(allele3);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.extend(allele4);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 4);
    
    Allele allele5 {parse_region("4:1000003-1000004", human), "T"};
    
    haplotype_tree.extend(allele5);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 4);
    
    Allele allele6 {parse_region("4:1000003-1000004", human), "A"};
    
    haplotype_tree.extend(allele6);
    
    REQUIRE(haplotype_tree.num_haplotypes() == 8);
}

TEST_CASE("haplotype tree survives serious pruning", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    CandidateVariantGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto a_region = parse_region("16:9299940-9300055", human);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    HaplotypeTree haplotype_tree {human};
    
    for (const auto& candidate : candidates) {
        haplotype_tree.extend(candidate.get_reference_allele());
        haplotype_tree.extend(candidate.get_alternative_allele());
    }
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    Haplotype the_reference_haplotype {human, a_region};
    
    auto er = std::equal_range(haplotypes.begin(), haplotypes.end(), the_reference_haplotype);
    
    haplotypes.erase(er.first, er.second);
    
    for (const auto& haplotype : haplotypes) {
        haplotype_tree.prune_all(haplotype);
    }
    
    REQUIRE(haplotype_tree.num_haplotypes() > 0);
    
    auto pruned_haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    for (const auto& pruned_haplotype : pruned_haplotypes) {
        REQUIRE(pruned_haplotype == the_reference_haplotype);
    }
    
    unique_least_complex(pruned_haplotypes);
    
    auto it = pruned_haplotypes.cbegin() + 1;
    
    std::for_each(it, pruned_haplotypes.cend(), [&haplotype_tree] (const auto& haplotype) {
        haplotype_tree.prune_all(haplotype);
    });
    
    REQUIRE(haplotype_tree.num_haplotypes() == 5);
    
    auto remaining_haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    REQUIRE(std::all_of(remaining_haplotypes.cbegin(), remaining_haplotypes.cend(),
                        [&pruned_haplotypes] (const auto& haplotype) {
                            return haplotype == pruned_haplotypes.front();
                        }));
    
    haplotype_tree.prune_all(remaining_haplotypes.front());
    
    REQUIRE(haplotype_tree.empty());
}

TEST_CASE("prune_unqiue leaves a single haplotype which contains the same alleles as the given haplotype", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    CandidateVariantGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto a_region = parse_region("16:9299940-9300055", human);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    HaplotypeTree haplotype_tree {human};
    
    for (const auto& candidate : candidates) {
        haplotype_tree.extend(candidate.get_reference_allele());
        haplotype_tree.extend(candidate.get_alternative_allele());
    }
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    Haplotype the_reference_haplotype {human, a_region};
    
    std::sort(haplotypes.begin(), haplotypes.end());
    
    auto er = std::equal_range(haplotypes.begin(), haplotypes.end(), the_reference_haplotype);
    
    std::sort(er.first, er.second, is_less_complex);
    
    REQUIRE(std::distance(er.first, er.second) == 5);
    
    auto haplotype_to_prune= *er.first;
    
    haplotype_tree.prune_unique(haplotype_to_prune);
    
    auto new_haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    std::sort(new_haplotypes.begin(), new_haplotypes.end());
    
    er = std::equal_range(new_haplotypes.begin(), new_haplotypes.end(), the_reference_haplotype);
    
    REQUIRE(std::distance(er.first, er.second) == 1);
}
