//
//  haplotype_tree_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <cstddef>
#include <algorithm>
#include <set>

#include "test_common.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "variant.hpp"
#include "candidate_variant_generator.hpp"
#include "alignment_candidate_variant_generator.hpp"
#include "haplotype.hpp"
#include "haplotype_tree.hpp"
#include "mappable_algorithms.hpp"
#include "read_filter.hpp"
#include "read_filters.hpp"

using std::cout;
using std::endl;

using Octopus::HaplotypeTree;

BOOST_AUTO_TEST_SUITE(Components)

static void sort(std::vector<Haplotype>& haplotypes)
{
    std::sort(std::begin(haplotypes), std::end(haplotypes));
}

BOOST_AUTO_TEST_CASE(haplotype_tree_does_not_bifurcate_on_alleles_positioned_past_the_leading_alleles)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    Allele allele1 {*parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {*parse_region("4:1000001-1000002", human), "C"};
    Allele allele3 {*parse_region("4:1000002-1000002", human), "GC"};
    Allele allele4 {*parse_region("4:1000005-1000007", human), ""};
    Allele allele5 {*parse_region("4:1000007-1000008", human), "G"};
    
    HaplotypeTree haplotype_tree {"4", human};
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    haplotype_tree.extend(allele4);
    haplotype_tree.extend(allele5);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 1);
}

BOOST_AUTO_TEST_CASE(haplotype_tree_ignores_duplicate_alleles_coming_from_same_allele)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    Allele allele1 {*parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {*parse_region("4:1000000-1000001", human), "C"};
    Allele allele3 {*parse_region("4:1000000-1000001", human), "A"};
    
    HaplotypeTree haplotype_tree {"4", human};
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
}

BOOST_AUTO_TEST_CASE(haplotype_tree_ignores_insertions_followed_immediatly_by_deletions_and_vice_versa)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    Allele allele1 {*parse_region("16:9300037-9300037", human), "TG"};
    Allele allele2 {*parse_region("16:9300037-9300051 ", human), ""};
    
    HaplotypeTree haplotype_tree {"16", human};
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 1);
    
    auto region = *parse_region("16:9300037-9300037", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(region);
    
    BOOST_CHECK(haplotypes[0].contains(allele1));
    BOOST_CHECK(!haplotypes[0].contains(allele2));
}

BOOST_AUTO_TEST_CASE(haplotype_tree_splits_overlapping_snps_into_different_branches)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    Allele allele1 {*parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {*parse_region("4:1000000-1000001", human), "C"};
    Allele allele3 {*parse_region("4:1000000-1000001", human), "G"};
    
    Allele allele4 {*parse_region("4:1000001-1000002", human), "G"};
    Allele allele5 {*parse_region("4:1000001-1000002", human), "C"};
    
    HaplotypeTree haplotype_tree {"4", human};
    
    haplotype_tree.extend(allele1);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 1);
    
    haplotype_tree.extend(allele2);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
    
    haplotype_tree.extend(allele3);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.extend(allele4);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.extend(allele5);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 6);
}

BOOST_AUTO_TEST_CASE(haplotype_tree_can_generate_haplotypes_in_a_region)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    Allele allele1 {*parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {*parse_region("4:1000002-1000003", human), "C"};
    Allele allele3 {*parse_region("4:1000002-1000003", human), "G"};
    Allele allele4 {*parse_region("4:1000004-1000005", human), "T"};
    
    HaplotypeTree haplotype_tree {"4", human};
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    haplotype_tree.extend(allele4);
    
    auto region = *parse_region("4:1000000-1000005", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(region);
    
    BOOST_CHECK(haplotypes.size() == 2);
    
    sort(haplotypes);
    
    BOOST_CHECK(haplotypes[0].get_sequence() == "ATCCT");
    BOOST_CHECK(haplotypes[1].get_sequence() == "ATGCT");
}

BOOST_AUTO_TEST_CASE(haplotype_tree_can_generate_haplotypes_ending_in_different_regions)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    HaplotypeTree haplotype_tree {"4", human};
    
    Allele allele1 {*parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {*parse_region("4:1000002-1000006", human), ""};
    Allele allele3 {*parse_region("4:1000002-1000003", human), "G"};
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
    
    const auto region1 = *parse_region("4:1000000-1000006", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(region1);
    
    BOOST_CHECK(haplotypes.size() == 2);
    
    const auto region2 = *parse_region("4:1000000-1000003", human);
    
    haplotypes = haplotype_tree.get_haplotypes(region2);
    
    BOOST_CHECK(haplotypes.size() == 2);
    
    sort(haplotypes);
    
    BOOST_CHECK(haplotypes[0].get_sequence() == "ATG");
    BOOST_CHECK(haplotypes[1].get_sequence() == "AT");
}

BOOST_AUTO_TEST_CASE(leading_haplotypes_can_be_removed_from_the_tree)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    HaplotypeTree haplotype_tree {"4", human};
    
    Allele allele1 {*parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {*parse_region("4:1000002-1000003", human), "C"};
    Allele allele3 {*parse_region("4:1000002-1000003", human), "G"};
    Allele allele4 {*parse_region("4:1000004-1000005", human), "T"};
    Allele allele5 {*parse_region("4:1000004-1000005", human), "C"};
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    haplotype_tree.extend(allele4);
    haplotype_tree.extend(allele5);
    
    const auto region = get_encompassing(allele1, allele5);
    
    auto haplotypes = haplotype_tree.get_haplotypes(region);
    
    BOOST_CHECK(haplotypes.size() == 4);
    
    sort(haplotypes);
    
    haplotype_tree.prune_all(haplotypes[0]); // ATCCC
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.prune_all(haplotypes[1]); // ATCCT
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
    
    haplotypes = haplotype_tree.get_haplotypes(region);
    
    sort(haplotypes);
    
    BOOST_CHECK(haplotypes[0].get_sequence() == "ATGCC");
    BOOST_CHECK(haplotypes[1].get_sequence() == "ATGCT");
}

BOOST_AUTO_TEST_CASE(pruned_branches_can_still_be_extended)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    HaplotypeTree haplotype_tree {"4", human};
    
    Allele allele1 {*parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {*parse_region("4:1000002-1000003", human), "C"};
    Allele allele3 {*parse_region("4:1000002-1000003", human), "G"};
    Allele allele4 {*parse_region("4:1000004-1000005", human), "T"};
    Allele allele5 {*parse_region("4:1000004-1000005", human), "C"};
    Allele allele6 {*parse_region("4:1000006-1000007", human), "A"};
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    haplotype_tree.extend(allele4);
    haplotype_tree.extend(allele5);
    haplotype_tree.extend(allele6);
    
    const auto region = *parse_region("4:1000000-1000007", human);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 4);
    
    Haplotype hap1 {region, human};
    hap1.push_back(allele3);
    hap1.push_back(allele4);
    hap1.push_back(allele6);
    
    haplotype_tree.prune_all(hap1);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
    
    auto haplotypes = haplotype_tree.get_haplotypes(region);
    
    sort(haplotypes);
    
    BOOST_CHECK(haplotypes[0].get_sequence() == "ATCCCAA");
    BOOST_CHECK(haplotypes[1].get_sequence() == "ATCCTAA");
    BOOST_CHECK(haplotypes[2].get_sequence() == "ATGCCAA");
    
    Haplotype hap2 {region, human};
    hap2.push_back(allele5);
    hap2.push_back(allele6);
    
    haplotype_tree.prune_all(hap2);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
    
    haplotypes = haplotype_tree.get_haplotypes(region);
    
    sort(haplotypes);
    
    BOOST_CHECK(haplotypes[0].get_sequence() == "ATCCTAA");
    BOOST_CHECK(haplotypes[1].get_sequence() == "ATGCCAG");
    
    haplotype_tree.extend(allele4);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
    
    haplotypes = haplotype_tree.get_haplotypes(region);
    
    sort(haplotypes);
    
    BOOST_CHECK(haplotypes[0].get_sequence() == "ATCCTAA");
    BOOST_CHECK(haplotypes[1].get_sequence() == "ATGCTAG");
}

BOOST_AUTO_TEST_CASE(extending_on_mnps_results_in_backtracked_bifurification)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    Allele allele1 {*parse_region("16:9300039-9300051", human), "TGTGTGTGCGTT"};
    Allele allele2 {*parse_region("16:9300039-9300051", human), ""};
    
    Allele allele3 {*parse_region("16:9300047-9300048", human), "C"};
    Allele allele4 {*parse_region("16:9300047-9300048", human), "T"};
    Allele allele5 {*parse_region("16:9300050-9300051", human), "T"};
    Allele allele6 {*parse_region("16:9300050-9300051", human), "G"};
    
    HaplotypeTree haplotype_tree {"16", human};
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
    
    haplotype_tree.extend(allele3);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.extend(allele4);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 4);
    
    haplotype_tree.extend(allele5);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 5);
    
    haplotype_tree.extend(allele6);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 8);
    
    auto a_region = *parse_region("16:9300039-9300051", human);
    
    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
    
    BOOST_CHECK(haplotypes.size() == 8);
    
    sort(haplotypes);
    
    auto it = std::unique(haplotypes.begin(), haplotypes.end());
    
    haplotypes.erase(it, haplotypes.end());
    
    BOOST_CHECK(haplotypes.size() == 5);
}

BOOST_AUTO_TEST_CASE(haplotype_tree_can_selectively_extend_branches)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    HaplotypeTree haplotype_tree {"4", human};
    
    Allele allele1 {*parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {*parse_region("4:1000000-1000003", human), ""};
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 2);
    
    Allele allele3 {*parse_region("4:1000001-1000002", human), "C"};
    Allele allele4 {*parse_region("4:1000002-1000003", human), "G"};
    
    haplotype_tree.extend(allele3);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 3);
    
    haplotype_tree.extend(allele4);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 4);
    
    Allele allele5 {*parse_region("4:1000003-1000004", human), "T"};
    
    haplotype_tree.extend(allele5);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 4);
    
    Allele allele6 {*parse_region("4:1000003-1000004", human), "A"};
    
    haplotype_tree.extend(allele6);
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 8);
}

BOOST_AUTO_TEST_CASE(haplotype_tree_survives_serious_pruning)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    ReadManager read_manager {NA12878_low_coverage};
    
    Octopus::CandidateVariantGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<Octopus::AlignmentCandidateVariantGenerator>(human, 0));
    
    const auto sample = read_manager.get_samples().front();
    
    const auto region = *parse_region("16:9299940-9300055", human);
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(region);
    
    HaplotypeTree haplotype_tree {"16", human};
    
    for (const auto& candidate : candidates) {
        haplotype_tree.extend(candidate.get_ref_allele());
        haplotype_tree.extend(candidate.get_alt_allele());
    }
    
    auto haplotypes = haplotype_tree.get_haplotypes(region);
    
    sort(haplotypes);
    
    Haplotype the_reference_haplotype {region, human};
    
    auto er = std::equal_range(haplotypes.begin(), haplotypes.end(), the_reference_haplotype);
    
    haplotypes.erase(er.first, er.second);
    
    for (const auto& haplotype : haplotypes) {
        haplotype_tree.prune_all(haplotype);
    }
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() > 0);
    
    auto pruned_haplotypes = haplotype_tree.get_haplotypes(region);
    
    for (const auto& pruned_haplotype : pruned_haplotypes) {
        BOOST_CHECK(pruned_haplotype == the_reference_haplotype);
    }
    
    unique_least_complex(pruned_haplotypes);
    
    auto it = pruned_haplotypes.cbegin() + 1;
    
    std::for_each(it, pruned_haplotypes.cend(), [&haplotype_tree] (const auto& haplotype) {
        haplotype_tree.prune_all(haplotype);
    });
    
    BOOST_CHECK(haplotype_tree.num_haplotypes() == 5);
    
    auto remaining_haplotypes = haplotype_tree.get_haplotypes(region);
    
    BOOST_CHECK(std::all_of(remaining_haplotypes.cbegin(), remaining_haplotypes.cend(),
                        [&pruned_haplotypes] (const auto& haplotype) {
                            return haplotype == pruned_haplotypes.front();
                        }));
    
    haplotype_tree.prune_all(remaining_haplotypes.front());
    
    BOOST_CHECK(haplotype_tree.empty());
}

BOOST_AUTO_TEST_CASE(prune_unqiue_leaves_a_single_haplotype_which_contains_the_same_alleles_as_the_given_haplotype)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    ReadManager read_manager {NA12878_low_coverage};
    
    const auto region = *parse_region("16:9299940-9300055", human);
    
    const auto sample = read_manager.get_samples().front();
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    Octopus::CandidateVariantGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<Octopus::AlignmentCandidateVariantGenerator>(human, 0));
    
    add_reads(reads, candidate_generator);
    
    auto candidates = candidate_generator.get_candidates(region);
    
    HaplotypeTree haplotype_tree {"16", human};
    
    for (const auto& candidate : candidates) {
        haplotype_tree.extend(candidate.get_ref_allele());
        haplotype_tree.extend(candidate.get_alt_allele());
    }
    
    auto haplotypes = haplotype_tree.get_haplotypes(region);
    
    Haplotype the_reference_haplotype {region, human};
    
    sort(haplotypes);
    
    auto er = std::equal_range(haplotypes.begin(), haplotypes.end(), the_reference_haplotype);
    
    std::sort(er.first, er.second, IsLessComplex());
    
    BOOST_CHECK(std::distance(er.first, er.second) == 5);
    
    auto haplotype_to_prune= *er.first;
    
    haplotype_tree.prune_unique(haplotype_to_prune);
    
    auto new_haplotypes = haplotype_tree.get_haplotypes(region);
    
    sort(new_haplotypes);
    
    er = std::equal_range(new_haplotypes.begin(), new_haplotypes.end(), the_reference_haplotype);
    
    BOOST_CHECK(std::distance(er.first, er.second) == 1);
}

BOOST_AUTO_TEST_CASE(contains_returns_true_if_the_given_haplotype_is_in_the_tree_in_any_form)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    Allele allele1 {*parse_region("4:1000000-1000001", human), "A"};
    Allele allele2 {*parse_region("4:1000000-1000001", human), "C"};
    Allele allele3 {*parse_region("4:1000000-1000001", human), "G"};
    Allele allele4 {*parse_region("4:1000001-1000002", human), "G"};
    Allele allele5 {*parse_region("4:1000001-1000002", human), "C"};
    
    const auto region = get_encompassing(allele1, allele5);
    
    HaplotypeTree haplotype_tree {"4", human};
    
    haplotype_tree.extend(allele1);
    haplotype_tree.extend(allele2);
    haplotype_tree.extend(allele3);
    haplotype_tree.extend(allele4);
    haplotype_tree.extend(allele5);
    
    Haplotype haplotype1 {region, human};
    haplotype1.push_back(allele1);
    haplotype1.push_back(allele4);
    
    BOOST_CHECK(haplotype_tree.contains(haplotype1));
    
    Haplotype haplotype2 {region, human};
    haplotype2.push_back(allele1);
    haplotype2.push_back(Allele {*parse_region("4:1000001-1000002", human), "A"});
    
    BOOST_CHECK(!haplotype_tree.contains(haplotype2));
}

BOOST_AUTO_TEST_CASE(is_unique_return_true_if_the_given_haplotype_occurs_extactly_once_in_the_tree)
{
    // TODO
}

BOOST_AUTO_TEST_SUITE_END()
