//
//  reference_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <future>
#include <random>

#include "test_common.hpp"
#include "reference_genome.hpp"
#include "mappable_algorithms.hpp"

#include "fasta.hpp"
#include "caching_fasta.hpp"

using std::cout;
using std::endl;

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(Components)
BOOST_AUTO_TEST_SUITE(IO)

BOOST_AUTO_TEST_CASE(ReferemceGenome_handles_basic_queries)
{
    BOOST_REQUIRE(test_file_exists(ecoli_reference_fasta));
    
    const auto ecoli = make_reference(ecoli_reference_fasta);
    
    BOOST_CHECK(ecoli.name() == "R00000042");
    BOOST_CHECK(ecoli.contains_region(GenomicRegion("R00000042", 10000, 2000000)));
    BOOST_CHECK(ecoli.contig_size("R00000042") == 5231428);
    BOOST_CHECK(!ecoli.has_contig("X"));
    BOOST_CHECK(ecoli.contig_region("R00000042") == GenomicRegion("R00000042", 0, 5231428));
    BOOST_CHECK(ecoli.fetch_sequence(GenomicRegion("R00000042", 0, 10)) == "AGCTTTTCAT"); // first line
    BOOST_CHECK(ecoli.fetch_sequence(GenomicRegion("R00000042", 69, 80)) == "CTTCTGAACTG"); // accross lines
    
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    BOOST_CHECK(human.name() == "human_g1k_v37");
    BOOST_CHECK(human.contains_region(GenomicRegion("1", 100, 10000)));
    BOOST_CHECK(!human.contains_region(GenomicRegion("1", 100, 3e8))); // too big
    BOOST_CHECK(human.contig_size("20") == 63025520);
    BOOST_CHECK(human.has_contig("X"));
    BOOST_CHECK(!human.has_contig("y"));
    BOOST_CHECK(human.contig_region("X") == GenomicRegion("X", 0, 155270560));
    BOOST_CHECK(human.fetch_sequence(GenomicRegion("15", 51265690, 51265700)) == "ACAATGTTGT");
    BOOST_CHECK(human.fetch_sequence(GenomicRegion("5", 100000, 100010)) == "AGGAAGTTTC");
}

BOOST_AUTO_TEST_CASE(ReferemceGenome_handles_edge_cases)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    BOOST_CHECK(human.fetch_sequence(GenomicRegion {"1", 100, 100}) == "");
}

BOOST_AUTO_TEST_CASE(parse_region_works_with_correctly_formatted_region_input)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    const auto r1 = parse_region("3", human);
    
    BOOST_CHECK(r1.contig_name() == "3");
    BOOST_CHECK(r1.begin() == 0);
    BOOST_CHECK(r1.end() == human.contig_size("3"));
    
    const auto r2 = parse_region("10:100-200", human);
    
    BOOST_CHECK(r2.contig_name() == "10");
    BOOST_CHECK(r2.begin() == 100);
    BOOST_CHECK(r2.end() == 200);
    
    const auto r3 = parse_region("18:102,029-102,029", human);
    
    BOOST_CHECK(r3.contig_name() == "18");
    BOOST_CHECK(r3.begin() == 102'029);
    BOOST_CHECK(r3.end() == 102'029);
    
    const auto r4 = parse_region("MT:100-", human);
    
    BOOST_CHECK(r4.contig_name() == "MT");
    BOOST_CHECK(r4.begin() == 100);
    BOOST_CHECK(r4.end() == human.contig_size("MT"));
    
    const auto r5 = parse_region("7:1,000,000", human);
    
    BOOST_CHECK(r5.contig_name() == "7");
    BOOST_CHECK(r5.begin() == 1'000'000);
    BOOST_CHECK(r5.end() == 1'000'001);
}

BOOST_AUTO_TEST_CASE(parse_region_throws_when_given_bad_region)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    bool all_throwed {true};
    
    try {
        parse_region("", human);
        all_throwed = false;
    } catch (...) {}
    try {
        parse_region("-", human);
        all_throwed = false;
    } catch (...) {}
    try {
        parse_region("5:100-99", human);
        all_throwed = false;
    } catch (...) {}
    try {
        parse_region("not_in_human", human);
        all_throwed = false;
    } catch (...) {}
    try {
        parse_region("0", human);
        all_throwed = false;
    } catch (...) {}
    try {
        parse_region("-1", human);
        all_throwed = false;
    } catch (...) {}
    try {
        parse_region("--1", human);
        all_throwed = false;
    } catch (...) {}
    try {
        parse_region("1:", human);
        all_throwed = false;
    } catch (...) {}
    try {
        parse_region("1:-", human);
        all_throwed = false;
    } catch (...) {}
    try {
        parse_region("2::0-100", human);
        all_throwed = false;
    } catch (...) {}
    
    
    BOOST_CHECK(all_throwed);
}

//BOOST_AUTO_TEST_CASE(cached_and_uncached_reference_genome_give_same_sequence)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    // First test simple walk along a single contig
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    std::string contig {"1"};
//    
//    std::random_device rd;
//    std::mt19937 g(rd());
//    
//    std::uniform_int_distribution<GenomicRegion::SizeType> dis(1, 50000);
//    
//    Fasta uncached_reference {human_reference_fasta};
//    CachingFasta cached_reference {human_reference_fasta, 10'000'000, 1.0, 1.0};
//    
//    std::string uncached_sequence {}, cached_sequence {};
//    uncached_sequence.reserve(uncached_reference.contig_size(contig));
//    cached_sequence.reserve(cached_reference.contig_size(contig));
//    
//    auto regions = decompose(parse_region(contig, human), 100);
//    
//    for (const auto& region : regions) {
//        uncached_sequence += uncached_reference.fetch_sequence(region);
//        cached_sequence += cached_reference.fetch_sequence(region);
//    }
//    
//    BOOST_CHECK(uncached_sequence == cached_sequence);
//    
////    // Next test walks along multiple contigs
////    
////    regions.clear();
////    
////    for (const auto& c : std::vector<std::string> {"1", "2", "3"}) {
////        auto more_regions = decompose(uncached_reference.contig_region(c), dis(g));
////        regions.insert(regions.end(), more_regions.begin(), more_regions.end());
////    }
////    
////    
//    
////    // Next test accesses to random points on a single contig
////    
////    uncached_sequence.clear();
////    cached_sequence.clear();
////    
////    std::shuffle(regions.begin(), regions.end(), g);
////    
////    for (const auto& region : regions) {
////        uncached_sequence += uncached_reference.fetch_sequence(region);
////        cached_sequence += cached_reference.fetch_sequence(region);
////    }
////    
////    BOOST_CHECK(uncached_sequence == cached_sequence);
////    
////    // Finally test random access to lots of contigs with different region sizes - may take a while
////    
////    //uncached_sequence.reserve(uncached_reference.g)
////    regions.clear();
////    regions.reserve(100'000'000); // ?
////    
////    for (const auto& region : get_all_contig_regions(reference)) {
////        auto more_regions = decompose(region, dis(g));
////        regions.insert(regions.end(), more_regions.begin(), more_regions.end());
////    }
////    
////    std::shuffle(regions.begin(), regions.end(), g);
////    
////    uncached_sequence.clear();
////    cached_sequence.clear();
////    
////    uncached_sequence.reserve(get_genome_size(reference));
////    cached_sequence.reserve(get_genome_size(reference));
////    
////    for (const auto& region : regions) {
////        uncached_sequence += uncached_reference.fetch_sequence(region);
////        cached_sequence += cached_reference.fetch_sequence(region);
////    }
////    
////    BOOST_CHECK(uncached_sequence == cached_sequence);
//}

BOOST_AUTO_TEST_CASE(CachingFasta_works_the_same_as_Fasta)
{
    BOOST_REQUIRE(test_file_exists(ecoli_reference_fasta));
    
    constexpr std::size_t max_cache_size {100000};
    
    const auto ecoli = make_reference(ecoli_reference_fasta, max_cache_size);
    
    BOOST_CHECK(ecoli.name() == "R00000042");
    BOOST_CHECK(ecoli.contains_region(GenomicRegion("R00000042", 10000, 2000000)));
    BOOST_CHECK(ecoli.contig_size("R00000042") == 5231428);
    BOOST_CHECK(!ecoli.has_contig("X"));
    BOOST_CHECK(ecoli.contig_region("R00000042") == GenomicRegion("R00000042", 0, 5231428));
    BOOST_CHECK(ecoli.fetch_sequence(GenomicRegion("R00000042", 0, 10)) == "AGCTTTTCAT"); // first line
    BOOST_CHECK(ecoli.fetch_sequence(GenomicRegion("R00000042", 69, 80)) == "CTTCTGAACTG"); // accross lines
    
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta, max_cache_size);
    
    BOOST_CHECK(human.name() == "human_g1k_v37");
    BOOST_CHECK(human.contains_region(GenomicRegion("1", 100, 10000)));
    BOOST_CHECK(!human.contains_region(GenomicRegion("1", 100, 3e8))); // too big
    BOOST_CHECK(human.contig_size("20") == 63025520);
    BOOST_CHECK(human.has_contig("X"));
    BOOST_CHECK(!human.has_contig("y"));
    BOOST_CHECK(human.contig_region("X") == GenomicRegion("X", 0, 155270560));
    BOOST_CHECK(human.fetch_sequence(GenomicRegion("15", 51265690, 51265700)) == "ACAATGTTGT");
    BOOST_CHECK(human.fetch_sequence(GenomicRegion("5", 100000, 100010)) == "AGGAAGTTTC");
}

BOOST_AUTO_TEST_CASE(ReferenceGenome_can_be_made_threadsafe)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta, 0, true);
    
    auto fut1 = std::async(std::launch::async, [&] () {
        return human.fetch_sequence(parse_region("1:1,000,000-1,000,005", human));
    });
    auto fut2 = std::async(std::launch::async, [&] () {
        return human.fetch_sequence(parse_region("2:1,000,000-1,000,005", human));
    });
    auto fut3 = std::async(std::launch::async, [&] () {
        return human.fetch_sequence(parse_region("3:1,000,000-1,000,005", human));
    });
    auto fut4 = std::async(std::launch::async, [&] () {
        return human.fetch_sequence(parse_region("4:1,000,000-1,000,005", human));
    });
    auto fut5 = std::async(std::launch::async, [&] () {
        return human.fetch_sequence(parse_region("5:1,000,000-1,000,005", human));
    });
    
    bool throwed {false};
    
    try {
        BOOST_CHECK(fut1.get() == "GGGCA");
        BOOST_CHECK(fut2.get() == "AAGAA");
        BOOST_CHECK(fut3.get() == "CCAAC");
        BOOST_CHECK(fut4.get() == "CTCCC");
        BOOST_CHECK(fut5.get() == "GTTTT");
    } catch (...) {
        throwed = true;
    }
    
    BOOST_CHECK(!throwed);
}

BOOST_AUTO_TEST_SUITE_END() // IO
BOOST_AUTO_TEST_SUITE_END() // Components

} // namespace test
} // namespace octopus
