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

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(ReferemceGenome_handles_basic_queries)
{
    auto ecoli = make_reference(ecoli_reference_fasta);
    
    BOOST_CHECK(ecoli.get_name() == "R00000042");
    BOOST_CHECK(ecoli.contains_region(GenomicRegion("R00000042", 10000, 2000000)));
    BOOST_CHECK(ecoli.get_contig_size("R00000042") == 5231428);
    BOOST_CHECK(!ecoli.has_contig("X"));
    BOOST_CHECK(ecoli.get_contig_region("R00000042") == GenomicRegion("R00000042", 0, 5231428));
    BOOST_CHECK(ecoli.get_sequence(GenomicRegion("R00000042", 0, 10)) == "AGCTTTTCAT"); // first line
    BOOST_CHECK(ecoli.get_sequence(GenomicRegion("R00000042", 69, 80)) == "CTTCTGAACTG"); // accross lines
    
    auto human = make_reference(human_reference_fasta);
    
    BOOST_CHECK(human.get_name() == "human_g1k_v37");
    BOOST_CHECK(human.contains_region(GenomicRegion("1", 100, 10000)));
    BOOST_CHECK(!human.contains_region(GenomicRegion("1", 100, 3e8))); // too big
    BOOST_CHECK(human.get_contig_size("20") == 63025520);
    BOOST_CHECK(human.has_contig("X"));
    BOOST_CHECK(!human.has_contig("y"));
    BOOST_CHECK(human.get_contig_region("X") == GenomicRegion("X", 0, 155270560));
    BOOST_CHECK(human.get_sequence(GenomicRegion("15", 51265690, 51265700)) == "ACAATGTTGT");
    BOOST_CHECK(human.get_sequence(GenomicRegion("5", 100000, 100010)) == "AGGAAGTTTC");
}

BOOST_AUTO_TEST_CASE(ReferemceGenome_handles_edge_cases)
{
    auto human = make_reference(human_reference_fasta);
    
    BOOST_CHECK(human.get_sequence(GenomicRegion {"1", 100, 100}) == "");
}

BOOST_AUTO_TEST_CASE(parse_region_works_with_correctly_formatted_region_input)
{
    auto human = make_reference(human_reference_fasta);
    
    auto r1 = parse_region("3", human);
    BOOST_CHECK(r1.get_contig_name() == "3");
    BOOST_CHECK(r1.get_begin() == 0);
    BOOST_CHECK(r1.get_end() == human.get_contig_size("3"));
    
    auto r2 = parse_region("10:100-200", human);
    BOOST_CHECK(r2.get_contig_name() == "10");
    BOOST_CHECK(r2.get_begin() == 100);
    BOOST_CHECK(r2.get_end() == 200);
    
    auto r3 = parse_region("18:102029", human);
    BOOST_CHECK(r3.get_contig_name() == "18");
    BOOST_CHECK(r3.get_begin() == 102029);
    BOOST_CHECK(r3.get_end() == 102029);
    
    auto r4 = parse_region("MT:100-", human);
    BOOST_CHECK(r4.get_contig_name() == "MT");
    BOOST_CHECK(r4.get_begin() == 100);
    BOOST_CHECK(r4.get_end() == human.get_contig_size("MT"));
}

BOOST_AUTO_TEST_CASE(parse_region_throws_when_region_is_not_formatted_correctly)
{
    auto human = make_reference(human_reference_fasta);
    
    bool throwed {};
    
    try {
        auto r1 = parse_region("-", human);
        throwed = false;
    } catch (...) {
        throwed = true;
    }
    
    BOOST_CHECK(throwed);
    
    try {
        auto r2 = parse_region("5:100-99", human);
        throwed = false;
    } catch (...) {
        throwed = true;
    }
    
    BOOST_CHECK(throwed);
    
//    try {
//        auto r3 = parse_region("2::0-100", human);
//        throwed = false;
//    } catch (...) {
//        throwed = true;
//    }
//    
//    BOOST_CHECK(throwed);
}

BOOST_AUTO_TEST_CASE(cached_and_uncached_reference_genome_give_same_sequence)
{
    // First test simple walk along a single contig
    
    std::string contig {"1"};
    
    std::random_device rd;
    std::mt19937 g(rd());
    
    std::uniform_int_distribution<GenomicRegion::SizeType> dis(1, 50000);
    
    auto uncached_reference = make_reference(human_reference_fasta, 0);
    auto cached_reference   = make_reference(human_reference_fasta, 10'000'000);
    
    std::string uncached_sequence {}, cached_sequence {};
    uncached_sequence.reserve(uncached_reference.get_contig_size(contig));
    cached_sequence.reserve(cached_reference.get_contig_size(contig));
    
    auto regions = decompose(parse_region(contig, uncached_reference), 100);
    
    for (const auto& region : regions) {
        uncached_sequence += uncached_reference.get_sequence(region);
        cached_sequence += cached_reference.get_sequence(region);
    }
    
    BOOST_CHECK(uncached_sequence == cached_sequence);
    
    exit(0);
    
//    // Next test walks along multiple contigs
//    
//    regions.clear();
//    
//    for (const auto& c : std::vector<std::string> {"1", "2", "3"}) {
//        auto more_regions = decompose(uncached_reference.get_contig_region(c), dis(g));
//        regions.insert(regions.end(), more_regions.begin(), more_regions.end());
//    }
//    
//    
    
    // Next test accesses to random points on a single contig
    
    uncached_sequence.clear();
    cached_sequence.clear();
    
    std::shuffle(regions.begin(), regions.end(), g);
    
    for (const auto& region : regions) {
        uncached_sequence += uncached_reference.get_sequence(region);
        cached_sequence += cached_reference.get_sequence(region);
    }
    
    BOOST_CHECK(uncached_sequence == cached_sequence);
    
    // Finally test random access to lots of contigs with different region sizes - may take a while
    
    //uncached_sequence.reserve(uncached_reference.g)
    regions.clear();
    regions.reserve(100'000'000); // ?
    
    for (const auto& region : get_all_contig_regions(uncached_reference)) {
        auto more_regions = decompose(region, dis(g));
        regions.insert(regions.end(), more_regions.begin(), more_regions.end());
    }
    
    std::shuffle(regions.begin(), regions.end(), g);
    
    uncached_sequence.clear();
    cached_sequence.clear();
    
    uncached_sequence.reserve(get_genome_size(uncached_reference));
    cached_sequence.reserve(get_genome_size(cached_reference));
    
    for (const auto& region : regions) {
        uncached_sequence += uncached_reference.get_sequence(region);
        cached_sequence += cached_reference.get_sequence(region);
    }
    
    BOOST_CHECK(uncached_sequence == cached_sequence);
}

//BOOST_AUTO_TEST_CASE(CachingFasta_works_the_same_as_Fasta)
//{
//    ReferenceGenomeFactory a_factory {};
//    
//    ReferenceGenome human_normal {a_factory.make(human_reference_fasta)};
//    ReferenceGenome human_cached {a_factory.make(human_reference_fasta, 200)};
//    
////    GenomicRegion region1 {"1", 1000000, 1000100};
////    GenomicRegion region2 {"1", 1000030, 1000080};
////    GenomicRegion region3 {"1",  999990, 1000110};
////    GenomicRegion region4 {"1",  1000200, 1000300};
////    GenomicRegion region5 {"1",  1000300, 1000400};
////    GenomicRegion region6 {"1",  1000350, 1000450};
//    
//    GenomicRegion region1 {"1", 1000000, 1000100};
//    GenomicRegion region2 {"1", 1000100, 1000200};
//    GenomicRegion region3 {"1", 1000200, 1000300};
//    
//    cout << human_cached.get_sequence(region1) << endl;
//    cout << human_cached.get_sequence(region2) << endl;
//    cout << human_cached.get_sequence(region3) << endl;
//    
////    cout << human_cached.get_sequence(region1) << endl;
////    cout << human_cached.get_sequence(region2) << endl;
////    cout << human_cached.get_sequence(region3) << endl;
////    cout << human_cached.get_sequence(region4) << endl;
////    cout << human_cached.get_sequence(region5) << endl;
////    cout << human_cached.get_sequence(region6) << endl;
////    
////    cout << human_cached.get_sequence(region1) << endl;
////    cout << human_cached.get_sequence(region2) << endl;
////    cout << human_cached.get_sequence(region3) << endl;
////    cout << human_cached.get_sequence(region4) << endl;
////    cout << human_cached.get_sequence(region5) << endl;
////    cout << human_cached.get_sequence(region6) << endl;
//}

//BOOST_AUTO_TEST_CASE(ReferenceGenome_can_be_made_threadsafe)
//{
//    auto reference = make_reference(human_reference_fasta, 0, true);
//    
//    auto fut1 = std::async(std::launch::async, [&reference] () {
//        return reference.get_sequence(parse_region("1:1,000,000-1,000,100", reference));
//    });
//    auto fut2 = std::async(std::launch::async, [&reference] () {
//        return reference.get_sequence(parse_region("2:1,000,000-1,000,100", reference));
//    });
//    auto fut3 = std::async(std::launch::async, [&reference] () {
//        return reference.get_sequence(parse_region("3:1,000,000-1,000,100", reference));
//    });
//    auto fut4 = std::async(std::launch::async, [&reference] () {
//        return reference.get_sequence(parse_region("4:1,000,000-1,000,100", reference));
//    });
//    auto fut5 = std::async(std::launch::async, [&reference] () {
//        return reference.get_sequence(parse_region("5:1,000,000-1,000,100", reference));
//    });
//    
//    cout << fut1.get() << endl;
//    cout << fut2.get() << endl;
//    cout << fut3.get() << endl;
//    cout << fut4.get() << endl;
//    cout << fut5.get() << endl;
//}

BOOST_AUTO_TEST_SUITE_END()
