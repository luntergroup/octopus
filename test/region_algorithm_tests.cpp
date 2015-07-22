//
//  region_algorithm_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include "mock_objects.h"
#include "genomic_region.h"
#include "mappable_algorithms.h"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(overlap_range_returns_a_filter_iterator_range_that_includes_all_overlapped_elements)
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    GenomicRegion test_region {"test", 5000, 6000};
    
    std::vector<GenomicRegion> true_overlaps {};
    std::copy_if(regions.cbegin(), regions.cend(), std::back_inserter(true_overlaps),
                 [&test_region] (const auto& region) { return overlaps(region, test_region); });
    
    auto overlapped = overlap_range(regions.cbegin(), regions.cend(), test_region);
    
    BOOST_CHECK(std::equal(true_overlaps.cbegin(), true_overlaps.cend(), overlapped.begin()));
}

BOOST_AUTO_TEST_CASE(bidirectionally_sorted_ranges_returns_a_set_of_ranges_that_are_themselves_bidirectionally_sorted_and_accounts_for_every_element)
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    auto bisorted_ranges = bidirectionally_sorted_ranges(regions.cbegin(), regions.cend());
    
    bool are_all_bisorted = std::all_of(bisorted_ranges.cbegin(), bisorted_ranges.cend(),
                                        [] (const auto& range) {
                                            return is_bidirectionally_sorted(range.begin(), range.end());
                                        });
    
    auto s = std::accumulate(bisorted_ranges.cbegin(), bisorted_ranges.cend(), 0,
                             [] (const auto& lhs, const auto& rhs) {
                                 return lhs + size(rhs);
                             });
    
    BOOST_CHECK(s == regions.size());
    BOOST_CHECK(are_all_bisorted);
}

BOOST_AUTO_TEST_CASE(overlap_range_returns_correct_range_if_regions_are_bidirectionally_sorted)
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    GenomicRegion test_region {"test", 5000, 6000};
    
    std::vector<GenomicRegion> true_overlaps {};
    std::copy_if(regions.cbegin(), regions.cend(), std::back_inserter(true_overlaps),
                 [&test_region] (const auto& region) { return overlaps(region, test_region); });
    
    auto bisorted_ranges = bidirectionally_sorted_ranges(regions.cbegin(), regions.cend());
    
    std::vector<GenomicRegion> overlapped {};
    for (const auto& bisorted_range : bisorted_ranges) {
        auto sorted_overlapped = overlap_range(bisorted_range.begin(), bisorted_range.end(), test_region);
        overlapped.insert(overlapped.end(), sorted_overlapped.begin(), sorted_overlapped.end());
    }
    
    BOOST_CHECK(std::equal(true_overlaps.cbegin(), true_overlaps.cend(), overlapped.begin()));
}

BOOST_AUTO_TEST_CASE(overlap_range_returns_correct_range_if_given_the_maximum_region_size)
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    GenomicRegion test_region {"test", 5000, 6000};
    
    std::vector<GenomicRegion> true_overlaps {};
    std::copy_if(regions.cbegin(), regions.cend(), std::back_inserter(true_overlaps),
                 [&test_region] (const auto& region) { return overlaps(region, test_region); });
    
    auto bisorted_ranges = bidirectionally_sorted_ranges(regions.cbegin(), regions.cend());
    
    auto max_region_size = size(*largest_element(regions.cbegin(), regions.cend()));
    
    auto overlapped = overlap_range(regions.cbegin(), regions.cend(), test_region, max_region_size);
    
    BOOST_CHECK(std::equal(true_overlaps.cbegin(), true_overlaps.cend(), overlapped.begin()));
}

BOOST_AUTO_TEST_CASE(contained_range_returns_a_range_of_iterators_that_span_all_elements_contained_by_a_region)
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    GenomicRegion test_region {"test", 5000, 6000};
    
    std::vector<GenomicRegion> true_contained {};
    std::copy_if(regions.cbegin(), regions.cend(), std::back_inserter(true_contained),
                 [&test_region] (const auto& region) { return contains(test_region, region); });
    
    auto contained = contained_range(regions.cbegin(), regions.cend(), test_region);
    
    BOOST_CHECK(std::equal(true_contained.cbegin(), true_contained.cend(), contained.begin()));
}

BOOST_AUTO_TEST_CASE(get_covered_regions_returns_regions_which_contain_all_elements_in_the_given_range_with_no_inter_range_overlaps)
{
    std::vector<GenomicRegion> regions {GenomicRegion {"18", 10000, 20000}, GenomicRegion {"18", 20000, 30000}, GenomicRegion {"18", 30000, 40000}, GenomicRegion {"18", 40001, 50000}, GenomicRegion {"18", 45000, 60000}};
    
    auto covered = get_covered_regions(regions.cbegin(), regions.cend());
    
    BOOST_CHECK(covered.size() == 2);
    BOOST_CHECK(covered[0] == GenomicRegion("18", 10000, 40000));
    BOOST_CHECK(covered[1] == GenomicRegion("18", 40001, 60000));
}

BOOST_AUTO_TEST_CASE(find_first_after_returns_the_first_element_that_is_after_the_given_region)
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    GenomicRegion test1 {"test", 500, 600};
    GenomicRegion test2 {"test", 2000, 2500};
    GenomicRegion test3 {"test", 7000, 9800};
    
    auto after1 = *find_first_after(regions.cbegin(), regions.cend(), test1);
    auto after2 = *find_first_after(regions.cbegin(), regions.cend(), test2);
    auto after3 = *find_first_after(regions.cbegin(), regions.cend(), test3);
    
    auto true_after1 = *std::find_if(regions.cbegin(), regions.cend(),
                                    [&test1] (const auto& region) {
                                        return is_after(region, test1);
                                    });
    
    auto true_after2 = *std::find_if(regions.cbegin(), regions.cend(),
                                    [&test2] (const auto& region) {
                                        return is_after(region, test2);
                                    });
    
    auto true_after3 = *std::find_if(regions.cbegin(), regions.cend(),
                                    [&test3] (const auto& region) {
                                        return is_after(region, test3);
                                    });
    
    BOOST_CHECK(after1 == true_after1);
    BOOST_CHECK(after2 == true_after2);
    BOOST_CHECK(after3 == true_after3);
}

BOOST_AUTO_TEST_SUITE_END()

//BOOST_AUTO_TEST_CASE(", "[region_algorithms]")
//{
//    
//}
//
//BOOST_AUTO_TEST_CASE(", "[region_algorithms]")
//{
//    
//}