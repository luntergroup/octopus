//
//  region_algorithm_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include "mock_objects.h"
#include "genomic_region.h"
#include "region_algorithms.h"

using std::cout;
using std::endl;

TEST_CASE("overlap_range returns a filter iterator range that includes all overlapped elements", "[region_algorithms]")
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    GenomicRegion test_region {"test", 5000, 6000};
    
    std::vector<GenomicRegion> true_overlaps {};
    std::copy_if(regions.cbegin(), regions.cend(), std::back_inserter(true_overlaps),
                 [&test_region] (const auto& region) { return overlaps(region, test_region); });
    
    auto overlapped = overlap_range(regions.cbegin(), regions.cend(), test_region);
    
    REQUIRE(std::equal(true_overlaps.cbegin(), true_overlaps.cend(), overlapped.begin()));
}

TEST_CASE("bidirectionally_sorted_ranges returns a set of ranges that are by themselves bidirectionally sorted and accounts for every element", "[region_algorithms]")
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
    
    REQUIRE(s == regions.size());
    REQUIRE(are_all_bisorted);
}

TEST_CASE("overlap_range returns correct range if regions are bidirectionally sorted", "[region_algorithms]")
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
    
    REQUIRE(std::equal(true_overlaps.cbegin(), true_overlaps.cend(), overlapped.begin()));
}

TEST_CASE("overlap_range returns correct range if given the maximum region size", "[region_algorithms]")
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    GenomicRegion test_region {"test", 5000, 6000};
    
    std::vector<GenomicRegion> true_overlaps {};
    std::copy_if(regions.cbegin(), regions.cend(), std::back_inserter(true_overlaps),
                 [&test_region] (const auto& region) { return overlaps(region, test_region); });
    
    auto bisorted_ranges = bidirectionally_sorted_ranges(regions.cbegin(), regions.cend());
    
    auto max_region_size = size(*largest(regions.cbegin(), regions.cend()));
    
    auto overlapped = overlap_range(regions.cbegin(), regions.cend(), test_region, max_region_size);
    
    REQUIRE(std::equal(true_overlaps.cbegin(), true_overlaps.cend(), overlapped.begin()));
}

TEST_CASE("contained_range returns a range of iterators that span all elements contained by a region", "[region_algorithms]")
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    GenomicRegion test_region {"test", 5000, 6000};
    
    std::vector<GenomicRegion> true_contained {};
    std::copy_if(regions.cbegin(), regions.cend(), std::back_inserter(true_contained),
                 [&test_region] (const auto& region) { return contains(test_region, region); });
    
    auto contained = contained_range(regions.cbegin(), regions.cend(), test_region);
    
    REQUIRE(std::equal(true_contained.cbegin(), true_contained.cend(), contained.begin()));
}

TEST_CASE("minimal_encompassing returns regions which contain all elements in the given range, with no inter-range-overlaps", "[region_algorithms]")
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    auto sub_regions = minimal_encompassing(regions.cbegin(), regions.cend());
    
    std::vector<GenomicRegion> contained_regions {};
    
    for (const auto& region : sub_regions) {
        auto contained = contained_range(regions.cbegin(), regions.cend(), region);
        contained_regions.insert(contained_regions.end(), contained.begin(), contained.end());
    }
    
    REQUIRE(std::equal(regions.cbegin(), regions.cend(), contained_regions.begin()));
}

TEST_CASE("find_first_after returns the first element that is_after the given region", "[region_algorithms]")
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
    
    REQUIRE(after1 == true_after1);
    REQUIRE(after2 == true_after2);
    REQUIRE(after3 == true_after3);
}

TEST_CASE("", "[region_algorithms]")
{
    
}

TEST_CASE("", "[region_algorithms]")
{
    
}
