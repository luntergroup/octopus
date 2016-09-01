// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include "basics/contig_region.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "utils/mappable_algorithms.hpp"
#include "mock_objects.hpp"

using std::cout;
using std::endl;

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(Components)
BOOST_AUTO_TEST_SUITE(MappableAlgorithms)

BOOST_AUTO_TEST_CASE(find_first_after_returns_the_first_element_that_is_after_the_given_region)
{
    const std::vector<ContigRegion> regions {
        ContigRegion {0, 1}, ContigRegion {0, 2}, ContigRegion {2, 5}, ContigRegion {3, 3},
        ContigRegion {3, 4}, ContigRegion {5, 5}, ContigRegion {5, 6}, ContigRegion {5, 8}
    };
    
    BOOST_CHECK_EQUAL(*find_first_after(regions, regions[0]), regions[2]);
    BOOST_CHECK_EQUAL(*find_first_after(regions, regions[1]), regions[2]);
    BOOST_CHECK_EQUAL(*find_first_after(regions, regions[3]), regions[4]);
    BOOST_CHECK_EQUAL(*find_first_after(regions, regions[4]), regions[5]);
    BOOST_CHECK_EQUAL(*find_first_after(regions, regions[5]), regions[6]);
    BOOST_CHECK(find_first_after(regions, regions[6]) == std::cend(regions));
    
    BOOST_CHECK_EQUAL(*find_first_after(regions, ContigRegion {2, 3}), regions[4]);
    BOOST_CHECK_EQUAL(*find_first_after(regions, ContigRegion {2, 4}), regions[5]);
}

BOOST_AUTO_TEST_CASE(overlap_range_returns_an_iterator_range_that_includes_all_overlapped_elements)
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    GenomicRegion test_region {"test", 5000, 6000};
    
    std::vector<GenomicRegion> true_overlaps {};
    std::copy_if(regions.cbegin(), regions.cend(), std::back_inserter(true_overlaps),
                 [&] (const auto& region) { return overlaps(region, test_region); });
    
    const auto overlapped = overlap_range(regions, test_region);
    
    BOOST_CHECK(std::equal(true_overlaps.cbegin(), true_overlaps.cend(), overlapped.begin()));
}

BOOST_AUTO_TEST_CASE(bidirectionally_sorted_ranges_returns_a_set_of_ranges_that_are_themselves_bidirectionally_sorted_and_accounts_for_every_element)
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    auto bisorted_ranges = extract_bidirectionally_sorted_ranges(regions);
    
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
    
    auto bisorted_ranges = extract_bidirectionally_sorted_ranges(regions);
    
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
    
    auto bisorted_ranges = extract_bidirectionally_sorted_ranges(regions);
    
    auto max_region_size = region_size(*largest_mappable(regions));
    
    auto overlapped = overlap_range(regions, test_region, max_region_size);
    
    BOOST_CHECK(std::equal(true_overlaps.cbegin(), true_overlaps.cend(), overlapped.begin()));
}

BOOST_AUTO_TEST_CASE(contained_range_returns_a_range_of_iterators_that_span_all_elements_contained_by_a_region)
{
    auto regions = generate_random_regions(10000, 100, 10000);
    
    GenomicRegion test_region {"test", 5000, 6000};
    
    std::vector<GenomicRegion> true_contained {};
    std::copy_if(regions.cbegin(), regions.cend(), std::back_inserter(true_contained),
                 [&test_region] (const auto& region) { return contains(test_region, region); });
    
    auto contained = contained_range(regions, test_region);
    
    BOOST_CHECK(std::equal(true_contained.cbegin(), true_contained.cend(), contained.begin()));
}

BOOST_AUTO_TEST_CASE(extract_covered_regions_returns_regions_which_contain_all_elements_in_the_given_range_with_no_inter_range_overlaps)
{
    std::vector<GenomicRegion> regions {
        GenomicRegion {"test", 10000, 20000}, GenomicRegion {"test", 20000, 30000},
        GenomicRegion {"test", 30000, 40000}, GenomicRegion {"test", 40001, 50000},
        GenomicRegion {"test", 45000, 60000}
    };
    
    const auto covered = extract_covered_regions(regions);
    
    BOOST_CHECK(covered.size() == 2);
    BOOST_CHECK(covered[0] == GenomicRegion("test", 10000, 40000));
    BOOST_CHECK(covered[1] == GenomicRegion("test", 40001, 60000));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
