// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <vector>
#include <iterator>
#include <algorithm>

#include "basics/contig_region.hpp"
#include "containers/mappable_flat_set.hpp"

namespace octopus { namespace test {

using octopus::MappableFlatSet;

BOOST_AUTO_TEST_SUITE(containers)
BOOST_AUTO_TEST_SUITE(mappable_flat_set)

BOOST_AUTO_TEST_CASE(emplace_works)
{
    MappableFlatSet<ContigRegion> set {};
    
    set.emplace(0, 0);
    BOOST_REQUIRE_EQUAL(set.size(), 1);
    set.emplace(0, 0);
    BOOST_CHECK_EQUAL(set.size(), 1);
    set.emplace(0, 1);
    BOOST_REQUIRE_EQUAL(set.size(), 2);
    set.emplace(0, 1);
    BOOST_CHECK_EQUAL(set.size(), 2);
    
    set.emplace(0, 3);
    set.emplace(1, 1);
    set.emplace(2, 4);
    set.emplace(4, 5);
    
    BOOST_CHECK_EQUAL(set.size(), 6);
    BOOST_CHECK(std::is_sorted(std::cbegin(set), std::cend(set)));
}

BOOST_AUTO_TEST_CASE(range_insert_acceptes_empty_ranges)
{
    MappableFlatSet<ContigRegion> set {};
    
    std::vector<ContigRegion> empty {};
    
    BOOST_CHECK_NO_THROW(set.insert(std::cbegin(empty), std::cend(empty)));
    BOOST_CHECK(set.empty());
    
    set.emplace(0, 1);
    
    BOOST_CHECK_NO_THROW(set.insert(std::cbegin(empty), std::cend(empty)));
    BOOST_CHECK_EQUAL(set.size(), 1);
}

BOOST_AUTO_TEST_CASE(range_insert_accepts_unsorted_ranges)
{
    const ContigRegion r1 {0, 1}, r2 {0, 2}, r3 {1, 1}, r4 {0, 0}, r5 {0, 4}, r6 {2, 2}, r7 {0, 3}, r8 {1, 2};
    
    std::vector<ContigRegion> unsorted {r1, r2, r3, r4, r5, r6, r7, r8};
    
    MappableFlatSet<ContigRegion> set {};
    
    set.insert(std::cbegin(unsorted), std::cend(unsorted));
    
    BOOST_CHECK_EQUAL(set.size(), unsorted.size());
    BOOST_CHECK(std::is_sorted(std::cbegin(set), std::cend(set)));
}

BOOST_AUTO_TEST_CASE(insert_hint_works)
{
    const ContigRegion r1 {0, 1}, r2 {0, 2}, r3 {1, 1}, r4 {0, 0}, r5 {0, 4},
                        r6 {2, 2}, r7 {1, 1}, r8 {0, 2}, r9 {0, 3}, r10 {1, 2};
    
    std::vector<ContigRegion> regions {r1, r2, r3, r4, r5, r6, r7, r8, r9, r10};
    
    MappableFlatSet<ContigRegion> set {};
    
    // std::inserter calls the hinted insert version
    std::copy(std::begin(regions), std::end(regions), std::inserter(set, std::begin(set)));
    
    auto unique_regions = regions;
    
    std::sort(std::begin(unique_regions), std::end(unique_regions));
    
    unique_regions.erase(std::unique(std::begin(unique_regions), std::end(unique_regions)),
                         std::end(unique_regions));
    
    BOOST_REQUIRE_EQUAL(set.size(), unique_regions.size());
    
    BOOST_CHECK(std::equal(std::cbegin(set), std::cend(set), std::cbegin(unique_regions)));
    
    set.clear();
    
    BOOST_REQUIRE(set.empty());
    
    std::copy(std::make_move_iterator(std::begin(regions)),
              std::make_move_iterator(std::end(regions)),
              std::inserter(set, std::begin(set)));
    
    BOOST_REQUIRE_EQUAL(set.size(), unique_regions.size());
    
    BOOST_CHECK(std::equal(std::cbegin(set), std::cend(set), std::cbegin(unique_regions)));
    
    auto it = set.insert(std::cbegin(set), r2);
    
    BOOST_CHECK_EQUAL(*it, r2);
    
    it = set.insert(std::cend(set), r2);
    
    BOOST_CHECK_EQUAL(*it, r2);
    
    it = set.insert(std::cbegin(set), r3);
    
    BOOST_CHECK_EQUAL(*it, r3);
    
    it = set.insert(std::cend(set), r3);
    
    BOOST_CHECK_EQUAL(*it, r3);
    
    const ContigRegion r11 {0, 5};
    
    it = set.insert(std::next(std::cbegin(set), 3), r11);
    
    BOOST_CHECK_EQUAL(set.size(), unique_regions.size() + 1);
    BOOST_CHECK_EQUAL(*it, r11);
    BOOST_CHECK(std::is_sorted(std::cbegin(set), std::cend(set)));
}

BOOST_AUTO_TEST_CASE(erase_maybe_removes_a_single_element)
{
    MappableFlatSet<ContigRegion> set {};
    
    ContigRegion r1 {}, r2 {0, 0}, r3 {0, 1}, r4 {0, 2}, r5 {1, 1}, r6 {1, 5};
    
    BOOST_CHECK_NO_THROW(set.erase(r1));
    
    set.emplace(0, 0);
    
    set.erase(r1);
    
    BOOST_CHECK(set.empty());
    
    set.insert(r2);
    set.erase(r2);
    
    BOOST_CHECK(set.empty());
    
    set.insert(r3);
    set.erase(r3);
    
    BOOST_CHECK(set.empty());
    
    set.insert(r2);
    set.insert(r3);
    set.insert(r4);
    set.insert(r5);
    
    set.erase(r4);
    
    BOOST_CHECK_EQUAL(set.size(), 3);
    
    set.insert(r4);
    set.insert(r6);
    
    set.erase(r5);
    BOOST_CHECK_EQUAL(set.size(), 4);
    BOOST_CHECK(std::is_sorted(std::cbegin(set), std::cend(set)));
    set.erase(r3);
    BOOST_CHECK_EQUAL(set.size(), 3);
    BOOST_CHECK(std::is_sorted(std::cbegin(set), std::cend(set)));
    set.erase(r6);
    BOOST_CHECK_EQUAL(set.size(), 2);
    BOOST_CHECK(std::is_sorted(std::cbegin(set), std::cend(set)));
    set.erase(r4);
    BOOST_CHECK_EQUAL(set.size(), 1);
    BOOST_CHECK(std::is_sorted(std::cbegin(set), std::cend(set)));
    set.erase(r2);
    BOOST_CHECK(set.empty());
}
    
BOOST_AUTO_TEST_CASE(erase_all_works_for_empty_range)
{
    MappableFlatSet<ContigRegion> set {};
    
    set.emplace(0, 0);
    set.emplace(0, 1);
    set.emplace(1, 1);
    
    std::vector<ContigRegion> empty {};
    
    BOOST_REQUIRE_NO_THROW(set.erase_all(std::cbegin(empty), std::cend(empty)));
    
    BOOST_CHECK_EQUAL(set.size(), 3);
}
    
BOOST_AUTO_TEST_CASE(erase_all_is_unchanged_if_there_are_no_matches)
{
    const ContigRegion r1 {0, 1}, r2 {0, 2}, r3 {1, 1}, r4 {0, 0}, r5 {0, 4},
    r6 {2, 2}, r7 {1, 1}, r8 {0, 2}, r9 {0, 3}, r10 {1, 2};
    
    std::vector<ContigRegion> regions {r1, r2, r3, r4, r5, r6, r7, r8, r9, r10};
    
    std::sort(std::begin(regions), std::end(regions));
    
    MappableFlatSet<ContigRegion> set {};
    
    auto n = set.erase_all(std::cbegin(regions), std::cend(regions));
    
    BOOST_CHECK_EQUAL(n, 0);
    BOOST_CHECK(set.empty());
    
    set.insert(std::cbegin(regions), std::cend(regions));
    
    std::vector<ContigRegion> others {ContigRegion {1, 3}, ContigRegion {3, 3}, ContigRegion {3, 4}};
    
    const auto s = set.size();
    
    n = set.erase_all(std::cbegin(others), std::cend(others));
    
    BOOST_CHECK_EQUAL(n, 0);
    BOOST_CHECK_EQUAL(set.size(), s);
}

BOOST_AUTO_TEST_CASE(erase_all_only_erases_given_elements_and_maintains_order)
{
    std::vector<ContigRegion> tmp {};
    
    MappableFlatSet<ContigRegion> set {};
    
    set.emplace(0, 0);
    
    tmp.emplace_back(0, 0);
    
    set.erase_all(std::cbegin(tmp), std::cend(tmp));
    
    BOOST_CHECK(set.empty());
    
    set.emplace(0, 0);
    set.emplace(0, 1);
    set.emplace(1, 1);
    set.emplace(0, 3);
    
    set.erase_all(std::cbegin(tmp), std::cend(tmp));
    
    BOOST_CHECK_EQUAL(set.size(), 3);
    
    tmp.emplace_back(0, 3);
    tmp.emplace_back(1, 1);
    tmp.emplace_back(1, 2);
    
    set.erase_all(std::cbegin(tmp), std::cend(tmp));
    
    BOOST_CHECK_EQUAL(set.size(), 1);
    
    // Now a tricky one
    
    set.clear();
    tmp.clear();
    
    BOOST_REQUIRE(set.empty());
    
    set.emplace(81, 82);
    set.emplace(136, 137);
    set.emplace(138, 138);
    set.emplace(163, 171);
    set.emplace(164, 164);
    set.emplace(164, 165);
    set.emplace(165, 166);
    set.emplace(166, 167);
    set.emplace(167, 168);
    set.emplace(168, 169);
    set.emplace(169, 170);
    set.emplace(170, 171);
    set.emplace(170, 179);
    set.emplace(171, 172);
    set.emplace(173, 174);
    
    tmp.emplace_back(163, 171);
    tmp.emplace_back(164, 164);
    tmp.emplace_back(164, 165);
    tmp.emplace_back(170, 171);
    tmp.emplace_back(170, 179);
    
    auto s = set.size();
    
    set.erase_all(std::cbegin(tmp), std::cend(tmp));
    
    BOOST_CHECK_EQUAL(set.size(), s - tmp.size());
    BOOST_CHECK(std::is_sorted(std::cbegin(set), std::cend(set)));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
