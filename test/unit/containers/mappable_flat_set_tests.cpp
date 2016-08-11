// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <vector>
#include <iterator>
#include <algorithm>

#include <basics/contig_region.hpp>
#include <containers/mappable_flat_set.hpp>

namespace octopus { namespace test {

using octopus::MappableFlatSet;

BOOST_AUTO_TEST_SUITE(containers)
BOOST_AUTO_TEST_SUITE(mappable_flat_set)

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

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
