// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <iterator>
#include <vector>

#include <boost/filesystem.hpp>

#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "io/read/read_manager.hpp"
#include "test_common.hpp"
#include "mock_objects.hpp"

namespace octopus { namespace test {

namespace fs = boost::filesystem;

BOOST_AUTO_TEST_SUITE(io)
BOOST_AUTO_TEST_SUITE(read_manager)

BOOST_AUTO_TEST_CASE(can_extract_reads_from_BAM_files)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    ReadManager read_manager {NA12878_low_coverage};
    
    BOOST_REQUIRE(read_manager.num_samples() == 1);
    
    const auto sample = read_manager.samples().front();
    
    const GenomicRegion region1 {"1", 9'990, 10'000};
    const GenomicRegion region2 {"10", 1'000'000, 1'000'100};
    const GenomicRegion region3 {"3", 100'000, 100'100};
    
    const auto reads1 = read_manager.fetch_reads(sample, region1);
    
    BOOST_CHECK(reads1.size() == 1);
    
    const auto reads2 = read_manager.fetch_reads(sample, region2);
    
    BOOST_CHECK(reads2.size() == 7);
    
    const auto reads3 = read_manager.fetch_reads(sample, region3);
    
    BOOST_CHECK(reads3.size() == 21);
}

BOOST_AUTO_TEST_CASE(can_extract_reads_from_CRAM_files)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage_cram));
    
    ReadManager read_manager {NA12878_low_coverage_cram};
    
    BOOST_REQUIRE(read_manager.num_samples() == 1);
    
    const auto sample = read_manager.samples().front();
    
    const GenomicRegion region1 {"1", 9'990, 10'000};
    const GenomicRegion region2 {"10", 1'000'000, 1'000'100};
    const GenomicRegion region3 {"3", 100'000, 100'100};
    
    const auto reads1 = read_manager.fetch_reads(sample, region1);
    
    BOOST_CHECK(reads1.size() == 1);
    
    const auto reads2 = read_manager.fetch_reads(sample, region2);
    
    BOOST_CHECK(reads2.size() == 7);
    
    const auto reads3 = read_manager.fetch_reads(sample, region3);
    
    BOOST_CHECK(reads3.size() == 21);
}

BOOST_AUTO_TEST_CASE(read_manager_multiple_files_below_max_file_limit_test)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    BOOST_REQUIRE(test_file_exists(HG00101));
    
    std::vector<fs::path> read_paths {NA12878_low_coverage, HG00101};
    
    constexpr unsigned maxOpenFiles {2};
    
    ReadManager read_manager(read_paths, maxOpenFiles);
    
    BOOST_REQUIRE(read_manager.num_samples() == 2);
    
    const auto samples = read_manager.samples();
    
    BOOST_REQUIRE(samples.size() == 2);
    
    const GenomicRegion a_big_region {"1", 2000000, 3000000};
    const GenomicRegion a_small_region {"10", 1000000, 1000100};
    
    const auto first_sample = samples.front();
    
    const auto big_region_reads1   = read_manager.fetch_reads(first_sample, a_big_region);
    const auto small_region_reads1 = read_manager.fetch_reads(first_sample, a_small_region);
    
    BOOST_CHECK(big_region_reads1.size() == 61225);
    BOOST_CHECK(small_region_reads1.size() == 10);
    
    const auto second_sample = samples.back();
    
    const auto big_region_reads2   = read_manager.fetch_reads(second_sample, a_big_region);
    const auto small_region_reads2 = read_manager.fetch_reads(second_sample, a_small_region);
    
    BOOST_CHECK(big_region_reads2.size() == 39200);
    BOOST_CHECK(small_region_reads2.size() == 7);
}

BOOST_AUTO_TEST_CASE(read_manager_multiple_files_above_max_file_limit_test)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    BOOST_REQUIRE(test_file_exists(HG00101));
    BOOST_REQUIRE(test_file_exists(HG00102));
    
    std::vector<fs::path> read_paths {NA12878_low_coverage, HG00101, HG00102};
    
    constexpr unsigned maxOpenFiles {2};
    
    ReadManager read_manager(read_paths, maxOpenFiles);
    
    BOOST_REQUIRE(read_manager.num_samples() == 3);
    
    const auto samples = read_manager.samples();
    
    const GenomicRegion a_big_region {"1", 2000000, 3000000};
    const GenomicRegion a_small_region {"10", 1000000, 1000100};
    
    const auto big_reads1 = read_manager.fetch_reads(samples[0], a_big_region);
    const auto big_reads2 = read_manager.fetch_reads(samples[1], a_big_region);
    const auto big_reads3 = read_manager.fetch_reads(samples[2], a_big_region);
    
    BOOST_CHECK(big_reads1.size() == 61225);
    BOOST_CHECK(big_reads2.size() == 49810);
    BOOST_CHECK(big_reads3.size() == 39200);
    
    const auto small_reads1 = read_manager.fetch_reads(samples[0], a_small_region);
    const auto small_reads2 = read_manager.fetch_reads(samples[1], a_small_region);
    const auto small_reads3 = read_manager.fetch_reads(samples[2], a_small_region);
    
    BOOST_CHECK(small_reads1.size() == 10);
    BOOST_CHECK(small_reads2.size() == 9);
    BOOST_CHECK(small_reads3.size() == 7);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
