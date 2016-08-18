// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <vector>
#include <iterator>
#include <algorithm>
#include <future>

#include <io/reference/reference_genome.hpp>
#include <io/reference/fasta.hpp>
#include <io/reference/caching_fasta.hpp>
#include <utils/mappable_algorithms.hpp>

#include "mock/mock_reference.hpp"

#include "test_common.hpp"

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(io)
BOOST_AUTO_TEST_SUITE(reference)

namespace {
    auto to_ints(const std::vector<std::string>& strings)
    {
        std::vector<int> result(strings.size());
        std::transform(std::cbegin(strings), std::cend(strings), std::begin(result),
                       [] (const auto& str) { return std::stoi(str); });
        return result;
    }
    
    template <typename Container>
    bool is_sorted(const Container& container)
    {
        return std::is_sorted(std::cbegin(container), std::cend(container));
    }
}

BOOST_AUTO_TEST_CASE(reference_genomes_can_be_fasta_files)
{
    BOOST_REQUIRE_NO_THROW(make_reference("/test/data/reference.fa"));
    
    const auto reference = make_reference("/test/data/reference.fa");
    
    auto contigs = reference.contig_names();
    
    BOOST_CHECK_EQUAL(contigs.size(), reference.num_contigs());
    
    BOOST_CHECK(is_sorted(to_ints(contigs)));
    
    BOOST_CHECK(std::all_of(std::cbegin(contigs), std::cend(contigs),
                            [&] (const auto& contig) {
                                return reference.has_contig(contig);
                            }));
}

BOOST_AUTO_TEST_CASE(contigs_are_reported_in_apperance_order)
{
    const auto reference = mock::make_reference();
    
    const auto contigs = reference.contig_names();
    
    BOOST_CHECK(is_sorted(to_ints(contigs)));
}

BOOST_AUTO_TEST_CASE(ReferenceGenome_handles_basic_queries)
{
    BOOST_REQUIRE(test_file_exists(ecoli_reference_fasta));
    
    const auto ecoli = make_reference(ecoli_reference_fasta);
    
    BOOST_CHECK(ecoli.name() == "R00000042");
    BOOST_CHECK(ecoli.contains(GenomicRegion("R00000042", 10000, 2000000)));
    BOOST_CHECK(ecoli.contig_size("R00000042") == 5231428);
    BOOST_CHECK(!ecoli.has_contig("X"));
    BOOST_CHECK(ecoli.contig_region("R00000042") == GenomicRegion("R00000042", 0, 5231428));
    BOOST_CHECK(ecoli.fetch_sequence(GenomicRegion("R00000042", 0, 10)) == "AGCTTTTCAT"); // first line
    BOOST_CHECK(ecoli.fetch_sequence(GenomicRegion("R00000042", 69, 80)) == "CTTCTGAACTG"); // accross lines
    
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    BOOST_CHECK(human.name() == "human_g1k_v37");
    BOOST_CHECK(human.contains(GenomicRegion("1", 100, 10000)));
    BOOST_CHECK(!human.contains(GenomicRegion("1", 100, 3e8))); // too big
    BOOST_CHECK(human.contig_size("20") == 63025520);
    BOOST_CHECK(human.has_contig("X"));
    BOOST_CHECK(!human.has_contig("y"));
    BOOST_CHECK(human.contig_region("X") == GenomicRegion("X", 0, 155270560));
    BOOST_CHECK(human.fetch_sequence(GenomicRegion("15", 51265690, 51265700)) == "ACAATGTTGT");
    BOOST_CHECK(human.fetch_sequence(GenomicRegion("5", 100000, 100010)) == "AGGAAGTTTC");
}

BOOST_AUTO_TEST_CASE(ReferenceGenome_handles_edge_cases)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    BOOST_CHECK(human.fetch_sequence(GenomicRegion {"1", 100, 100}) == "");
}

BOOST_AUTO_TEST_CASE(CachingFasta_works_the_same_as_Fasta)
{
    BOOST_REQUIRE(test_file_exists(ecoli_reference_fasta));
    
    constexpr std::size_t max_cache_size {100000};
    
    const auto ecoli = make_reference(ecoli_reference_fasta, max_cache_size);
    
    BOOST_CHECK(ecoli.name() == "R00000042");
    BOOST_CHECK(ecoli.contains(GenomicRegion("R00000042", 10000, 2000000)));
    BOOST_CHECK(ecoli.contig_size("R00000042") == 5231428);
    BOOST_CHECK(!ecoli.has_contig("X"));
    BOOST_CHECK(ecoli.contig_region("R00000042") == GenomicRegion("R00000042", 0, 5231428));
    BOOST_CHECK(ecoli.fetch_sequence(GenomicRegion("R00000042", 0, 10)) == "AGCTTTTCAT"); // first line
    BOOST_CHECK(ecoli.fetch_sequence(GenomicRegion("R00000042", 69, 80)) == "CTTCTGAACTG"); // accross lines
    
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta, max_cache_size);
    
    BOOST_CHECK(human.name() == "human_g1k_v37");
    BOOST_CHECK(human.contains(GenomicRegion("1", 100, 10000)));
    BOOST_CHECK(!human.contains(GenomicRegion("1", 100, 3e8))); // too big
    BOOST_CHECK(human.contig_size("20") == 63025520);
    BOOST_CHECK(human.has_contig("X"));
    BOOST_CHECK(!human.has_contig("y"));
    BOOST_CHECK(human.contig_region("X") == GenomicRegion("X", 0, 155270560));
    BOOST_CHECK(human.fetch_sequence(GenomicRegion("15", 51265690, 51265700)) == "ACAATGTTGT");
    BOOST_CHECK(human.fetch_sequence(GenomicRegion("5", 100000, 100010)) == "AGGAAGTTTC");
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
