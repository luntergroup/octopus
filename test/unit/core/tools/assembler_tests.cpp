// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <exception>

#include <core/tools/vargen/utils/assembler.hpp>

namespace octopus { namespace test {

using octopus::coretools::Assembler;

BOOST_AUTO_TEST_SUITE(core)
BOOST_AUTO_TEST_SUITE(assembler)

BOOST_AUTO_TEST_CASE(assembler_can_be_constructed_with_reference_sequence)
{
    const Assembler::NucleotideSequence reference {"AAAAACCCCC"};
    
    constexpr unsigned kmer_size {5};
    
    Assembler assembler {kmer_size, reference};
    
    BOOST_CHECK(!assembler.is_empty());
    BOOST_CHECK(assembler.is_all_reference());
}

BOOST_AUTO_TEST_CASE(reference_sequence_can_be_inserted_into_an_empty_assembler)
{
    const Assembler::NucleotideSequence reference {"AAAAACCCCC"};
    
    constexpr unsigned kmer_size {5};
    
    Assembler assembler {kmer_size};
    
    BOOST_REQUIRE(assembler.is_empty());
    
    BOOST_REQUIRE_NO_THROW(assembler.insert_reference(reference));
    
    BOOST_CHECK(!assembler.is_empty());
    BOOST_CHECK(assembler.is_all_reference());
}

BOOST_AUTO_TEST_CASE(assemblers_can_be_cleared_and_reused)
{
    const Assembler::NucleotideSequence reference {"AAAAACCCCC"};
    
    constexpr unsigned kmer_size {5};
    
    Assembler assembler {kmer_size, reference};
    
    BOOST_REQUIRE(!assembler.is_empty());
    
    BOOST_REQUIRE_NO_THROW(assembler.clear());
    
    BOOST_CHECK(assembler.is_empty());
    
    BOOST_REQUIRE_NO_THROW(assembler.insert_reference(reference));
    
    BOOST_CHECK(!assembler.is_empty());
}

BOOST_AUTO_TEST_CASE(assembler_throws_if_reference_sequence_is_inserted_twice)
{
    const Assembler::NucleotideSequence reference {"AAAAACCCCC"};
    
    constexpr unsigned kmer_size {5};
    
    Assembler assembler {kmer_size, reference};
    
    BOOST_CHECK_THROW(assembler.insert_reference(reference), std::exception);
}



BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
    
} // namespace test
} // namespace octopus
