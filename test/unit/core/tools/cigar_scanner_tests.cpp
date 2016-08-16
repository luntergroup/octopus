// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <core/tools/vargen/variant_generator.hpp>

namespace octopus { namespace test {

namespace {
    auto make_cigar_scanner()
    {
        using Builder = VariantGenerator::Builder;
        return Builder().add_generator(Builder::Generator::Alignment);
    }
}

BOOST_AUTO_TEST_SUITE(core)
BOOST_AUTO_TEST_SUITE(cigar_scanner)



BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
    
} // namespace test
} // namespace octopus
