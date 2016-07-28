//
//  candidate_variant_generation_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <algorithm>

#include "test_common.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "variant.hpp"
#include "candidate_generators.hpp"
#include "vcf_reader.hpp"

using std::cout;
using std::endl;

using octopus::CandidateVariantGenerator;
using octopus::CandidateGeneratorBuilder;

BOOST_AUTO_TEST_SUITE(Components)
BOOST_AUTO_TEST_SUITE(CandidateGenerators)

BOOST_AUTO_TEST_CASE(CandidateGeneratorBuilder_can_construct_candidate_generators)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    CandidateGeneratorBuilder builder {};
    builder.set_reference(human);
    builder.set_min_base_quality(0);
    builder.add_generator(CandidateGeneratorBuilder::Generator::Alignment);
    //builder.add_generator(CandidateGeneratorBuilder::Generator::Assembler);
    
    auto generator = builder.build();
    
    BOOST_CHECK(generator.requires_reads());
    
    // TODO: more checks
}

BOOST_AUTO_TEST_CASE(generate_candidates_returns_sorted_and_unique_candidates)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    ReadManager read_manager {NA12878_low_coverage};
    
    const auto sample = read_manager.samples().front();
    
    const auto region = parse_region("1:0-2000000", human);
    
    CandidateGeneratorBuilder builder;
    builder.set_reference(human);
    builder.set_min_base_quality(0);
    builder.add_generator(CandidateGeneratorBuilder::Generator::Alignment);
    
    auto candidate_generator = builder.build();
    
    const auto reads = read_manager.fetch_reads(sample, region);
    
    add_reads(reads, candidate_generator);
    
    auto candidates = candidate_generator.generate_candidates(region);
    
    BOOST_REQUIRE(std::is_sorted(std::cbegin(candidates), std::cend(candidates)));
    BOOST_CHECK(std::unique(std::begin(candidates), std::end(candidates)) == std::end(candidates));
}

BOOST_AUTO_TEST_CASE(AlignmentCandidateVariantGenerator_ignores_snps_with_low_base_qualities)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    const auto region = parse_region("1:22,298,915-22,299,027", human);
    
    ReadManager read_manager {NA12878_low_coverage};
    
    const auto sample = read_manager.samples().front();
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    octopus::AlignmentCandidateVariantGenerator candidate_generator {
        human,
        octopus::AlignmentCandidateVariantGenerator::Options {14}
    };
    
    add_reads(reads, candidate_generator);
    
    auto candidates = candidate_generator.generate_candidates(region);
    
    BOOST_CHECK(candidates.empty());
}

BOOST_AUTO_TEST_CASE(can_specify_the_minimum_number_of_reads_that_must_support_a_candidate)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    // TODO
}

BOOST_AUTO_TEST_CASE(can_specify_the_maximum_size_of_candidates)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    constexpr unsigned max_variant_size {5};
    
    auto candidate_generator = CandidateGeneratorBuilder().set_reference(human)
    .add_generator(CandidateGeneratorBuilder::Generator::Alignment).set_max_variant_size(max_variant_size).build();
    
    ReadManager read_manager {NA12878_low_coverage};
    
    const auto region = parse_region("16:9290000-9300000", human);
    
    const auto sample = read_manager.samples().front();
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    add_reads(reads, candidate_generator);
    
    auto candidates = candidate_generator.generate_candidates(region);
    
    BOOST_CHECK(std::all_of(std::cbegin(candidates), std::cend(candidates),
                            [=] (const auto& candidate) { return region_size(candidate) <= max_variant_size; }));
}

BOOST_AUTO_TEST_CASE(only_insertions_are_included_when_the_max_variant_size_is_zero)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    auto candidate_generator = CandidateGeneratorBuilder().set_reference(human)
    .add_generator(CandidateGeneratorBuilder::Generator::Alignment).set_max_variant_size(0).build();
    
    ReadManager read_manager {NA12878_low_coverage};
    
    const auto region = parse_region("16:9299940-9300055", human);
    
    const auto sample = read_manager.samples().front();
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    add_reads(reads, candidate_generator);
    
    auto candidates = candidate_generator.generate_candidates(region);
    
    BOOST_CHECK(std::all_of(std::cbegin(candidates), std::cend(candidates),
                            [] (const auto& candidate) { return is_insertion(candidate); }));
}

BOOST_AUTO_TEST_CASE(AlignmentCandidateVariantGenerator_includes_all_alleles_in_the_same_region)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    const auto region = parse_region("7:122579662-122579817", human);
    
    ReadManager read_manager {NA12878_low_coverage};
    
    const auto sample = read_manager.samples().front();
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    octopus::AlignmentCandidateVariantGenerator candidate_generator {
        human,
        octopus::AlignmentCandidateVariantGenerator::Options {10}
    };
    
    add_reads(reads, candidate_generator);
    
    auto candidates = candidate_generator.generate_candidates(region);
}

BOOST_AUTO_TEST_CASE(OnlineCandidateVariantGenerator_can_fetch_variants_from_online_web_service)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    const auto region = parse_region("X:10000-10500", human);
    
    octopus::OnlineCandidateVariantGenerator candidate_generator {human};
    
    auto candidates = candidate_generator.generate_candidates(region);
    
    // TODO
}

BOOST_AUTO_TEST_CASE(CandidateVariantGenerator_combines_multiple_generators)
{
    // TODO
}

BOOST_AUTO_TEST_CASE(ExternalCandidateVariantGenerator_gets_candidates_from_vcf)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    octopus::ExternalCandidateVariantGenerator generator {sample_vcf};
    
    const auto region = parse_region("X:10,095,000-10,100,000", human);
    
    auto candidates = generator.generate_candidates(region);
    
    BOOST_CHECK(candidates.size() == 16);
}

BOOST_AUTO_TEST_SUITE_END() // CandidateGenerators
BOOST_AUTO_TEST_SUITE_END() // Components
