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
#include <memory>

#include "test_common.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "variant.hpp"
#include "candidate_generators.hpp"
#include "vcf_reader.hpp"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(CandidateVariantGenerator_does_not_give_duplicate_candidates)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    ReadManager read_manager {NA12878_low_coverage};
    
    const auto sample = read_manager.get_samples().front();
    
    const auto region = *parse_region("7:122579662-122579817", human);
    
    Octopus::CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<Octopus::AlignmentCandidateVariantGenerator>(human, 0));
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    add_reads(reads, candidate_generator);
    
    auto candidates = candidate_generator.get_candidates(region);
    
    BOOST_CHECK(candidates.size() == 15);
}

BOOST_AUTO_TEST_CASE(AlignmentCandidateVariantGenerator_ignores_snps_with_low_base_qualities)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    const auto region = *parse_region("7:122579662-122579817", human);
    
    ReadManager read_manager {NA12878_low_coverage};
    
    const auto sample = read_manager.get_samples().front();
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    Octopus::AlignmentCandidateVariantGenerator candidate_generator {human, 10};
    
    add_reads(reads, candidate_generator);
    
    auto candidates = candidate_generator.get_candidates(region);
    
    BOOST_CHECK(candidates.size() == 12);
}

BOOST_AUTO_TEST_CASE(AlignmentCandidateVariantGenerator_includes_all_alleles_in_the_same_region)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    const auto human = make_reference(human_reference_fasta);
    
    const auto region = *parse_region("7:122579662-122579817", human);
    
    ReadManager read_manager {NA12878_low_coverage};
    
    const auto sample = read_manager.get_samples().front();
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    Octopus::AlignmentCandidateVariantGenerator candidate_generator {human, 10};
    
    add_reads(reads, candidate_generator);
    
    auto candidates = candidate_generator.get_candidates(region);
}

BOOST_AUTO_TEST_CASE(OnlineCandidateVariantGenerator_can_fetch_variants_from_online_web_service)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    const auto human = make_reference(human_reference_fasta);
    
    const auto region = *parse_region("X:10000-10500", human);
    
    Octopus::OnlineCandidateVariantGenerator candidate_generator {human};
    
    auto candidates = candidate_generator.get_candidates(region);
    
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
    
    Octopus::ExternalCandidateVariantGenerator generator {sample_vcf};
    
    const auto region = *parse_region("X:10,095,000-10,100,000", human);
    
    auto candidates = generator.get_candidates(region);
    
    BOOST_CHECK(candidates.size() == 16);
}

BOOST_AUTO_TEST_CASE(CandidateGeneratorBuilder_can_construct_candidate_generators)
{
    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
    
    using Octopus::CandidateGeneratorBuilder;
    
    Octopus::CandidateGeneratorBuilder builder {};
    
    builder.add_generator(CandidateGeneratorBuilder::Generator::Alignment);
    
    // TODO
}

BOOST_AUTO_TEST_SUITE_END()
