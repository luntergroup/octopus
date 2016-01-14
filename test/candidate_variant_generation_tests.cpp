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
    auto human = make_reference(human_reference_fasta);
    
    ReadManager a_read_manager {HG00101};
    
    Octopus::CandidateVariantGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<Octopus::AlignmentCandidateVariantGenerator>(human, 0));
    
    auto sample_ids = a_read_manager.get_samples();
    auto the_sample_id = sample_ids.at(0);
    
    auto a_region = *parse_region("7:122579662-122579817", human);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    BOOST_CHECK(candidates.size() == 15);
}

BOOST_AUTO_TEST_CASE(AlignmentCandidateVariantGenerator_ignores_snps_with_low_base_qualities)
{
    auto human = make_reference(human_reference_fasta);
    
    ReadManager a_read_manager {HG00101};
    
    Octopus::AlignmentCandidateVariantGenerator candidate_generator {human, 10};
    
    auto a_region = *parse_region("7:122579662-122579817", human);
    
    auto sample_ids = a_read_manager.get_samples();
    auto the_sample_id = sample_ids.at(0);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    BOOST_CHECK(candidates.size() == 12);
}

BOOST_AUTO_TEST_CASE(AlignmentCandidateVariantGenerator_includes_all_alleles_in_the_same_region)
{
    auto human = make_reference(human_reference_fasta);
    
    ReadManager a_read_manager {HG00101};
    
    Octopus::AlignmentCandidateVariantGenerator candidate_generator {human, 10};
    
    auto a_region = *parse_region("7:122579662-122579817", human);
    
    auto sample_ids = a_read_manager.get_samples();
    auto the_sample_id = sample_ids.at(0);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
}

BOOST_AUTO_TEST_CASE(OnlineCandidateVariantGenerator_can_fetch_variants_from_online_web_service)
{
    auto human = make_reference(human_reference_fasta);
    
    Octopus::OnlineCandidateVariantGenerator candidate_generator {human};
    
    GenomicRegion a_region {"X", 10000, 10500};
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    for (const auto& candidate : candidates) {
        cout << candidate << endl;
    }
}

//BOOST_AUTO_TEST_CASE(CandidateVariantGenerator_combines_multiple_generators)
//{
//    const constexpr unsigned kmer_size {15};
//    
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
//    
//    CandidateVariantGenerator candidate_generator {};
//    
//    candidate_generator.register_generator(
//            std::make_unique<AssemblerCandidateVariantGenerator>(human, kmer_size, 0)
//    );
//    
//    auto sample_ids = a_read_manager.get_sample_ids();
//    auto the_sample_id = sample_ids.at(0);
//    
//    auto a_region = *parse_region("10:1000000-1000100", human);
//    
//    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
//    
//    for (const auto& read : reads) {
//        candidate_generator.add_read(read);
//    }
//    
//    auto candidates = candidate_generator.get_candidates(a_region);
//    
//    //std::cout << "Num assembler candidates: " << candidates.size() << std::endl;
//    
//    //    for (const auto& candidate : candidates) {
//    //        std::cout << candidate.get_removed_region() << " " << candidate.get_sequence_removed() << " "
//    //        << candidate.get_sequence_added() << std::endl;
//    //    }
//}

BOOST_AUTO_TEST_CASE(ExternalCandidateVariantGenerator_gets_candidates_from_vcf)
{
    auto reference = make_reference(human_reference_fasta);
    
    Octopus::ExternalCandidateVariantGenerator generator {sample_vcf};
    
    auto region = *parse_region("X:10,095,000-10,100,000", reference);
    
    auto candidates = generator.get_candidates(region);
    
    BOOST_CHECK(candidates.size() == 16);
}

BOOST_AUTO_TEST_SUITE_END()
