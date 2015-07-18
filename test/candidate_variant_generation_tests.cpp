//
//  candidate_variant_generation_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <memory>

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "read_manager.h"
#include "variant.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "assembler_candidate_variant_generator.h"
#include "online_candidate_variant_generator.h"
#include "external_variant_candidates.h"
#include "variant_file_factory.h"
#include "variant_file_reader.h"

using std::cout;
using std::endl;

TEST_CASE("CandidateVariantGenerator does not give duplicate candidates", "[candidates]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    CandidateVariantGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto a_region = parse_region("7:122579662-122579817", human);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    REQUIRE(candidates.size() == 15);
}

TEST_CASE("AlignmentCandidateVariantGenerator ignores snps with low base qualities", "[candidates]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    AlignmentCandidateVariantGenerator candidate_generator {human, 10};
    
    auto a_region = parse_region("7:122579662-122579817", human);
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    REQUIRE(candidates.size() == 12);
}

TEST_CASE("AlignmentCandidateVariantGenerator includes all alleles in the same region", "[candidates]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    AlignmentCandidateVariantGenerator candidate_generator {human, 10};
    
    auto a_region = parse_region("7:122579662-122579817", human);
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
}

TEST_CASE("OnlineCandidateVariantGenerator can fetch variants from online web service", "[candidates]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    OnlineCandidateVariantGenerator candidate_generator {human};
    
    GenomicRegion a_region {"X", 10000, 10500};
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    for (const auto& candidate : candidates) {
        cout << candidate << endl;
    }
}

//TEST_CASE("CandidateVariantGenerator combines multiple generators", "[candidates]")
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
//    auto a_region = parse_region("10:1000000-1000100", human);
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
