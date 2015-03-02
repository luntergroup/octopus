//
//  variant_candidate_generation_test.cpp
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
#include "variant_candidate_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "assembler_candidate_variant_generator.h"
#include "external_variant_candidates.h"
#include "variant_file_factory.h"
#include "variant_file.h"

TEST_CASE("alignment_variant_candidate_generator_test", "[variant_candidate]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam});
    
    VariantCandidateGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human));
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto a_region = parse_region("10:1000000-1000100", human);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    for (const auto& read : reads) {
        candidate_generator.add_read(read);
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    std::cout << "Num alignment candidates: " << candidates.size() << std::endl;
    
//    for (const auto& candidate : candidates) {
//        std::cout << candidate.get_removed_region() << " " << candidate.get_sequence_removed() << " "
//        << candidate.get_sequence_added() << std::endl;
//    }
}

TEST_CASE("assembler_variant_candidate_generator_test", "[variant_candidate]")
{
    const constexpr unsigned kmer_size {15};
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam});
    
    VariantCandidateGenerator candidate_generator {};
    
    candidate_generator.register_generator(
            std::make_unique<AssemblerCandidateVariantGenerator>(human, kmer_size)
    );
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto a_region = parse_region("10:1000000-1000100", human);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    for (const auto& read : reads) {
        candidate_generator.add_read(read);
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    std::cout << "Num assembler candidates: " << candidates.size() << std::endl;
    
    //    for (const auto& candidate : candidates) {
    //        std::cout << candidate.get_removed_region() << " " << candidate.get_sequence_removed() << " "
    //        << candidate.get_sequence_added() << std::endl;
    //    }
}

//TEST_CASE("external_variant_candidate_generator_test", "[variant_candidate]")
//{
//    VariantFileFactory a_variant_file_factory {};
//    VariantFile a_variant_file {a_variant_file_factory.make_impl("some.vcf")};
//    VariantCandidateGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<ExternalVariantCandidates>(a_variant_file));
//    
//    GenomicRegion a_region {"10", 1000000, 1000100};
//    auto candidates = candidate_generator.get_candidates(a_region);
//    std::cout << "Num external candidates: " << candidates.size() << std::endl;
//}

TEST_CASE("all_variant_candidate_generator_test", "[variant_candidate]")
{
    const constexpr unsigned kmer_size {15};
    
    //VariantFileFactory a_variant_file_factory {};
    //VariantFile a_variant_file {a_variant_file_factory.make_impl("some.vcf")};
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam});
    
    VariantCandidateGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human));
    candidate_generator.register_generator(
        std::make_unique<AssemblerCandidateVariantGenerator>(human, kmer_size));
    //candidate_generator.register_generator(std::make_unique<ExternalVariantCandidates>(a_variant_file));
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto a_region = parse_region("10:1000000-1000100", human);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    for (const auto& read : reads) {
        candidate_generator.add_read(read);
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    std::cout << "Total candidates: " << candidates.size() << std::endl;
}
