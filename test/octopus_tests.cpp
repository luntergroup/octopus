//
//  octopus_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

/*
    This file contains some general tests for sections of the overall Octopus workflow
 */

#include "catch.hpp"

#include <iostream>
#include <string>
#include <iterator>
#include <vector>
#include <memory>
#include <algorithm>

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "test_common.h"
#include "htslib_facade.h"
#include "read_manager.h"
#include "mock_objects.h"
#include "read_filter.h"
#include "read_filters.h"
#include "read_transform.h"
#include "read_transformations.h"
#include "variant.h"
#include "variant_factory.h"
#include "variant_candidate_generator.h"
#include "variant_utils.h"
#include "alignment_candidate_variant_generator.h"
#include "assembler_candidate_variant_generator.h"
#include "external_variant_candidates.h"
#include "variant_file_factory.h"
#include "variant_file.h"

TEST_CASE("read_filter_transform_generate_left_align_test", "[octopus]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam});
    VariantFactory a_variant_factory {};
    
    std::vector<GenomicRegion> some_regions {};
    some_regions.emplace_back(parse_region("1:10000000-20000000", human));
    some_regions.emplace_back(parse_region("18:389289-399289", human));
    some_regions.emplace_back(parse_region("X:5000000-5010000", human));
    some_regions.emplace_back(parse_region("Y", human));
    
    auto samples = a_read_manager.get_sample_ids();
    
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    ReadFilter<ReadIterator> a_read_filter {};
    a_read_filter.register_filter(is_not_secondary_alignment);
    a_read_filter.register_filter([] (const AlignedRead& the_read) {
        return is_good_mapping_quality(the_read, 20);
    });
    a_read_filter.register_filter([] (const AlignedRead& the_read) {
        return has_sufficient_good_quality_bases(the_read, 20, 10);
    });
    a_read_filter.register_filter(is_not_duplicate<ReadIterator>);
    
    ReadTransform a_read_transform {};
    a_read_transform.register_transform(trim_adapters);
    a_read_transform.register_transform(trim_soft_clipped);
    
    unsigned kmer_size {15};
    VariantCandidateGenerator candidate_generator {};
    candidate_generator.register_generator(
            std::make_unique<AlignmentCandidateVariantGenerator>(human, a_variant_factory));
//    candidate_generator.register_generator(
//            std::make_unique<AssemblerCandidateVariantGenerator>(human, kmer_size));
    
    std::vector<Variant> aligned_candidates {};
    
    for (const auto& region : some_regions) {
        std::cout << "Starting region " << region << std::endl;
        
        std::vector<AlignedRead> reads_in_region {};
        
        for (const auto& sample : samples) {
            auto&& sample_reads = a_read_manager.fetch_reads(sample, region);
            reads_in_region.insert(reads_in_region.end(),
                                   std::make_move_iterator(sample_reads.begin()),
                                   std::make_move_iterator(sample_reads.end()));
        }
        
        std::cout << "Found " << reads_in_region.size() << " reads" << std::endl;
        
        std::sort(reads_in_region.begin(), reads_in_region.end());
        
        std::vector<AlignedRead> good_reads {}, bad_reads {};
        good_reads.reserve(reads_in_region.size());
        bad_reads.reserve(reads_in_region.size());
        a_read_filter.filter_reads(std::make_move_iterator(reads_in_region.begin()),
                                   std::make_move_iterator(reads_in_region.end()),
                                   ContextBackInserter(good_reads),
                                   ContextBackInserter(bad_reads));
        good_reads.shrink_to_fit();
        bad_reads.shrink_to_fit();
        
        std::cout << "Found " << good_reads.size() << " good reads" << std::endl;
        
        a_read_transform.transform_reads(good_reads.begin(), good_reads.end());
        
        std::vector<Variant> variants_in_region {};
        
        for (const auto& read : good_reads) {
            candidate_generator.add_read(read);
        }
        
        auto candidates_in_region = candidate_generator.get_candidates(region);
        
        candidate_generator.clear();
        
        std::cout << "Found " << candidates_in_region.size() << " candidate variants" << std::endl;
        
        REQUIRE(std::is_sorted(candidates_in_region.cbegin(), candidates_in_region.cend()));
        
        for (auto& candidate : candidates_in_region) {
            aligned_candidates.emplace_back(left_align(candidate, human, a_variant_factory));
        }
    }
    
    std::cout << "Found " << aligned_candidates.size() << " total candidates" << std::endl;
}
