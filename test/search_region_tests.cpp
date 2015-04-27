//
//  search_region_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 27/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>
#include <iterator>
#include <vector>
#include <memory>
#include <algorithm>

#include "test_common.h"
#include "reference_genome.h"
#include "region_utils.h"
#include "reference_genome_factory.h"
#include "test_common.h"
#include "read_manager.h"
#include "read_filter.h"
#include "read_filters.h"
#include "read_utils.h"
#include "allele.h"
#include "variant.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "haplotype_tree.h"
#include "read_model.h"
#include "search_regions.h"

using std::cout;
using std::endl;

//TEST_CASE("how many regions are dense?", "[search_regions]")
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager {std::vector<std::string> {human_1000g_bam1}};
//    
//    auto a_region = parse_region("3:900000-1000000", human);
//    
//    auto samples = a_read_manager.get_sample_ids();
//    
//    auto reads = a_read_manager.fetch_reads(samples, a_region);
//    
//    using ReadIterator = std::vector<AlignedRead>::const_iterator;
//    ReadFilter<ReadIterator> a_read_filter {};
//    a_read_filter.register_filter([] (const AlignedRead& the_read) {
//        return is_good_mapping_quality(the_read, 5);
//    });
//    
//    auto good_reads = filter_reads(std::move(reads), a_read_filter).first;
//    
//    CandidateVariantGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 10));
//    
//    for (auto& sample_reads : good_reads) {
//        std::sort(sample_reads.second.begin(), sample_reads.second.end());
//        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
//    }
//    
//    auto candidates = candidate_generator.get_candidates(a_region);
//    
//    unsigned num_dense {};
//    
//    for (const auto& candidate : candidates) {
//        if (is_dense_region(candidate.get_region(), good_reads, candidates, 14)) {
//            cout << candidate << endl;
//            ++num_dense;
//        }
//    }
//    
//    cout << num_dense << endl;
//}

TEST_CASE("next_sub_region's bounds are respected", "[search_regions]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    ReadManager a_read_manager {std::vector<std::string> {human_1000g_bam1}};
    
    //auto a_region = parse_region("16:9299885-9300100", human);
    auto a_region = parse_region("16:9000000-9400000", human);
    
    auto samples = a_read_manager.get_sample_ids();
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    ReadFilter<ReadIterator> a_read_filter {};
    a_read_filter.register_filter([] (const AlignedRead& the_read) {
        return is_good_mapping_quality(the_read, 5);
    });

    auto good_reads = filter_reads(std::move(reads), a_read_filter).first;
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    for (auto& sample_reads : good_reads) {
        std::sort(sample_reads.second.begin(), sample_reads.second.end());
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    auto sub_region = parse_region("16:9000000-9000000", human);
    
    sub_region = parse_region("16:9000453-9000563", human);
    auto next_region = next_sub_region(a_region, sub_region, good_reads, candidates, 5, 0);
    cout << "next sub-region " << next_region << endl;
    cout << count_overlapped(candidates.cbegin(), candidates.cend(), next_region) << endl;
//    next_region = next_sub_region(a_region, next_region, good_reads, candidates, 8, 0);
//    cout << "next sub-region " << next_region << endl;
    
//    bool bound_respected {true};
//    
//    while (ends_before(sub_region, a_region)) {
//        cout << sub_region << endl;
//        sub_region = next_sub_region(a_region, sub_region, good_reads, candidates, 5, 0);
//        
//        auto num_variants_in_sub_region = count_overlapped(candidates.cbegin(), candidates.cend(), sub_region);
//        
//        if (num_variants_in_sub_region > 5) {
//            bound_respected = false;
//            cout << num_variants_in_sub_region << endl;
//            cout << sub_region << endl;
//            break;
//        }
//    }
//    
//    REQUIRE(bound_respected);
}
