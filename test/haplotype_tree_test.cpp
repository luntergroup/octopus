//
//  haplotype_tree_test.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>
#include <cstddef>

#include "test_common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "read_manager.h"
#include "variant.h"
#include "variant_candidate_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "variational_bayes_genotype_model.h"
#include "haplotype_tree.h"
#include "region_utils.h"
#include "read_filter.h"
#include "read_filters.h"

using std::cout;
using std::endl;

TEST_CASE("haplotype_extension_test", "[haplotype_tree")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager {std::vector<std::string> {human_1000g_bam1}};
    
    VariantCandidateGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto a_region = parse_region("13:70634415-70635400", human);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    std::sort(reads.begin(), reads.end());
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
}

TEST_CASE("haplotype_tree_single_sample_test", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    VariantCandidateGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto a_region = parse_region("16:78203360-78203450", human);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    ReadFilter<ReadIterator> a_read_filter {};
    a_read_filter.register_filter(is_mapped);
    
    std::vector<AlignedRead> good_reads {}, bad_reads {};
    good_reads.reserve(reads.size());
    bad_reads.reserve(reads.size());
    a_read_filter.filter_reads(std::make_move_iterator(reads.begin()), std::make_move_iterator(reads.end()),
                               ContextBackInserter(good_reads), ContextBackInserter(bad_reads));
    
    candidate_generator.add_reads(good_reads.cbegin(), good_reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    HaplotypeTree haplotype_tree {human};
    
    for (const auto& candidate : candidates) {
        //cout << candidate << endl;
        //haplotype_tree.extend_haplotypes(candidate);
    }
    
    Allele allele1 {parse_region("1:100-103", human), "AGA"};
    Allele allele2 {parse_region("1:101-102", human), "C"};
    
    cout << contains(allele1, allele2) << endl;
}
