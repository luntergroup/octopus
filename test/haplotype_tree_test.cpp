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
#include "haplotype_tree.h"
#include "genotype_model.h"

TEST_CASE("haplotype_tree_single_sample_test", "[haplotype_tree]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    VariantFactory a_variant_factory {};
    
    VariantCandidateGenerator candidate_generator {};
    
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, a_variant_factory, 0));
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    unsigned ploidy {2};
    std::size_t max_num_haplotypes {10};
    double min_posterior {0.0001};
    
    GenotypeModel the_genotype_model {ploidy, max_num_haplotypes};
    
    HaplotypeTree the_haplotype_tree {human, a_read_manager, the_genotype_model, sample_ids,
                                        max_num_haplotypes, min_posterior};
    
    auto a_region = parse_region("11:1000000-1000120", human);
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(reads.cbegin(), reads.cend());
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    auto haplotypes = the_haplotype_tree.get_haplotypes(candidates);
}
