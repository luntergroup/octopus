//
//  genotype_model_tests.cpp
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
#include "variant_utils.h"
#include "variant_candidate_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"
#include "standard_genotype_model.h"

TEST_CASE("haploid_standard_genotype_model_test", "[standard_genotype_model]")
{
    unsigned ploidy {1};
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome ecoli(a_factory.make(ecoli_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {ecoli_bam});
    
    VariantCandidateGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(ecoli, 0));
    
    auto a_region = parse_region("R00000042:99640-99745", ecoli);
    
    auto reference_sequence = ecoli.get_sequence(a_region);
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(some_reads.begin(), some_reads.end());
    
    auto variants = candidate_generator.get_candidates(a_region);
    
    Haplotype reference_haplotype {ecoli};
    reference_haplotype.emplace_back(a_region, std::move(reference_sequence));
    
    Haplotype best_haplotype {ecoli};
    for (const auto& variant : variants) {
        if (is_snp(variant)) {
            add_to_back(variant, best_haplotype);
        }
    }
    
    Haplotype okay_haplotype {ecoli};
    add_to_back(variants[0], okay_haplotype);
    add_to_back(variants[1], okay_haplotype);
    add_to_back(variants[3], okay_haplotype);
    add_to_back(variants[4], okay_haplotype);
    add_to_back(variants[5], okay_haplotype);
    add_to_back(variants[6], okay_haplotype);
    add_to_back(variants[11], okay_haplotype);
    
    unsigned num_haplotypes {3};
    std::vector<Haplotype> haplotypes {reference_haplotype, best_haplotype, okay_haplotype};
    
    auto genotypes = get_all_genotypes(haplotypes, ploidy);
    
    REQUIRE(genotypes.size() == num_genotypes(num_haplotypes, ploidy));
    
    ReadModel a_read_model {ploidy};
    StandardGenotypeModel the_model {a_read_model, ploidy};
    
    
}

TEST_CASE("diploid_standard_genotype_model_test", "[standard_genotype_model]")
{
    
}
