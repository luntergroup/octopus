//
//  read_model_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
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

TEST_CASE("haploid_read_model_test", "[read_model]")
{
    unsigned ploidy {1};
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome ecoli(a_factory.make(ecoli_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {ecoli_bam});
    
    VariantFactory a_variant_factory {};
    VariantCandidateGenerator candidate_generator {};
    candidate_generator.register_generator(
        std::make_unique<AlignmentCandidateVariantGenerator>(ecoli, a_variant_factory, 0));
    
    auto a_region = parse_region("R00000042:99640-99745", ecoli);
    
    auto reference_sequence = ecoli.get_sequence(a_region);
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(some_reads.begin(), some_reads.end());
    
    auto variants = candidate_generator.get_candidates(a_region);
    
    Haplotype reference_haplotype {ecoli, a_region}; // single fully supporting read
    
    Haplotype best_haplotype {ecoli, a_region}; // most reads fully support this
    for (const auto& variant : variants) {
        if (is_snp(variant)) {
            best_haplotype.emplace_back(variant);
        }
    }
    
    Haplotype okay_haplotype {ecoli, a_region}; // Bad insertion and 3 missing snps
    okay_haplotype.emplace_back(variants[0]); okay_haplotype.emplace_back(variants[1]);
    okay_haplotype.emplace_back(variants[3]); okay_haplotype.emplace_back(variants[4]);
    okay_haplotype.emplace_back(variants[5]); okay_haplotype.emplace_back(variants[6]);
    okay_haplotype.emplace_back(variants[11]);
    
    Haplotype worst_haplotype {ecoli, a_region}; // no fully supporting reads
    for (const auto& variant : variants) {
        if (is_indel(variant)) {
            worst_haplotype.emplace_back(variant);
        }
    }
    
    unsigned num_haplotypes {4};
    std::vector<Haplotype> haplotypes {reference_haplotype, best_haplotype, okay_haplotype, worst_haplotype};
    
    auto genotypes = get_all_genotypes(haplotypes, ploidy);
    
    REQUIRE(genotypes.size() == num_genotypes(num_haplotypes, ploidy));
    
    ReadModel the_model {ploidy};
    
    auto ref_hap_log_prob   = the_model.log_probability(some_reads, haplotypes[0]);
    auto best_hap_log_prob  = the_model.log_probability(some_reads, haplotypes[1]);
    auto okay_hap_log_prob  = the_model.log_probability(some_reads, haplotypes[2]);
    auto worst_hap_log_prob = the_model.log_probability(some_reads, haplotypes[3]);
    
    REQUIRE(best_hap_log_prob > ref_hap_log_prob);
    REQUIRE(best_hap_log_prob > okay_hap_log_prob);
    REQUIRE(best_hap_log_prob > worst_hap_log_prob);
    REQUIRE(okay_hap_log_prob > ref_hap_log_prob);
    REQUIRE(okay_hap_log_prob > worst_hap_log_prob);
}

TEST_CASE("diploid_read_model_test", "[read_model]")
{
    unsigned ploidy {2};
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    VariantFactory a_variant_factory {};
    VariantCandidateGenerator candidate_generator {};
    candidate_generator.register_generator(
                std::make_unique<AlignmentCandidateVariantGenerator>(human, a_variant_factory, 0));
    
    auto a_region = parse_region("2:104142870-104142884", human);
    
    auto reference_sequence = human.get_sequence(a_region);
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(some_reads.begin(), some_reads.end());
    
    auto variants = candidate_generator.get_candidates(a_region);
    
    REQUIRE(variants.size() == 3);
    
    Haplotype reference_haplotype {human, a_region}; // there are no reads completely supporting the reference
    
    Haplotype hap1 {human, a_region};
    hap1.emplace_back(variants[0]); // high quality insert
    hap1.emplace_back(variants[2]); // high quality snp
    
    Haplotype hap2 {human, a_region};
    hap2.emplace_back(variants[1]); // this is a low quality snp
    
    Haplotype hap3 {human, a_region};
    hap3.emplace_back(variants[0]);
    hap3.emplace_back(variants[1]);
    hap3.emplace_back(variants[2]);
    
    unsigned num_haplotypes {4};
    std::vector<Haplotype> haplotypes {reference_haplotype, hap1, hap2, hap3};
    
    ReadModel the_model {ploidy};
    
    auto ref_log_prob  = the_model.log_probability(some_reads, haplotypes[0]);
    auto hap1_log_prob = the_model.log_probability(some_reads, haplotypes[1]);
    auto hap2_log_prob = the_model.log_probability(some_reads, haplotypes[2]);
    auto hap3_log_prob = the_model.log_probability(some_reads, haplotypes[3]);
    
    REQUIRE(ref_log_prob > hap2_log_prob); // As the snp in hap2 is not good evidence
    REQUIRE(hap1_log_prob > ref_log_prob);
    REQUIRE(hap1_log_prob > hap2_log_prob);
    REQUIRE(hap1_log_prob > hap3_log_prob);
    REQUIRE(hap3_log_prob > ref_log_prob);
    REQUIRE(hap3_log_prob > hap2_log_prob);
}
