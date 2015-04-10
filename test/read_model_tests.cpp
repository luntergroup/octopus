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
#include <unordered_map>

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
    
    VariantCandidateGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(ecoli, 0));
    
    auto a_region = parse_region("R00000042:99640-99745", ecoli);
    
    auto reference_sequence = ecoli.get_sequence(a_region);
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.at(0);
    
    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    candidate_generator.add_reads(some_reads.begin(), some_reads.end());
    
    auto variants = candidate_generator.get_candidates(a_region);
    
    Haplotype reference_haplotype {ecoli, a_region}; // single fully supporting read
    
    Haplotype best_haplotype {ecoli}; // most reads fully support this
    for (const auto& variant : variants) {
        if (is_snp(variant)) {
            add_to_back(variant, best_haplotype);
        }
    }
    
    Haplotype okay_haplotype {ecoli}; // Bad insertion and 3 missing snps
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
    
    ReadModel the_model {ploidy};
    unsigned sample {0};
    
    std::unordered_map<Genotype, double> genotype_log_probabilities {};
    
    for (const auto& genotype : genotypes) {
        genotype_log_probabilities[genotype] = the_model.log_probability(some_reads, genotype, sample);
    }
    
    std::sort(genotypes.begin(), genotypes.end(), [&genotype_log_probabilities] (const auto& g1, const auto& g2) {
        return genotype_log_probabilities[g1] > genotype_log_probabilities[g2];
    });
    
    REQUIRE(genotypes.at(0).at(0) == best_haplotype);
    REQUIRE(genotypes.at(1).at(0) == okay_haplotype);
    REQUIRE(genotypes.at(2).at(0) == reference_haplotype);
}

TEST_CASE("diploid_read_model_test", "[read_model]")
{
    unsigned ploidy {2};
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    VariantCandidateGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
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
    add_to_back(variants[0], hap1); // well supported insert
    add_to_back(variants[2], hap1); // well supported snp
    
    Haplotype hap2 {human, a_region};
    add_to_back(variants[1], hap2); // this is a low quality snp
    
    Haplotype hap3 {human, a_region};
    add_to_back(variants[0], hap3);
    add_to_back(variants[1], hap3);
    add_to_back(variants[2], hap3);
    
    unsigned num_haplotypes {4};
    std::vector<Haplotype> haplotypes {reference_haplotype, hap1, hap2, hap3};
    
    ReadModel the_model {ploidy};
    unsigned sample {0};
    
//    const auto& hap2_supporting_read =  some_reads.at(1);
//    auto ref_log_prob  = the_model.log_probability(hap2_supporting_read, reference_haplotype, sample);
//    auto hap2_log_prob = the_model.log_probability(hap2_supporting_read, hap2, sample);
    
    auto genotypes = get_all_genotypes(haplotypes, ploidy);
    
    std::unordered_map<Genotype, double> genotype_log_probabilities {};
    
    for (const auto& genotype : genotypes) {
        genotype_log_probabilities[genotype] = the_model.log_probability(some_reads, genotype, sample);
    }
    
    std::sort(genotypes.begin(), genotypes.end(), [&genotype_log_probabilities] (const auto& g1, const auto& g2) {
        return genotype_log_probabilities[g1] > genotype_log_probabilities[g2];
    });
    
//    for (const auto& genotype : genotypes) {
//        std::cout << genotype << " " << genotype_log_probabilities[genotype] << std::endl;
//    }
    
    REQUIRE(genotypes.at(0).num_occurences(hap1) == 1);
    REQUIRE(genotypes.at(0).num_occurences(hap2) == 1);
    
    REQUIRE(genotypes.at(1).num_occurences(hap1) == 1);
    REQUIRE(genotypes.at(1).num_occurences(reference_haplotype) == 1);
}

TEST_CASE("two_sample_diploid_read_model_test", "[read_model]")
{
    unsigned ploidy {2};
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1, human_1000g_bam2});
    
    VariantCandidateGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    auto a_region = parse_region("2:104142870-104142884", human);
    
    auto reference_sequence = human.get_sequence(a_region);
    
    auto sample_ids = a_read_manager.get_sample_ids();
    
    auto some_reads = a_read_manager.fetch_reads(sample_ids, a_region);
    
    candidate_generator.add_reads(some_reads[sample_ids[0]].cbegin(), some_reads[sample_ids[0]].cend());
    candidate_generator.add_reads(some_reads[sample_ids[1]].cbegin(), some_reads[sample_ids[1]].cend());
    
    auto variants = candidate_generator.get_candidates(a_region);
    
    REQUIRE(variants.size() == 3);
    
    Haplotype reference_haplotype {human, a_region}; // there are no reads completely supporting the reference
    
    Haplotype hap1 {human, a_region};
    add_to_back(variants[0], hap1); // well supported insert
    add_to_back(variants[2], hap1); // well supported snp
    
    Haplotype hap2 {human, a_region};
    add_to_back(variants[1], hap2); // this is a low quality snp
    
    Haplotype hap3 {human, a_region};
    add_to_back(variants[0], hap3);
    add_to_back(variants[1], hap3);
    add_to_back(variants[2], hap3);
    
    unsigned num_haplotypes {4};
    std::vector<Haplotype> haplotypes {reference_haplotype, hap1, hap2, hap3};
    
    auto genotypes = get_all_genotypes(haplotypes, ploidy);
    
    REQUIRE(genotypes.size() == num_genotypes(num_haplotypes, ploidy));
    
    ReadModel a_read_model {ploidy};
    
    std::unordered_map<Genotype, double> genotype_log_probabilities0 {};
    
    for (const auto& genotype : genotypes) {
        genotype_log_probabilities0[genotype] = a_read_model.log_probability(some_reads[sample_ids[0]], genotype, 0);
    }
    
    std::unordered_map<Genotype, double> genotype_log_probabilities1 {};
    
    for (const auto& genotype : genotypes) {
        genotype_log_probabilities1[genotype] = a_read_model.log_probability(some_reads[sample_ids[1]], genotype, 1);
    }
    
//    for (const auto& genotype : genotypes) {
//        std::cout << genotype << " " << genotype_log_probabilities0[genotype] << " " << genotype_log_probabilities1[genotype] << std::endl;
//    }
}
