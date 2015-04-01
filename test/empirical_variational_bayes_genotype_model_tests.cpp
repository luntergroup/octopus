//
//  empirical_variational_bayes_genotype_model_tests.cpp
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
#include "empirical_variational_bayes_genotype_model.h"

TEST_CASE("haploid_empirical_variational_bayes_genotype_model", "[empirical_variational_bayes_genotype_model]")
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
    
    REQUIRE(variants.size() == 14);
    
    Haplotype reference_haplotype {ecoli};
    reference_haplotype.emplace_back(a_region, std::move(reference_sequence));
    
    Haplotype best_haplotype {ecoli};
    for (const auto& variant : variants) {
        if (is_snp(variant)) {
            best_haplotype.emplace_back(variant);
        }
    }
    
    Haplotype okay_haplotype {ecoli};
    okay_haplotype.emplace_back(variants[1]); okay_haplotype.emplace_back(variants[3]);
    okay_haplotype.emplace_back(variants[5]); okay_haplotype.emplace_back(variants[8]);
    okay_haplotype.emplace_back(variants[11]); okay_haplotype.emplace_back(variants[12]);
    
    Haplotype worst_haplotype {ecoli};
    for (const auto& variant : variants) {
        if (is_indel(variant)) {
            worst_haplotype.emplace_back(variant);
        }
    }
    
    unsigned num_haplotypes {4};
    std::vector<Haplotype> haplotypes {reference_haplotype, best_haplotype, okay_haplotype, worst_haplotype};
    
    auto genotypes = get_all_genotypes(haplotypes, ploidy);
    
    REQUIRE(genotypes.size() == num_genotypes(num_haplotypes, ploidy));
    
    ReadModel a_read_model {ploidy};
    EmpiricalVariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    EmpiricalVariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
    pseudo_counts[reference_haplotype] = 2;
    pseudo_counts[best_haplotype]      = 1;
    pseudo_counts[okay_haplotype]      = 1;
    pseudo_counts[worst_haplotype]     = 1;
    
    EmpiricalVariationalBayesGenotypeModel::SampleGenotypeResponsabilities responsabilities {};
    
    for (const auto& genotype : genotypes) {
        responsabilities[genotype] = the_model.genotype_responsability(genotype, some_reads, pseudo_counts, 0, genotypes);
    }
    
    auto reference_haplotype_expected_count = the_model.expected_haplotype_count(reference_haplotype, responsabilities);
    auto best_haplotype_expected_count      = the_model.expected_haplotype_count(best_haplotype, responsabilities);
    auto okay_haplotype_expected_count      = the_model.expected_haplotype_count(okay_haplotype, responsabilities);
    auto worst_haplotype_expected_count     = the_model.expected_haplotype_count(worst_haplotype, responsabilities);
    
    std::cout << reference_haplotype_expected_count << std::endl;
    std::cout << best_haplotype_expected_count << std::endl;
    std::cout << okay_haplotype_expected_count << std::endl;
    std::cout << worst_haplotype_expected_count << std::endl;
}
