//
//  variational_bayes_genotype_model_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>
#include <cstddef>
#include <vector>
#include <cmath>

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
#include "variational_bayes_genotype_model.h"
#include "maths.h"

using std::cout;
using std::endl;

bool is_close_to_one(double val)
{
    return std::abs(val - 1) < 0.0000000000001;
}

TEST_CASE("haplotype_normalisation_test", "[variational_bayes_genotype_model]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    GenomicRegion a_region {"3", 1000000, 1000001};
    
    Haplotype haplotype1 {human};
    haplotype1.emplace_back(a_region, "C");
    Haplotype haplotype2 {human};
    haplotype2.emplace_back(a_region, "G");
    
    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
    pseudo_counts[haplotype1] = 5;
    pseudo_counts[haplotype2] = 3;
    
    unsigned ploidy {2};
    ReadModel a_read_model {ploidy};
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    double haplotype_posterior_sum {};
    haplotype_posterior_sum += the_model.posterior_haplotype_probability(haplotype1, pseudo_counts);
    haplotype_posterior_sum += the_model.posterior_haplotype_probability(haplotype2, pseudo_counts);
    
    REQUIRE(is_close_to_one(haplotype_posterior_sum));
}

TEST_CASE("genotype_normalisation_test", "[variational_bayes_genotype_model]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    GenomicRegion region1 {"3", 1000000, 1000001};
    GenomicRegion region2 {"3", 1000010, 1000011};
    
    Haplotype haplotype1 {human};
    haplotype1.emplace_back(region1, "A");
    haplotype1.emplace_back(region2, "A");
    Haplotype haplotype2 {human};
    haplotype2.emplace_back(region1, "C");
    haplotype2.emplace_back(region2, "C");
    Haplotype haplotype3 {human};
    haplotype3.emplace_back(region1, "G");
    haplotype3.emplace_back(region2, "G");
    Haplotype haplotype4 {human};
    haplotype4.emplace_back(region1, "A");
    haplotype4.emplace_back(region2, "C");
    Haplotype haplotype5 {human};
    haplotype5.emplace_back(region1, "C");
    haplotype5.emplace_back(region2, "G");
    Haplotype haplotype6 {human};
    haplotype6.emplace_back(region1, "G");
    haplotype6.emplace_back(region2, "C");
    
    std::vector<Haplotype> haplotypes {haplotype1, haplotype2, haplotype3, haplotype4, haplotype5, haplotype6};
    
    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
    pseudo_counts[haplotype1] = 1;
    pseudo_counts[haplotype2] = 3.2;
    pseudo_counts[haplotype3] = 2;
    pseudo_counts[haplotype4] = 1.5;
    pseudo_counts[haplotype5] = 5.6;
    pseudo_counts[haplotype6] = 1.1;

    unsigned ploidy {1};
    ReadModel a_read_model1 {ploidy};
    VariationalBayesGenotypeModel the_model1 {a_read_model1, ploidy};
    
    auto genotypes = get_all_genotypes(haplotypes, ploidy);
    
    double genotype_posterior_sum {0};
    for (const auto& genotype : genotypes) {
        genotype_posterior_sum += the_model1.posterior_predictive_probability(genotype, pseudo_counts);
    }
    
    REQUIRE(is_close_to_one(genotype_posterior_sum));
    
    ploidy = 2;
    genotypes = get_all_genotypes(haplotypes, ploidy);
    ReadModel a_read_model2 {ploidy};
    VariationalBayesGenotypeModel the_model2 {a_read_model2, ploidy};
    
    genotype_posterior_sum = 0;
    for (const auto& genotype : genotypes) {
        genotype_posterior_sum += the_model2.posterior_predictive_probability(genotype, pseudo_counts);
    }
    
    REQUIRE(is_close_to_one(genotype_posterior_sum));
    
    ploidy = 3;
    genotypes = get_all_genotypes(haplotypes, ploidy);
    ReadModel a_read_model3 {ploidy};
    VariationalBayesGenotypeModel the_model3 {a_read_model3, ploidy};
    
    genotype_posterior_sum = 0;
    for (const auto& genotype : genotypes) {
        genotype_posterior_sum += the_model3.posterior_predictive_probability(genotype, pseudo_counts);
    }
    
    REQUIRE(is_close_to_one(genotype_posterior_sum));
}

TEST_CASE("single_sample_haploid_variational_bayes_genotype_model", "[variational_bayes_genotype_model]")
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
    
    ReadModel a_read_model {ploidy};
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
    pseudo_counts[reference_haplotype] = 2;
    pseudo_counts[best_haplotype]      = 1;
    pseudo_counts[okay_haplotype]      = 1;
    
    VariationalBayesGenotypeModel::SampleGenotypeResponsabilities responsabilities {};
    
    for (const auto& genotype : genotypes) {
        responsabilities[genotype] = the_model.genotype_responsability(genotype, some_reads, pseudo_counts, 0, genotypes);
    }
    
    std::sort(genotypes.begin(), genotypes.end(), [&responsabilities] (const auto& g1, const auto& g2) {
        return responsabilities[g1] > responsabilities[g2];
    });
    
    REQUIRE(genotypes.at(0).num_occurences(best_haplotype) == 1);
    
//    cout << reference_haplotype_expected_count << endl;
//    cout << best_haplotype_expected_count << endl;
//    cout << okay_haplotype_expected_count << endl;
//    cout << worst_haplotype_expected_count << endl;
}

TEST_CASE("single_sample_diploid_variational_bayes_genotype_model", "[variational_bayes_genotype_model]")
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
    
    auto genotypes = get_all_genotypes(haplotypes, ploidy);
    
    REQUIRE(genotypes.size() == num_genotypes(num_haplotypes, ploidy));
    
    ReadModel a_read_model {ploidy};
    
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
    pseudo_counts[reference_haplotype] = 400; // set so we accept snps with qual > 21
    pseudo_counts[hap1]                = 1;
    pseudo_counts[hap2]                = 1;
    pseudo_counts[hap3]                = 1;
    
    VariationalBayesGenotypeModel::SampleGenotypeResponsabilities responsabilities {};
    VariationalBayesGenotypeModel::HaplotypePseudoCounts posterior_pseudo_counts {pseudo_counts};
    
    for (unsigned i {}; i < 10; ++i) {
        for (const auto& genotype : genotypes) {
            responsabilities[genotype] = the_model.genotype_responsability(genotype, some_reads,
                                                                           posterior_pseudo_counts, 0, genotypes);
        }
        
        posterior_pseudo_counts[reference_haplotype] = the_model.posterior_haplotype_pseudo_count(reference_haplotype,
                                                                              pseudo_counts[reference_haplotype],
                                                                              { responsabilities });
        posterior_pseudo_counts[hap1] = the_model.posterior_haplotype_pseudo_count(hap1, pseudo_counts[hap1],
                                                                              { responsabilities });
        posterior_pseudo_counts[hap2] = the_model.posterior_haplotype_pseudo_count(hap2, pseudo_counts[hap2],
                                                                              { responsabilities });
        posterior_pseudo_counts[hap3] = the_model.posterior_haplotype_pseudo_count(hap3, pseudo_counts[hap3],
                                                                              { responsabilities });
    }
    
    std::sort(genotypes.begin(), genotypes.end(), [&responsabilities] (const auto& g1, const auto& g2) {
        return responsabilities[g1] > responsabilities[g2];
    });
    
//    for (auto& g : genotypes) {
//        cout << g << " " << responsabilities[g] << endl;
//    }
    
    REQUIRE(genotypes.at(0).num_occurences(hap1) == 1);
    REQUIRE(genotypes.at(0).num_occurences(reference_haplotype) == 1);
    
    REQUIRE(genotypes.at(1).num_occurences(hap1) == 1);
    REQUIRE(genotypes.at(1).num_occurences(hap2) == 1);
}

TEST_CASE("two_samples_diploid_variational_bayes_genotype_model1", "[variational_bayes_genotype_model]")
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
    
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
    pseudo_counts[reference_haplotype] = 1;
    pseudo_counts[hap1]                = 1;
    pseudo_counts[hap2]                = 1;
    pseudo_counts[hap3]                = 1;
    
    SamplesReads the_reads {};
    the_reads.emplace_back(std::move(some_reads[sample_ids[0]]));
    the_reads.emplace_back(std::move(some_reads[sample_ids[1]]));
    
    auto results = update_parameters(the_model, genotypes, pseudo_counts, the_reads, 10);
    auto responsabilities        = results.first;
    auto posterior_pseudo_counts = results.second;
    
//    for (auto& p : responsabilities[0]) {
//        cout << p.first << " " << p.second << endl;
//    }
//    cout << endl;
//    for (auto& p : responsabilities[1]) {
//        cout << p.first << " " << p.second << endl;
//    }
//    cout << endl;
//    for (const auto& c : posterior_pseudo_counts) {
//        cout << c.first << " " << c.second << endl;
//    }
}

TEST_CASE("two_samples_diploid_variational_bayes_genotype_model2", "[variational_bayes_genotype_model]")
{
    unsigned ploidy {2};
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1, human_1000g_bam3});
    
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
    
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
    pseudo_counts[reference_haplotype] = 1;
    pseudo_counts[hap1]                = 1;
    pseudo_counts[hap2]                = 1;
    pseudo_counts[hap3]                = 1;
    
    SamplesReads the_reads {};
    the_reads.emplace_back(std::move(some_reads[sample_ids[0]]));
    the_reads.emplace_back(std::move(some_reads[sample_ids[1]]));
    
    auto results = update_parameters(the_model, genotypes, pseudo_counts, the_reads, 10);
    auto responsabilities        = results.first;
    auto posterior_pseudo_counts = results.second;
    
//    for (auto& p : responsabilities[0]) {
//        cout << p.first << " " << p.second << endl;
//    }
//    cout << endl;
//    for (auto& p : responsabilities[1]) {
//        cout << p.first << " " << p.second << endl;
//    }
//    cout << endl;
//    for (const auto& c : posterior_pseudo_counts) {
//        cout << c.first << " " << c.second << endl;
//    }
}

TEST_CASE("three_samples_diploid_variational_bayes_genotype_model", "[variational_bayes_genotype_model]")
{
    unsigned ploidy {2};
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1, human_1000g_bam2, human_1000g_bam3});
    
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
    
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
    pseudo_counts[reference_haplotype] = 1;
    pseudo_counts[hap1]                = 1;
    pseudo_counts[hap2]                = 1;
    pseudo_counts[hap3]                = 1;
    
    SamplesReads the_reads {};
    the_reads.emplace_back(std::move(some_reads[sample_ids[0]]));
    the_reads.emplace_back(std::move(some_reads[sample_ids[1]]));
    the_reads.emplace_back(std::move(some_reads[sample_ids[2]]));
    
    auto results = update_parameters(the_model, genotypes, pseudo_counts, the_reads, 10);
    auto responsabilities        = results.first;
    auto posterior_pseudo_counts = results.second;
    
//    for (auto& p : responsabilities[0]) {
//        cout << p.first << " " << p.second << endl;
//    }
//    cout << endl;
//    for (auto& p : responsabilities[1]) {
//        cout << p.first << " " << p.second << endl;
//    }
//    cout << endl;
//    for (auto& p : responsabilities[2]) {
//        cout << p.first << " " << p.second << endl;
//    }
}
