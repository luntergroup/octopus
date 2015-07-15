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
#include <unordered_map>
#include <cmath>

#include "test_common.h"
#include "test_utils.h"
#include "common.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "read_manager.h"
#include "allele.h"
#include "variant.h"
#include "variant_utils.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"
#include "haplotype_prior_model.h"
#include "bayesian_genotype_model.h"
#include "variational_bayes_genotype_model.h"
#include "maths.h"
#include "read_filter.h"
#include "read_filters.h"
#include "haplotype_tree.h"

using std::cout;
using std::endl;

using Octopus::ReadModel;

using Octopus::BayesianGenotypeModel::VariationalBayesGenotypeModel;
using Octopus::BayesianGenotypeModel::haplotype_population_probability;
using Octopus::BayesianGenotypeModel::posterior_predictive_probability;

TEST_CASE("haplotype posteriors sum to one", "[variational_bayes_genotype_model]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    GenomicRegion a_region {"3", 1000000, 1000001};
    
    Haplotype haplotype1 {human};
    haplotype1.push_back(Allele {a_region, "C"});
    Haplotype haplotype2 {human};
    haplotype2.push_back(Allele {a_region, "G"});
    
    Octopus::BayesianGenotypeModel::HaplotypePseudoCounts<Octopus::ProbabilityType> pseudo_counts {};
    pseudo_counts[haplotype1] = 5;
    pseudo_counts[haplotype2] = 3;
    
    unsigned ploidy {2};
    ReadModel a_read_model {ploidy};
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    Octopus::ProbabilityType haplotype_posterior_sum {};
    haplotype_posterior_sum += haplotype_population_probability(haplotype1, pseudo_counts);
    haplotype_posterior_sum += haplotype_population_probability(haplotype2, pseudo_counts);
    
    REQUIRE(is_close_to_one(haplotype_posterior_sum));
}

TEST_CASE("genotype posteriors sum to one", "[variational_bayes_genotype_model]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    GenomicRegion region1 {"3", 1000000, 1000001};
    GenomicRegion region2 {"3", 1000010, 1000011};
    
    Haplotype haplotype1 {human};
    haplotype1.push_back(Allele {region1, "A"});
    haplotype1.push_back(Allele {region2, "A"});
    Haplotype haplotype2 {human};
    haplotype2.push_back(Allele {region1, "C"});
    haplotype2.push_back(Allele {region2, "C"});
    Haplotype haplotype3 {human};
    haplotype3.push_back(Allele {region1, "G"});
    haplotype3.push_back(Allele {region2, "G"});
    Haplotype haplotype4 {human};
    haplotype4.push_back(Allele {region1, "A"});
    haplotype4.push_back(Allele {region2, "C"});
    Haplotype haplotype5 {human};
    haplotype5.push_back(Allele {region1, "C"});
    haplotype5.push_back(Allele {region2, "G"});
    Haplotype haplotype6 {human};
    haplotype6.push_back(Allele {region1, "G"});
    haplotype6.push_back(Allele {region2, "C"});
    
    std::vector<Haplotype> haplotypes {haplotype1, haplotype2, haplotype3, haplotype4, haplotype5, haplotype6};
    
    Octopus::BayesianGenotypeModel::HaplotypePseudoCounts<Octopus::ProbabilityType> pseudo_counts {};
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
    
    Octopus::ProbabilityType genotype_posterior_sum {0};
    for (const auto& genotype : genotypes) {
        genotype_posterior_sum += posterior_predictive_probability(genotype, pseudo_counts);
    }
    
    REQUIRE(is_close_to_one(genotype_posterior_sum));
    
    ploidy = 2;
    genotypes = get_all_genotypes(haplotypes, ploidy);
    ReadModel a_read_model2 {ploidy};
    VariationalBayesGenotypeModel the_model2 {a_read_model2, ploidy};
    
    genotype_posterior_sum = 0;
    for (const auto& genotype : genotypes) {
        genotype_posterior_sum += posterior_predictive_probability(genotype, pseudo_counts);
    }
    
    REQUIRE(is_close_to_one(genotype_posterior_sum));
    
    ploidy = 3;
    genotypes = get_all_genotypes(haplotypes, ploidy);
    ReadModel a_read_model3 {ploidy};
    VariationalBayesGenotypeModel the_model3 {a_read_model3, ploidy};
    
    genotype_posterior_sum = 0;
    for (const auto& genotype : genotypes) {
        genotype_posterior_sum += posterior_predictive_probability(genotype, pseudo_counts);
    }
    
    REQUIRE(is_close_to_one(genotype_posterior_sum));
}

TEST_CASE("what happens when two heterzygous genotypes are equally likely", "[variational_bayes_genotype_model]")
{
    cout.precision(10);
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    auto samples = a_read_manager.get_sample_ids();
    
    auto a_region = parse_region("4:79282976-79283139", human);
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 10));
    
    for (const auto& sample_reads : reads) {
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    unsigned ploidy {2};
    ReadModel a_read_model {ploidy};
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    Haplotype hap1 {human, a_region};
    hap1.push_back(candidates.at(0).get_reference_allele());
    hap1.push_back(candidates.at(2).get_reference_allele());
    hap1.push_back(candidates.at(3).get_reference_allele());
    
    Haplotype hap2 {human, a_region};
    hap2.push_back(candidates.at(0).get_alternative_allele());
    hap2.push_back(candidates.at(2).get_alternative_allele());
    hap2.push_back(candidates.at(3).get_alternative_allele());
    
    Haplotype hap3 {human, a_region};
    hap3.push_back(candidates.at(0).get_reference_allele());
    hap3.push_back(candidates.at(2).get_alternative_allele());
    hap3.push_back(candidates.at(3).get_alternative_allele());
    
    Haplotype hap4 {human, a_region};
    hap4.push_back(candidates.at(0).get_alternative_allele());
    hap4.push_back(candidates.at(2).get_reference_allele());
    hap4.push_back(candidates.at(3).get_reference_allele());
    
    std::vector<Haplotype> haplotypes {hap1, hap2, hap3, hap4};
    
    Genotype<Haplotype> g1 {hap1, hap2};
    Genotype<Haplotype> g2 {hap3, hap4};
    
    std::vector<Genotype<Haplotype>> gens {g1, g2};
    
//    Octopus::BayesianGenotypeModel::ReadRanges<ReadManager::SampleIdType,
//    decltype(reads)::mapped_type::iterator> read_ranges {};
//    for (const auto& sample : samples) {
//        read_ranges.emplace(sample, std::make_pair(reads[sample].begin(), reads[sample].end()));
//    }
    
    auto priors = Octopus::get_haplotype_prior_probabilities<double>(haplotypes, candidates.begin(), candidates.end());
    
//    priors[hap1] = 0.9;
//    priors[hap2] = 0.1;
//    priors[hap3] = 0.3;
//    priors[hap4] = 0.3;
    
    auto prior_counts = Octopus::BayesianGenotypeModel::get_haplotype_prior_pseudo_counts(priors, hap1, 50.0);
    
//    prior_counts[hap1] = 0.9;
//    prior_counts[hap2] = 0.1;
//    prior_counts[hap3] = 0.3;
//    prior_counts[hap4] = 0.3;
    
//    std::vector<double> p {0.9, 0.1, 0.3, 0.3};
//    //std::vector<double> p {priors[hap1], priors[hap2], priors[hap3], priors[hap4]};
//    //std::vector<double> p(4, 0.5);
//    
//    //std::vector<double> a = maximum_likelihood_dirichlet_params(p, 1.0, 100);
//    std::vector<double> a {0.9, 0.1, 0.3, 0.3};
//    
//    cout << "a: ";
//    for (auto x : a) cout << x << " ";
//    cout << endl;
//    
//    std::vector<double> z1 {1, 1, 0, 0};
//    std::vector<double> z2 {0, 0, 1, 1};
//    
////    cout << multinomial_pdf(z1, p) << endl;
////    cout << multinomial_pdf(z2, p) << endl;
//    
//    cout << "prior: " << endl;
//    cout << dirichlet_multinomial(z1, a) << endl;
//    cout << dirichlet_multinomial(z2, a) << endl;
//    cout << "diff: " << (dirichlet_multinomial(z1, a) - dirichlet_multinomial(z2, a)) << endl;
//    
//    for (auto& x : a) x += 0.5;
//    
//    cout << "posterior: " << endl;
//    cout << dirichlet_multinomial(z1, a) << endl;
//    cout << dirichlet_multinomial(z2, a) << endl;
//    cout << "diff: " << (dirichlet_multinomial(z1, a) - dirichlet_multinomial(z2, a)) << endl;
    
//    cout << "prior probabilities: " << endl;
//    cout << "hap1: " << priors[hap1] << endl;
//    cout << "hap2: " << priors[hap2] << endl;
//    cout << "hap3: " << priors[hap3] << endl;
//    cout << "hap4: " << priors[hap4] << endl;
//    cout << endl;
//    
//    cout << "prior counts: " << endl;
//    cout << "hap1: " << prior_counts[hap1] << endl;
//    cout << "hap2: " << prior_counts[hap2] << endl;
//    cout << "hap3: " << prior_counts[hap3] << endl;
//    cout << "hap4: " << prior_counts[hap4] << endl;
//    cout << endl;
//    
//    cout << "log genotype likelihoods (should be equal): " << endl;
//    cout << "g1: " << a_read_model.log_probability(reads.at(samples.front()).cbegin(), reads.at(samples.front()).cend(), g1, samples.front()) << endl;
//    cout << "g2: " << a_read_model.log_probability(reads.at(samples.front()).cbegin(), reads.at(samples.front()).cend(), g2, samples.front()) << endl;
//    cout << endl;
//    
//    cout << "log_expected_genotype_probability: " << endl;
//    cout << "g1: " << the_model.log_expected_genotype_probability(g1, prior_counts) << endl;
//    cout << "g2: " << the_model.log_expected_genotype_probability(g2, prior_counts) << endl;
//    cout << endl;
//    
//    cout << "log_rho: " << endl;
//    cout << "g1: " << the_model.log_rho(g1, prior_counts, read_ranges.at(samples.front()).first, read_ranges.at(samples.front()).second, samples.front()) << endl;
//    cout << "g2: " << the_model.log_rho(g2, prior_counts, read_ranges.at(samples.front()).first, read_ranges.at(samples.front()).second, samples.front()) << endl;
//    cout << endl;
//    
//    auto latent_posteriors = Octopus::BayesianGenotypeModel::update_latents(the_model, gens, prior_counts, read_ranges, 1);
//    
//    auto expected_occurences = Octopus::BayesianGenotypeModel::expected_haplotype_occurences(prior_counts, latent_posteriors.haplotype_pseudo_counts);
//    
//    cout << "effective occurences (1 iteration): " << endl;
//    cout << "hap1: "  << expected_occurences.at(hap1) << endl;
//    cout << "hap2: "  << expected_occurences.at(hap2) << endl;
//    cout << "hap3: "  << expected_occurences.at(hap3) << endl;
//    cout << "hap4: "  << expected_occurences.at(hap4) << endl;
//    cout << endl;
//    
//    cout << "posterior counts (1 iteration): " << endl;
//    cout << "hap1: "  << latent_posteriors.haplotype_pseudo_counts.at(hap1) << endl;
//    cout << "hap2: "  << latent_posteriors.haplotype_pseudo_counts.at(hap2) << endl;
//    cout << "hap3: "  << latent_posteriors.haplotype_pseudo_counts.at(hap3) << endl;
//    cout << "hap4: "  << latent_posteriors.haplotype_pseudo_counts.at(hap4) << endl;
//    cout << endl;
//    
//    cout << "genotype posteriors (1 iteration): " << endl;
//    cout << "g1: "  << latent_posteriors.genotype_probabilities.at(samples.front()).at(g1) << endl;
//    cout << "g2: "  << latent_posteriors.genotype_probabilities.at(samples.front()).at(g2) << endl;
//    cout << endl;
//    
//    unsigned n_iterations {100};
//    
//    latent_posteriors = Octopus::BayesianGenotypeModel::update_latents(the_model, gens, prior_counts, read_ranges, n_iterations);
//    
//    expected_occurences = Octopus::BayesianGenotypeModel::expected_haplotype_occurences(prior_counts, latent_posteriors.haplotype_pseudo_counts);
//    
//    cout << "effective occurences (" + std::to_string(n_iterations) + " iteration): " << endl;
//    cout << "hap1: "  << expected_occurences.at(hap1) << endl;
//    cout << "hap2: "  << expected_occurences.at(hap2) << endl;
//    cout << "hap3: "  << expected_occurences.at(hap3) << endl;
//    cout << "hap4: "  << expected_occurences.at(hap4) << endl;
//    cout << endl;
//    
//    cout << "posterior counts (" + std::to_string(n_iterations) + " iteration): " << endl;
//    cout << "hap1: "  << latent_posteriors.haplotype_pseudo_counts.at(hap1) << endl;
//    cout << "hap2: "  << latent_posteriors.haplotype_pseudo_counts.at(hap2) << endl;
//    cout << "hap3: "  << latent_posteriors.haplotype_pseudo_counts.at(hap3) << endl;
//    cout << "hap4: "  << latent_posteriors.haplotype_pseudo_counts.at(hap4) << endl;
//    cout << endl;
//    
//    cout << "posteriors (" + std::to_string(n_iterations) + " iteration): " << endl;
//    cout << "g1: "   << latent_posteriors.genotype_probabilities.at(samples.front()).at(g1) << endl;
//    cout << "g2: "   << latent_posteriors.genotype_probabilities.at(samples.front()).at(g2) << endl;
}

TEST_CASE("what happens when a unlikely haplotype has a large liklihood", "[variational_bayes_genotype_model]")
{
    cout.precision(10);
    
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    auto samples = a_read_manager.get_sample_ids();
    
    auto a_region = parse_region("11:67503118-67503253", human);
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 10));
    
    for (const auto& sample_reads : reads) {
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    unsigned ploidy {2};
    ReadModel a_read_model {ploidy};
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    Haplotype hap1 {human, a_region};
    hap1.push_back(candidates.at(1).get_reference_allele());
    hap1.push_back(candidates.at(2).get_reference_allele());
    
    Haplotype hap2 {human, a_region};
    hap2.push_back(candidates.at(1).get_alternative_allele());
    hap2.push_back(candidates.at(2).get_alternative_allele());
    
    std::vector<Haplotype> haplotypes {hap1, hap2};
    
    Genotype<Haplotype> g1 {hap1, hap1};
    Genotype<Haplotype> g2 {hap1, hap2};
    
    std::vector<Genotype<Haplotype>> gens {g1, g2};
    
//    Octopus::BayesianGenotypeModel::ReadRanges<ReadManager::SampleIdType,
//    decltype(reads)::mapped_type::iterator> read_ranges {};
//    for (const auto& sample : samples) {
//        read_ranges.emplace(sample, std::make_pair(reads[sample].begin(), reads[sample].end()));
//    }
    
    auto priors = Octopus::get_haplotype_prior_probabilities<double>(haplotypes, candidates.begin(), candidates.end());
    auto prior_counts = Octopus::BayesianGenotypeModel::get_haplotype_prior_pseudo_counts(priors, hap1, 2.0);
    
//    cout << "prior probabilities: " << endl;
//    cout << "hap1: " << priors[hap1] << endl;
//    cout << "hap2: " << priors[hap2] << endl;
//    cout << endl;
//    
//    cout << "prior counts: " << endl;
//    cout << "hap1: " << prior_counts[hap1] << endl;
//    cout << "hap2: " << prior_counts[hap2] << endl;
//    cout << endl;
//    
//    cout << "log genotype likelihoods: " << endl;
//    cout << "g1: " << a_read_model.log_probability(reads.at(samples.front()).cbegin(), reads.at(samples.front()).cend(), g1, samples.front()) << endl;
//    cout << "g2: " << a_read_model.log_probability(reads.at(samples.front()).cbegin(), reads.at(samples.front()).cend(), g2, samples.front()) << endl;
//    cout << endl;
//    
//    cout << "log_expected_genotype_probability: " << endl;
//    cout << "g1: " << the_model.log_expected_genotype_probability(g1, prior_counts) << endl;
//    cout << "g2: " << the_model.log_expected_genotype_probability(g2, prior_counts) << endl;
//    cout << endl;
//    
//    cout << "log_rho: " << endl;
//    cout << "g1: " << the_model.log_rho(g1, prior_counts, read_ranges.at(samples.front()).first, read_ranges.at(samples.front()).second, samples.front()) << endl;
//    cout << "g2: " << the_model.log_rho(g2, prior_counts, read_ranges.at(samples.front()).first, read_ranges.at(samples.front()).second, samples.front()) << endl;
//    cout << endl;
//    
//    auto latent_posteriors = Octopus::BayesianGenotypeModel::update_latents(the_model, gens, prior_counts, read_ranges, 1);
//    
//    auto expected_occurences = Octopus::BayesianGenotypeModel::expected_haplotype_occurences(prior_counts, latent_posteriors.haplotype_pseudo_counts);
//    
//    cout << "effective occurences (1 iteration): " << endl;
//    cout << "hap1: "  << expected_occurences.at(hap1) << endl;
//    cout << "hap2: "  << expected_occurences.at(hap2) << endl;
//    cout << endl;
//    
//    cout << "posterior counts (1 iteration): " << endl;
//    cout << "hap1: "  << latent_posteriors.haplotype_pseudo_counts.at(hap1) << endl;
//    cout << "hap2: "  << latent_posteriors.haplotype_pseudo_counts.at(hap2) << endl;
//    cout << endl;
//    
//    cout << "genotype posteriors (1 iteration): " << endl;
//    cout << "g1: "  << latent_posteriors.genotype_probabilities.at(samples.front()).at(g1) << endl;
//    cout << "g2: "  << latent_posteriors.genotype_probabilities.at(samples.front()).at(g2) << endl;
//    cout << endl;
//    
//    unsigned n_iterations {50};
//    
//    latent_posteriors = Octopus::BayesianGenotypeModel::update_latents(the_model, gens, prior_counts, read_ranges, n_iterations);
//    
//    expected_occurences = Octopus::BayesianGenotypeModel::expected_haplotype_occurences(prior_counts, latent_posteriors.haplotype_pseudo_counts);
//    
//    cout << "effective occurences (" + std::to_string(n_iterations) + " iteration): " << endl;
//    cout << "hap1: "  << expected_occurences.at(hap1) << endl;
//    cout << "hap2: "  << expected_occurences.at(hap2) << endl;
//    cout << endl;
//    
//    cout << "posterior counts (" + std::to_string(n_iterations) + " iteration): " << endl;
//    cout << "hap1: "  << latent_posteriors.haplotype_pseudo_counts.at(hap1) << endl;
//    cout << "hap2: "  << latent_posteriors.haplotype_pseudo_counts.at(hap2) << endl;
//    cout << endl;
//    
//    cout << "posteriors (" + std::to_string(n_iterations) + " iteration): " << endl;
//    cout << "g1: "   << latent_posteriors.genotype_probabilities.at(samples.front()).at(g1) << endl;
//    cout << "g2: "   << latent_posteriors.genotype_probabilities.at(samples.front()).at(g2) << endl;
}
