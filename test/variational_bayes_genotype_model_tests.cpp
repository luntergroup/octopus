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
