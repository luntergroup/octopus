//
//  ecoli_test1.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <cstddef>
#include <vector>
#include <unordered_map>

#include "test_common.h"
#include "test_utils.h"
#include "reference_genome.h"
#include "read_manager.h"
#include "allele.h"
#include "variant.h"
#include "variant_utils.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"
#include "variational_bayes_genotype_model.h"
#include "maths.h"
#include "read_filter.h"
#include "read_filters.h"
#include "haplotype_tree.h"

using std::cout;
using std::endl;

//BOOST_AUTO_TEST_CASE(ecoli test 1 : R00000042:99640-99745)
//{
//    unsigned ploidy {1};
//
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome ecoli(a_factory.make(ecoli_reference_fasta));
//
//    ReadManager a_read_manager(std::vector<std::string> {ecoli_bam});
//
//    CandidateVariantGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(ecoli, 0));
//
//    auto a_region = parse_region("R00000042:99640-99745", ecoli);
//
//    auto reference_sequence = ecoli.get_sequence(a_region);
//
//    auto sample_ids = a_read_manager.get_sample_ids();
//    auto the_sample_id = sample_ids.at(0);
//
//    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
//
//    candidate_generator.add_reads(some_reads.begin(), some_reads.end());
//
//    auto variants = candidate_generator.get_candidates(a_region);
//
//    Haplotype reference_haplotype {ecoli, a_region}; // no fully supporting read, just a read with all N's
//
//    Haplotype best_haplotype {ecoli, a_region}; // all reads fully support this
//    for (const auto& variant : variants) {
//        if (is_snp(variant)) {
//            add_to_back(variant, best_haplotype);
//        }
//    }
//
//    Haplotype okay_haplotype {ecoli, a_region}; // Bad insertion and 3 missing snps
//    add_to_back(variants[0], okay_haplotype);
//    add_to_back(variants[1], okay_haplotype);
//    add_to_back(variants[3], okay_haplotype);
//    add_to_back(variants[4], okay_haplotype);
//    add_to_back(variants[5], okay_haplotype);
//    add_to_back(variants[6], okay_haplotype);
//    add_to_back(variants[11], okay_haplotype);
//
//    unsigned num_haplotypes {3};
//    std::vector<Haplotype> haplotypes {reference_haplotype, best_haplotype, okay_haplotype};
//
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//
//    BOOST_CHECK(genotypes.size() == num_genotypes(num_haplotypes, ploidy));
//
//    ReadModel a_read_model {ploidy};
//    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
//
//    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
//    pseudo_counts[reference_haplotype] = 1;
//    pseudo_counts[best_haplotype]      = 1;
//    pseudo_counts[okay_haplotype]      = 1;
//
//    VariationalBayesGenotypeModel::SampleGenotypeResponsabilities responsabilities {};
//
//    for (const auto& genotype : genotypes) {
//        responsabilities[genotype] = the_model.genotype_responsability(genotype, some_reads.cbegin(),
//                                                                       some_reads.cend(), pseudo_counts, 0, genotypes);
//    }
//
//    std::sort(genotypes.begin(), genotypes.end(), [&responsabilities] (const auto& g1, const auto& g2) {
//        return responsabilities[g1] > responsabilities[g2];
//    });
//
//    BOOST_CHECK(genotypes.at(0).num_occurences(best_haplotype) == 1);
//
////    cout << reference_haplotype_expected_count << endl;
////    cout << best_haplotype_expected_count << endl;
////    cout << okay_haplotype_expected_count << endl;
////    cout << worst_haplotype_expected_count << endl;
//}
