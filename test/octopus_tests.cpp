//
//  octopus_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <iterator>
#include <vector>
#include <memory>
#include <algorithm>

#include "test_common.h"
#include "reference_genome.h"
#include "mappable_algorithms.h"
#include "test_common.h"
#include "read_manager.h"
#include "mock_objects.h"
#include "read_filter.h"
#include "read_filters.h"
#include "read_transform.h"
#include "read_transformations.h"
#include "allele.h"
#include "variant.h"
#include "candidate_variant_generator.h"
#include "variant_utils.h"
#include "alignment_candidate_variant_generator.h"
#include "assembler_candidate_variant_generator.h"
#include "external_variant_candidates.h"
#include "haplotype.h"
#include "genotype.h"
#include "haplotype_tree.h"
#include "read_model.h"
#include "variational_bayes_genotype_model.h"
#include "variant_file_factory.h"
#include "variant_file_reader.h"
#include "octopus.h"

using std::cout;
using std::endl;

//BOOST_AUTO_TEST_CASE(can call in complex region)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager {std::vector<std::string> {human_1000g_bam3}};
//    
//    //auto a_region = parse_region("16:9299900-9300055", human);
//    //auto a_region = parse_region("16:9299850-9299970", human);
//    //auto a_region = parse_region("16:9300000-9300100", human);
//    //auto a_region = parse_region("16:9299900-9300038", human);
//    
//    auto samples = a_read_manager.get_sample_ids();
//    
//    using ReadIterator = std::vector<AlignedRead>::const_iterator;
//    ReadFilter<ReadIterator> a_read_filter {};
//    a_read_filter.register_filter(is_not_secondary_alignment);
//    a_read_filter.register_filter([] (const AlignedRead& the_read) {
//        return is_good_mapping_quality(the_read, 20);
//    });
//    a_read_filter.register_filter([] (const AlignedRead& the_read) {
//        return has_sufficient_good_quality_bases(the_read, 20, 10);
//    });
//    a_read_filter.register_filter(is_not_duplicate<ReadIterator>);
//    
//    ReadTransform a_read_transform {};
//    a_read_transform.register_transform(trim_adapters);
//    a_read_transform.register_transform(trim_soft_clipped);
//    
//    VariantCandidateGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//    
//    auto reads = a_read_manager.fetch_reads(samples[0], a_region);
//    
//    cout << "there are " << reads.size() << " reads" << endl;
//    
//    std::sort(reads.begin(), reads.end());
//    
//    std::vector<AlignedRead> good_reads {}, bad_reads {};
//    good_reads.reserve(reads.size());
//    bad_reads.reserve(reads.size());
//    a_read_filter.filter_reads(std::make_move_iterator(reads.begin()), std::make_move_iterator(reads.end()),
//                               ContextBackInserter(good_reads), ContextBackInserter(bad_reads));
//    reads.clear();
//    good_reads.shrink_to_fit();
//    bad_reads.shrink_to_fit();
//    
//    cout << "there are " << good_reads.size() << " good reads" << endl;
//    
//    a_read_transform.transform_reads(good_reads.begin(), good_reads.end());
//    
//    candidate_generator.add_reads(good_reads.cbegin(), good_reads.cend());
//    
//    auto candidates = candidate_generator.get_candidates(a_region);
//    
//    candidate_generator.clear();
//    
//    for (auto& candidate : candidates) {
//        left_align(candidate, human);
//    }
//    
//    cout << "there are " << candidates.size() << " candidates" << endl;
//    
//    HaplotypeTree haplotype_tree {human};
//    
//    for (const auto& candidate : candidates) {
//        haplotype_tree.extend(candidate.get_reference_allele());
//        haplotype_tree.extend(candidate.get_alternative_allele());
//    }
//    
//    auto haplotypes = haplotype_tree.get_haplotypes(a_region);
//    
//    cout << "there are " << haplotypes.size() << " haplotypes" << endl;
//    
//    unique_least_complex(haplotypes);
//    
//    for (const auto& haplotype : haplotypes) {
//        haplotype_tree.prune_unique(haplotype);
//    }
//    
//    cout << "there are " << haplotypes.size() << " unique haplotypes" << endl;
//    
//    Haplotype reference_haplotype {human, a_region};
//    
//    VariationalBayesGenotypeModel::HaplotypePseudoCounts pseudo_counts {};
//    
//    for (const auto& haplotype : haplotypes) {
//        if (haplotype == reference_haplotype) {
//            pseudo_counts[haplotype] = 1e6;
//        } else {
//            pseudo_counts[haplotype] = 15.0;
//        }
//    }
//    
//    unsigned ploidy {2};
//    
//    auto genotypes = get_all_genotypes(haplotypes, ploidy);
//    
//    cout << "there are " << genotypes.size() << " genotypes" << endl;
//    
//    ReadModel a_read_model {ploidy};
//    
//    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
//    
//    SamplesReads the_reads {};
//    the_reads.push_back({good_reads.cbegin(), good_reads.cend()});
//    
//    auto results = update_parameters(the_model, genotypes, pseudo_counts, the_reads, 3);
//    auto responsabilities        = results.first;
//    auto posterior_pseudo_counts = results.second;
//    
//    auto& sample_responsabilities = responsabilities[0];
//    
//    std::sort(genotypes.begin(), genotypes.end(), [&sample_responsabilities] (const auto& g1, const auto& g2) {
//        return sample_responsabilities[g1] > sample_responsabilities[g2];
//    });
//    
//    cout << genotypes[0] << endl;
//    
//    for (const auto& variant : candidates) {
//        cout << variant << " "
//        << the_model.posterior_probability_allele_in_sample(variant.get_reference_allele(), haplotypes,
//                                                            sample_responsabilities, genotypes)
//        << " "
//        << the_model.posterior_probability_allele_in_sample(variant.get_alternative_allele(), haplotypes,
//                                                            sample_responsabilities, genotypes)
//        << endl;
//    }
//}