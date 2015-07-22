//
//  hiv_test.cpp
//  Octopus
//
//  Created by Daniel Cooke on 07/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>

#include "test_common.h"
#include "reference_genome.h"
#include "mappable_algorithms.h"
#include "reference_genome_factory.h"
#include "test_common.h"
#include "read_manager.h"
#include "read_filter.h"
#include "read_filters.h"
#include "read_utils.h"
#include "allele.h"
#include "variant.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"
#include "haplotype_prior_model.h"
#include "bayesian_genotype_model.h"
#include "variational_bayes_genotype_model.h"
#include "haplotype_phaser.h"
#include "mappable_set.h"

using std::cout;
using std::endl;

using Octopus::HaplotypeTree;
using Octopus::ReadModel;
using Octopus::HaplotypePhaser;
using Octopus::VariationalBayesGenotypeModel;

//TEST_CASE("HIV test 1)
//{
//    cout << "starting HIV test" << endl;
//    
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome hiv {a_factory.make(hiv_reference)};
//    
//    ReadManager a_read_manager(std::vector<std::string> {hiv_bam1});
//    
//    auto samples = a_read_manager.get_sample_ids();
//    
//    auto sample = samples.front();
//    
//    auto contig_name = hiv.get_contig_names().front();
//    
//    //GenomicRegion a_region {contig_name, 840, 900};
//    GenomicRegion a_region {contig_name, 0, 10000};
//    
//    auto unsorted_reads = a_read_manager.fetch_reads(sample, a_region);
//    
//    MappableSet<AlignedRead> reads {std::make_move_iterator(unsorted_reads.begin()),
//        std::make_move_iterator(unsorted_reads.end())};
//    
//    cout << "there are " << reads.size() << " reads" << endl;
//    
//    using ReadIterator = decltype(reads)::const_iterator;
//    ReadFilter<ReadIterator> a_read_filter {};
//    a_read_filter.register_filter([] (const AlignedRead& the_read) {
//        return is_good_mapping_quality(the_read, 10);
//    });
//    a_read_filter.register_filter(is_not_duplicate<ReadIterator>);
//    
//    auto good_reads = filter_reads(std::move(reads), a_read_filter).first;
//    
//    cout << "there are " << good_reads.size() << " good reads" << endl;
//    
//    CandidateVariantGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(hiv, 10));
//    candidate_generator.add_reads(good_reads.cbegin(), good_reads.cend());
//    
//    auto candidates = candidate_generator.get_candidates(a_region);
//    
//    cout << "there are " << candidates.size() << " candidates in good reads" << endl;
//    
//    cout << "min coverage before downsample " << min_coverage(good_reads, a_region) << endl;
//    cout << "max coverage before downsample " << max_coverage(good_reads, a_region) << endl;
//    cout << "mean coverage before downsample " << mean_coverage(good_reads, a_region) << endl;
//    
//    auto downsampled_reads = downsample(good_reads, 500, 200);
//    
//    cout << "there are " << downsampled_reads.size() << " downsampled reads" << endl;
//    
//    cout << "min coverage after downsample " << min_coverage(downsampled_reads, a_region) << endl;
//    cout << "max coverage after downsample " << max_coverage(downsampled_reads, a_region) << endl;
//    cout << "mean coverage after downsample " << mean_coverage(downsampled_reads, a_region) << endl;
//    
//    candidate_generator.clear();
//    candidate_generator.add_reads(downsampled_reads.cbegin(), downsampled_reads.cend());
//    
//    auto downsampled_candidates = candidate_generator.get_candidates(a_region);
//    
//    cout << "there are " << downsampled_candidates.size() << " candidates in downsampled reads" << endl;
//    
//    exit(0);
//    
////    unsigned ploidy {1};
////    ReadModel a_read_model {ploidy};
////    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
////    
////    unsigned max_haplotypes {128};
////    HaplotypePhaser phaser {hiv, the_model, ploidy, max_haplotypes};
////    
////    cout << "max genotypes: " << num_genotypes(max_haplotypes, ploidy) << endl;
////    //exit(0);
////    
////    Octopus::BayesianGenotypeModel::ReadRanges<ReadManager::SampleIdType,
////                std::move_iterator<decltype(good_reads)::mapped_type::iterator>> read_ranges {};
////    read_ranges.emplace(samples.front(), std::make_pair(std::make_move_iterator(downsampled_reads.begin()),
////                                                        std::make_move_iterator(downsampled_reads.end())));
//////    for (const auto& sample : samples) {
//////        read_ranges.emplace(sample, std::make_pair(std::make_move_iterator(good_reads[sample].begin()),
//////                                                   std::make_move_iterator(good_reads[sample].end())));
//////    }
////    
////    phaser.put_data(read_ranges, downsampled_candidates.cbegin(), downsampled_candidates.cend());
////    auto phased_regions = phaser.get_phased_regions(HaplotypePhaser::SteamingStatus::Finished);
////    
////    cout << "phased into " << phased_regions.size() << " sections" << endl;
////    
////    for (const auto& haplotype_count : phased_regions.front().the_latent_posteriors.haplotype_pseudo_counts) {
////        if (haplotype_count.second > 0.5) {
////            cout << haplotype_count.first << endl;
////            haplotype_count.first.print_explicit_alleles();
////            cout << haplotype_count.second << endl;
////        }
////    }
//}
