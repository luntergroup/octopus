//
//  variant_caller_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>

#include "test_common.h"
#include "reference_genome.h"
#include "region_algorithms.h"
#include "reference_genome_factory.h"
#include "test_common.h"
#include "read_manager.h"
#include "read_filter.h"
#include "read_filters.h"
#include "read_utils.h"
#include "allele.h"
#include "variant.h"
#include "variant_utils.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"
#include "variational_bayes_genotype_model.h"
#include "haplotype_phaser.h"
#include "variant_caller.h"

#include "maths.h"

using std::cout;
using std::endl;

using Octopus::ReadModel;
using Octopus::ReadModel;
using Octopus::VariationalBayesGenotypeModel;
using Octopus::HaplotypePhaser;

TEST_CASE("can call", "variant_caller")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    auto samples = a_read_manager.get_sample_ids();
    
    //auto a_region = parse_region("16:62646800-62647150", human);
    auto a_region = parse_region("11:67503118-67503253", human);
    //auto a_region = parse_region("11:67503125-67503164", human);
    //auto a_region = parse_region("11:67502976-67503491", human);
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    ReadFilter<ReadIterator> a_read_filter {};
    a_read_filter.register_filter([] (const AlignedRead& the_read) {
        return is_good_mapping_quality(the_read, 10);
    });
    
    auto good_reads = filter_reads(std::move(reads), a_read_filter).first;
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 10));
    
    for (const auto& sample_reads : good_reads) {
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    unsigned ploidy {2};
    ReadModel a_read_model {ploidy};
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    unsigned max_haplotypes {64};
    HaplotypePhaser phaser {human, the_model, ploidy, max_haplotypes};
    
    Octopus::BayesianGenotypeModel::ReadRanges<ReadManager::SampleIdType,
        std::move_iterator<decltype(good_reads)::mapped_type::iterator>> read_ranges {};
    for (const auto& sample : samples) {
        read_ranges.emplace(sample, std::make_pair(std::make_move_iterator(good_reads[sample].begin()),
                                                   std::make_move_iterator(good_reads[sample].end())));
    }
    
    phaser.put_data(read_ranges, candidates.cbegin(), candidates.cend());

    auto phased_regions = phaser.get_phased_regions(true);
    
    auto allele_posteriors = Octopus::VariantCaller::get_allele_posteriors(samples, phased_regions, candidates);
    
    for (const auto& haplotype_count : phased_regions.front().the_latent_posteriors.haplotype_pseudo_counts) {
        if (haplotype_count.second > 1.0) {
            cout << haplotype_count.first << endl;
            haplotype_count.first.print_explicit_alleles();
            cout << haplotype_count.second << endl;
        }
    }
    
    for (const auto& sample : samples) {
        for (const auto& allele_posterior : allele_posteriors.at(sample)) {
            cout << sample << " " << allele_posterior.first << " " << allele_posterior.second << endl;
        }
    }
}

TEST_CASE("can call2", "variant_caller")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    auto samples = a_read_manager.get_sample_ids();
    
    //auto a_region = parse_region("6:29909020-29909269", human);
    //auto a_region = parse_region("6:29909330-29909670", human);
    //auto a_region = parse_region("6:29912399-29912528", human);
    auto a_region = parse_region("6:29913709-29914116", human);
    //auto a_region = parse_region("6:29912084-29912570", human);
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    ReadFilter<ReadIterator> a_read_filter {};
    a_read_filter.register_filter([] (const AlignedRead& the_read) {
        return is_good_mapping_quality(the_read, 10);
    });
    
    auto good_reads = filter_reads(std::move(reads), a_read_filter).first;
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 10));
    
    for (const auto& sample_reads : good_reads) {
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    unsigned ploidy {3};
    ReadModel a_read_model {ploidy};
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    unsigned max_haplotypes {64};
    HaplotypePhaser phaser {human, the_model, ploidy, max_haplotypes};
    
    Octopus::BayesianGenotypeModel::ReadRanges<ReadManager::SampleIdType,
    std::move_iterator<decltype(good_reads)::mapped_type::iterator>> read_ranges {};
    for (const auto& sample : samples) {
        read_ranges.emplace(sample, std::make_pair(std::make_move_iterator(good_reads[sample].begin()),
                                                   std::make_move_iterator(good_reads[sample].end())));
    }
    
//    phaser.put_data(read_ranges, candidates.cbegin(), candidates.cend());
//    
//    auto phased_regions = phaser.get_phased_regions(true);
//    
//    auto allele_posteriors = Octopus::VariantCaller::get_allele_posteriors(samples, phased_regions, candidates);
//    
//    for (const auto& haplotype_count : phased_regions.front().the_latent_posteriors.haplotype_pseudo_counts) {
//        if (haplotype_count.second > 1.0) {
//            cout << haplotype_count.first << endl;
//            haplotype_count.first.print_explicit_alleles();
//            cout << haplotype_count.second << endl;
//        }
//    }
//    
//    for (const auto& sample : samples) {
//        for (const auto& allele_posterior : allele_posteriors.at(sample)) {
//            cout << sample << " " << allele_posterior.first << " " << allele_posterior.second << endl;
//        }
//    }
}

TEST_CASE("reference allele posteriors in regions with no reads are the reference haplotype priors", "[variant_caller]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    auto samples = a_read_manager.get_sample_ids();
    
    auto a_region = parse_region("16:24740884-24741782", human);
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    ReadFilter<ReadIterator> a_read_filter {};
    a_read_filter.register_filter([] (const AlignedRead& the_read) {
        return is_good_mapping_quality(the_read, 10);
    });
    
    auto good_reads = filter_reads(std::move(reads), a_read_filter).first;
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 10));
    
    for (const auto& sample_reads : good_reads) {
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    unsigned ploidy {2};
    ReadModel a_read_model {ploidy};
    VariationalBayesGenotypeModel the_model {a_read_model, ploidy};
    
    unsigned max_haplotypes {64};
    HaplotypePhaser phaser {human, the_model, ploidy, max_haplotypes};
    
    Octopus::BayesianGenotypeModel::ReadRanges<ReadManager::SampleIdType,
    std::move_iterator<decltype(good_reads)::mapped_type::iterator>> read_ranges {};
    for (const auto& sample : samples) {
        read_ranges.emplace(sample, std::make_pair(std::make_move_iterator(good_reads[sample].begin()),
                                                   std::make_move_iterator(good_reads[sample].end())));
    }
    
//    phaser.put_data(read_ranges, candidates.cbegin(), candidates.cend());
//    
//    auto phased_regions = phaser.get_phased_regions(true);
//    
//    auto alleles_between_candidates = get_reference_alleles_between_variants(candidates, human);
//    
//    auto allele_posteriors = Octopus::VariantCaller::get_allele_posteriors(samples, phased_regions,
//                                                                           alleles_between_candidates);
//    
//    for (const auto& sample : samples) {
//        for (const auto& allele_posterior : allele_posteriors.at(sample)) {
//            cout << sample << " " << allele_posterior.first << " " << allele_posterior.second << endl;
//        }
//    }
}
