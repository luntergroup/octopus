//
//  search_region_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 27/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <iterator>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>

#include "test_common.hpp"
#include "mock_objects.hpp"
#include "reference_genome.hpp"
#include "mappable_algorithms.hpp"
#include "test_common.hpp"
#include "read_manager.hpp"
#include "read_filter.hpp"
#include "read_filters.hpp"
#include "read_utils.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "candidate_variant_generator.hpp"
#include "alignment_candidate_variant_generator.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "haplotype_tree.hpp"
#include "read_model.hpp"
#include "search_regions.hpp"

using std::cout;
using std::endl;

using Octopus::advance_region;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(advance_region_always_gives_a_region_more_advanced_than_the_previous_region)
{
    auto human = make_reference(human_reference_fasta);
    
    ReadManager a_read_manager {human_1000g_bam1, human_1000g_bam2, human_1000g_bam3};
    
    auto a_region = parse_region("16:62646800-62647030", human);
    
    auto samples = a_read_manager.get_samples();
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    for (auto& sample_reads : reads) {
        std::sort(sample_reads.second.begin(), sample_reads.second.end());
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    auto sub_region = parse_region("16:62646915-62646941", human);
    
    // This case was causing trouble
    auto next_region = advance_region(sub_region, reads, candidates, 4, 3);
    
    BOOST_CHECK(ends_before(sub_region, next_region));
}

BOOST_AUTO_TEST_CASE(search_regions_contains_all_variants_in_list_exactly_once_when_max_indicators_is_0)
{
    auto human = make_reference(human_reference_fasta);
    
    ReadManager a_read_manager {human_1000g_bam1};
    
    auto a_region = parse_region("16:8999900-9400000", human);
    
    auto samples = a_read_manager.get_samples();
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    ReadFilter<ReadIterator> a_read_filter {};
    a_read_filter.register_filter([] (const AlignedRead& the_read) {
        return is_good_mapping_quality(the_read, 5);
    });
    
    auto good_reads = Octopus::filter_reads(std::move(reads), a_read_filter).first;
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    for (auto& sample_reads : good_reads) {
        std::sort(sample_reads.second.begin(), sample_reads.second.end());
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    auto sub_region = parse_region("16:8999900-8999900", human);
    
    std::set<Variant> included {};
    
    unsigned max_variants {3};
    unsigned max_indicators {0};
    bool bound_respected {true};
    
    while (ends_before(sub_region, a_region)) {
        sub_region = advance_region(sub_region, good_reads, candidates, max_variants, max_indicators);
        auto overlapped = overlap_range(candidates.cbegin(), candidates.cend(), sub_region);
        std::for_each(overlapped.begin(), overlapped.end(), [&included, &bound_respected] (const auto& variant) {
            if (included.count(variant) > 0) {
                bound_respected = false;
            } else {
                included.emplace(variant);
            }
        });
        
        if (!bound_respected) break;
    }
    
    BOOST_CHECK(bound_respected);
    
    for (const auto& variant : candidates) {
        if (included.count(variant) == 0) {
            bound_respected = false;
            break;
        }
    }
    
    BOOST_CHECK(bound_respected);
    
    included.clear();
    
    sub_region = parse_region("16:8999900-8999900", human);
    max_variants = 5;
    
    while (ends_before(sub_region, a_region)) {
        sub_region = advance_region(sub_region, good_reads, candidates, max_variants, max_indicators);
        auto overlapped = overlap_range(candidates.cbegin(), candidates.cend(), sub_region);
        std::for_each(overlapped.begin(), overlapped.end(), [&included, &bound_respected] (const auto& variant) {
            if (included.count(variant) > 0) {
                bound_respected = false;
            } else {
                included.emplace(variant);
            }
        });
        
        if (!bound_respected) break;
    }
    
    BOOST_CHECK(bound_respected);
    
    for (const auto& variant : candidates) {
        if (included.count(variant) == 0) {
            bound_respected = false;
            break;
        }
    }
    
    BOOST_CHECK(bound_respected);
    
    included.clear();
    
    sub_region = parse_region("16:8999900-8999900", human);
    max_variants = 8;
    
    while (ends_before(sub_region, a_region)) {
        sub_region = advance_region(sub_region, good_reads, candidates, max_variants, max_indicators);
        auto overlapped = overlap_range(candidates.cbegin(), candidates.cend(), sub_region);
        std::for_each(overlapped.begin(), overlapped.end(), [&included, &bound_respected] (const auto& variant) {
            if (included.count(variant) > 0) {
                bound_respected = false;
            } else {
                included.emplace(variant);
            }
        });
        
        if (!bound_respected) break;
    }
    
    BOOST_CHECK(bound_respected);
    
    for (const auto& variant : candidates) {
        if (included.count(variant) == 0) {
            bound_respected = false;
            break;
        }
    }
    
    BOOST_CHECK(bound_respected);
}

BOOST_AUTO_TEST_CASE(advance_regions_bounds_are_respected)
{
    auto human = make_reference(human_reference_fasta);
    
    ReadManager a_read_manager {human_1000g_bam1};
    
    auto a_region = parse_region("16:8999900-9400000", human);
    
    auto samples = a_read_manager.get_samples();
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    ReadFilter<ReadIterator> a_read_filter {};
    a_read_filter.register_filter([] (const AlignedRead& the_read) {
        return is_good_mapping_quality(the_read, 5);
    });
    
    auto good_reads = Octopus::filter_reads(std::move(reads), a_read_filter).first;
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    for (auto& sample_reads : good_reads) {
        std::sort(sample_reads.second.begin(), sample_reads.second.end());
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    auto sub_region = parse_region("16:8999900-8999900", human);
    
    unsigned max_variants {3};
    unsigned max_indicators {0};
    bool bound_respected {true};
    
    while (ends_before(sub_region, a_region)) {
        sub_region = advance_region(sub_region, good_reads, candidates, max_variants, max_indicators);
        auto num_variants_in_sub_region = count_overlapped(candidates.cbegin(), candidates.cend(), sub_region);
        if (num_variants_in_sub_region > max_variants) {
            auto overlapped = overlap_range(candidates.cbegin(), candidates.cend(), sub_region);
            auto rightmost = rightmost_mappable(overlapped.begin(), overlapped.end());
            if (!overlaps(*rightmost, overlapped.back())) {
                bound_respected = false;
                break;
            }
        }
    }
    
    BOOST_CHECK(bound_respected);

    sub_region = parse_region("16:8999900-8999900", human);
    max_variants = 5;
    
    while (ends_before(sub_region, a_region)) {
        sub_region = advance_region(sub_region, good_reads, candidates, max_variants, max_indicators);
        auto num_variants_in_sub_region = count_overlapped(candidates.cbegin(), candidates.cend(), sub_region);
        if (num_variants_in_sub_region > max_variants) {
            auto overlapped = overlap_range(candidates.cbegin(), candidates.cend(), sub_region);
            auto rightmost = rightmost_mappable(overlapped.begin(), overlapped.end());
            if (!overlaps(*rightmost, overlapped.back())) {
                bound_respected = false;
                break;
            }
        }
    }
    
    BOOST_CHECK(bound_respected);

    sub_region = parse_region("16:8999900-8999900", human);
    max_variants = 8;
    
    while (ends_before(sub_region, a_region)) {
        sub_region = advance_region(sub_region, good_reads, candidates, max_variants, max_indicators);
        auto num_variants_in_sub_region = count_overlapped(candidates.cbegin(), candidates.cend(), sub_region);
        if (num_variants_in_sub_region > max_variants) {
            auto overlapped = overlap_range(candidates.cbegin(), candidates.cend(), sub_region);
            auto rightmost = rightmost_mappable(overlapped.begin(), overlapped.end());
            if (!overlaps(*rightmost, overlapped.back())) {
                bound_respected = false;
                break;
            }
        }
    }
    
    BOOST_CHECK(bound_respected);

    sub_region = parse_region("16:8999900-8999900", human);
    max_variants   = 3;
    max_indicators = 1;
    
    while (ends_before(sub_region, a_region)) {
        sub_region = advance_region(sub_region, good_reads, candidates, max_variants, max_indicators);
        auto num_variants_in_sub_region = count_overlapped(candidates.cbegin(), candidates.cend(), sub_region);
        if (num_variants_in_sub_region > max_variants) {
            auto overlapped = overlap_range(candidates.cbegin(), candidates.cend(), sub_region);
            auto rightmost = rightmost_mappable(overlapped.begin(), overlapped.end());
            if (!overlaps(*rightmost, overlapped.back())) {
                bound_respected = false;
                break;
            }
        }
    }
    
    BOOST_CHECK(bound_respected);

    sub_region = parse_region("16:8999900-8999900", human);
    max_variants   = 5;
    max_indicators = 2;
    
    while (ends_before(sub_region, a_region)) {
        sub_region = advance_region(sub_region, good_reads, candidates, max_variants, max_indicators);
        auto num_variants_in_sub_region = count_overlapped(candidates.cbegin(), candidates.cend(), sub_region);
        if (num_variants_in_sub_region > max_variants) {
            auto overlapped = overlap_range(candidates.cbegin(), candidates.cend(), sub_region);
            auto rightmost = rightmost_mappable(overlapped.begin(), overlapped.end());
            if (!overlaps(*rightmost, overlapped.back())) {
                bound_respected = false;
                break;
            }
        }
    }
    
    BOOST_CHECK(bound_respected);

    sub_region = parse_region("16:8999900-8999900", human);
    max_variants   = 8;
    max_indicators = 3;
    
    while (ends_before(sub_region, a_region)) {
        sub_region = advance_region(sub_region, good_reads, candidates, max_variants, max_indicators);
        auto num_variants_in_sub_region = count_overlapped(candidates.cbegin(), candidates.cend(), sub_region);
        if (num_variants_in_sub_region > max_variants) {
            auto overlapped = overlap_range(candidates.cbegin(), candidates.cend(), sub_region);
            auto rightmost = rightmost_mappable(overlapped.begin(), overlapped.end());
            if (!overlaps(*rightmost, overlapped.back())) {
                bound_respected = false;
                break;
            }
        }
    }
    
    BOOST_CHECK(bound_respected);
}

BOOST_AUTO_TEST_CASE(setting_max_included_to_zero_in_advance_region_results_in_the_largest_variant_free_reference_region)
{
    auto human = make_reference(human_reference_fasta);
    
    ReadManager a_read_manager {human_1000g_bam1};
    
    auto a_region = parse_region("16:8999900-9400000", human);
    
    auto samples = a_read_manager.get_samples();
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    ReadFilter<ReadIterator> a_read_filter {};
    a_read_filter.register_filter([] (const AlignedRead& the_read) {
        return is_good_mapping_quality(the_read, 5);
    });
    
    auto good_reads = Octopus::filter_reads(std::move(reads), a_read_filter).first;
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    for (auto& sample_reads : good_reads) {
        std::sort(sample_reads.second.begin(), sample_reads.second.end());
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    auto sub_region = parse_region("16:9198984-9199045", human);
    
    unsigned max_variants {0};
    unsigned max_indicators {0};
    
    sub_region = advance_region(sub_region, reads, candidates, max_variants, max_indicators);
    
    BOOST_CHECK(sub_region == parse_region("16:9199045-9199306", human));
    
    sub_region = parse_region("16:9006622-9006661", human);
    sub_region = advance_region(sub_region, reads, candidates, max_variants, max_indicators);
    
    BOOST_CHECK(sub_region == parse_region("16:9006661-9006671", human));
}

BOOST_AUTO_TEST_CASE(search_regions_includes_all_possible_indicators)
{
    auto human = make_reference(human_reference_fasta);
    
    ReadManager a_read_manager {human_1000g_bam1};
    
    auto a_region = parse_region("16:62646780-62647130", human);
    
    auto samples = a_read_manager.get_samples();
    
    auto reads = a_read_manager.fetch_reads(samples, a_region);
    
    CandidateVariantGenerator candidate_generator {};
    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
    
    for (auto& sample_reads : reads) {
        std::sort(sample_reads.second.begin(), sample_reads.second.end());
        candidate_generator.add_reads(sample_reads.second.cbegin(), sample_reads.second.cend());
    }
    
    auto candidates = candidate_generator.get_candidates(a_region);
    
    auto sub_region = parse_region("16:62646834-62646972", human);
    
    auto next_region = advance_region(sub_region, reads, candidates, 12, 7);
    
    BOOST_CHECK(next_region == parse_region("16:62646834-62646990", human));
}

//BOOST_AUTO_TEST_CASE(advance_region_works_on_unusual_and_messy_data)
//{
//    auto candidates = generate_random_regions(10000, 5, 10000);
//    
//    GenomicRegion region {"test", 0, 10000};
//    unsigned max_variants {10};
//    unsigned max_indicators {5};
//    
//    std::unordered_map<std::string, std::vector<GenomicRegion>> reads {};
//    reads.emplace("t", generate_random_regions(10000, 100, 10000));
//    
//    auto regions = Octopus::cover_region(region, reads, candidates, max_variants, max_indicators);
//    
//    for (const auto& r : regions) {
//        cout << r << endl;
//    }
//}

BOOST_AUTO_TEST_SUITE_END()
