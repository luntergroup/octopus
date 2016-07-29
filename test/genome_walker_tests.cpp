//
//  genome_walker_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 27/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

//#include <iostream>
//#include <iterator>
//#include <vector>
//
//#include "test_common.hpp"
//#include "reference_genome.hpp"
//#include "read_manager.hpp"
//#include "candidate_variant_generators.hpp"
//#include "mappable_flat_multi_set.hpp"
//#include "genome_walker.hpp"
//
//using std::cout;
//using std::endl;
//
//using octopus::Composer::Builder;
//using octopus::GenomeWalker;

BOOST_AUTO_TEST_SUITE(Components)

//static MappableFlatMultiSet<Variant> make_mappable_set(const std::vector<Variant>& variants)
//{
//    return MappableFlatMultiSet<Variant> {std::cbegin(variants), std::cend(variants)};
//}
//
//BOOST_AUTO_TEST_CASE(walkers_return_empty_regions_when_done)
//{
//    // TODO
//}
//
//BOOST_AUTO_TEST_CASE(advance_region_always_gives_a_region_more_advanced_than_the_previous_region)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    constexpr unsigned max_variant_size {5};
//    
//    auto candidate_generator = Composer::Builder().set_reference(human)
//    .add_generator(Composer::Builder::Generator::Alignment).set_max_variant_size(max_variant_size).build();
//    
//    ReadManager read_manager {NA12878_low_coverage};
//    
//    const auto region = parse_region("6:29,723,537-29,725,747", human);
//    
//    const auto reads = read_manager.fetch_reads(region);
//    
//    add_reads(reads, candidate_generator);
//    
//    const auto candidates = make_mappable_set(candidate_generator.generate_candidates(region));
//    
//    GenomeWalker walker1 {5};
//    
//    //auto next = walker1.walk(region, reads, candidates);
//    
//    // TODO
//}
//
//BOOST_AUTO_TEST_CASE(a_walker_with_no_indicators_never_give_regions_with_shared_candidates)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    constexpr unsigned max_variant_size {5};
//    
//    auto candidate_generator = Composer::Builder().set_reference(human)
//    .add_generator(Composer::Builder::Generator::Alignment).set_max_variant_size(max_variant_size).build();
//    
//    ReadManager read_manager {NA12878_low_coverage};
//    
//    const auto region = parse_region("6:29,723,537-29,725,747", human);
//    
//    const auto reads = read_manager.fetch_reads(region);
//    
//    add_reads(reads, candidate_generator);
//    
//    const auto candidates = make_mappable_set(candidate_generator.generate_candidates(region));
//    
//    GenomeWalker walker1 {1};
//    
//    //auto cur_region = walker1.walk(region, reads, candidates);
//    
//    // TODO
//}

BOOST_AUTO_TEST_SUITE_END()
