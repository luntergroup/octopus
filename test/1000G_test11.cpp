//
//  1000G_test11.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <cstddef>
#include <vector>
#include <unordered_map>

#include "test_common.hpp"
#include "test_utils.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "variant_utils.hpp"
#include "candidate_variant_generator.hpp"
#include "alignment_candidate_variant_generator.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "read_model.hpp"
#include "maths.hpp"
#include "read_filter.hpp"
#include "read_filters.hpp"
#include "haplotype_tree.hpp"

using std::cout;
using std::endl;

//BOOST_AUTO_TEST_CASE(1000G test 11 9:21725129-21725168)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human(a_factory.make(human_reference_fasta));
//    
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
//    
//    CandidateVariantGenerator candidate_generator {};
//    
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//    
//    auto sample_ids = a_read_manager.get_sample_ids();
//    auto the_sample_id = sample_ids.at(0);
//    
//    auto a_region = parse_region("9:21725129-21725168", human);
//    
//}
