//
//  haplotype_prior_model_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

//#include <iostream>
//#include <string>
//
//#include "test_common.hpp"
//#include "reference_genome.hpp"
//#include "read_manager.hpp"
//#include "allele.hpp"
//#include "variant.hpp"
//#include "variant_utils.hpp"
//#include "candidate_variant_generator.hpp"
//#include "alignment_candidate_variant_generator.hpp"
//#include "haplotype.hpp"
//#include "genotype.hpp"
//#include "haplotype_prior_model.hpp"
//
//using std::cout;
//using std::endl;

//BOOST_AUTO_TEST_CASE(haplotype_priors_can_be_evaluated)
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
//    
//    auto a_region = parse_region("11:27282186-27282290", human);
//    
//    auto samples = a_read_manager.get_sample_ids();
//    
//    auto reads = a_read_manager.fetch_reads(samples[0], a_region);
//    
//    CandidateVariantGenerator candidate_generator {};
//    candidate_generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(human, 0));
//    
//    candidate_generator.add_reads(reads.begin(), reads.end());
//    
//    auto candidates = candidate_generator.get_candidates(a_region);
//    
//    Haplotype reference_haplotype {human, a_region};
//    
//    Haplotype true_haplotype {human, a_region};
//    true_haplotype.push_back(Allele {parse_region("11:27282267-27282268", human), "G"});
//    
//    Haplotype false_haplotype1 {human, a_region};
//    false_haplotype1.push_back(Allele {parse_region("11:27282258-27282259", human), "C"});
//    false_haplotype1.push_back(Allele {parse_region("11:27282267-27282268", human), "G"});
//    
//    Haplotype false_haplotype2 {human, a_region};
//    false_haplotype2.push_back(Allele {parse_region("11:27282199-27282200", human), "C"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282207-27282208", human), "C"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282218-27282219", human), "C"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282239-27282240", human), "T"});
//    false_haplotype2.push_back(Allele {parse_region("11:27282267-27282268", human), "G"});
//    
////    cout << prior_probability<double>(reference_haplotype, candidates) << endl;
////    cout << prior_probability<double>(true_haplotype, candidates) << endl;
////    cout << prior_probability<double>(false_haplotype1, candidates) << endl;
////    cout << prior_probability<double>(false_haplotype2, candidates) << endl;
//}
