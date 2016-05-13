////
////  cancer_caller_tests.cpp
////  Octopus
////
////  Created by Daniel Cooke on 11/02/2016.
////  Copyright © 2016 Oxford University. All rights reserved.
////
//
//#define BOOST_TEST_DYN_LINK
//
//#include <boost/test/unit_test.hpp>
//
//#include <iostream>
//#include <string>
//#include <cstddef>
//#include <set>
//#include <vector>
//
//#include "test_common.hpp"
//#include "reference_genome.hpp"
//#include "read_manager.hpp"
//#include "variant.hpp"
//#include "candidate_variant_generator.hpp"
//#include "alignment_candidate_variant_generator.hpp"
//#include "haplotype.hpp"
//#include "haplotype_tree.hpp"
//#include "genotype.hpp"
//
//#include "cancer_genotype.hpp"
//
//using std::cout;
//using std::endl;
//
//static Haplotype make_haplotype(Allele allele, const ReferenceGenome& reference)
//{
//    Haplotype result {allele.mapped_region(), reference};
//    result.push_back(std::move(allele));
//    return result;
//}
//
//BOOST_AUTO_TEST_SUITE(Components)
//BOOST_AUTO_TEST_SUITE(CancerCalling)
//
//BOOST_AUTO_TEST_CASE(cancer_genotype_can_be_indexed_with_the_last_element_being_the_cancer_element)
//{
//    BOOST_REQUIRE(test_file_exists(human_reference_fasta));
//    
//    const auto human = make_reference(human_reference_fasta);
//    
//    const auto region = parse_region("1:1000000-1000001", human);
//    
//    Allele allele1 {region, "A"}, allele2 {region, "T"};
//    
//    Genotype<Allele> germline_genotype {allele1, allele2};
//    
//    CancerGenotype<Allele> cancer_genotype1 {germline_genotype, allele1};
//    
//    const auto hap1 = make_haplotype(allele1, human);
//    const auto hap2 = make_haplotype(allele2, human);
//    
//    CancerGenotype<Haplotype> cancer_genotype2 {};
//}
//
//BOOST_AUTO_TEST_CASE(can_generate_all_cancer_genotypes_from_a_set_of_haplotypes)
//{
//    
//}
//
//BOOST_AUTO_TEST_SUITE_END() // CancerCalling
//BOOST_AUTO_TEST_SUITE_END() // Components
