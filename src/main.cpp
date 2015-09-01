//
//  main.cpp
//  Octopus
//
//  Created by Daniel Cooke on 03/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MODULE Main
//#include <boost/test/unit_test.hpp>

//#include <iostream>
//#include <cstdlib>
//
//#include "program_options.h"
//#include "octopus.h"

#include <iostream>
#include <string>
#include <chrono>
#include <algorithm>

#include "test_common.h"
#include "genomic_region.h"
#include "variant.h"
#include "reference_genome.h"
#include "read_manager.h"
#include "mappable_set.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "mappable_set.h"
#include "mappable_map.h"
#include "haplotype_tree.h"
#include "genotype_model.h"
#include "population_genotype_model.h"

#include "sequence_utils.h"

using std::cout;
using std::endl;

void test()
{
    auto reference = make_reference(human_reference_fasta);
    
    ReadManager read_manager {human_1000g_bam1};
    
    auto region = parse_region("2:104142870-104142884", reference);
    
    CandidateVariantGenerator generator {};
    generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(reference, 10));
    
    auto sample = read_manager.get_sample_ids().front();
    
    auto reads = read_manager.fetch_reads(region);
    
    cout << "there are " << reads[sample].size() << " reads" << endl;
    
    MappableMap<std::string, AlignedRead> read_map {{sample, MappableSet<AlignedRead> {reads[sample].begin(), reads[sample].end()}}};
    
    generator.add_reads(reads[sample].cbegin(), reads[sample].cend());
    
    auto candidates = generator.get_candidates(region);
    
    cout << "there are " << candidates.size() << " candidates" << endl;
    
    Octopus::HaplotypeTree tree {reference};
    
    for (const auto& v : candidates) {
        tree.extend(v.get_reference_allele());
        tree.extend(v.get_alternative_allele());
    }
    
    auto haplotypes = tree.get_haplotypes(region);
    
    cout << "there are " << haplotypes.size() << " haplotypes" << endl;
    
    auto genotype_model = std::make_unique<Octopus::PopulationGenotypeModel>(1, 2);
    
    auto genotype_posteriors = genotype_model->evaluate(haplotypes, read_map);
    
    auto it = std::max_element(genotype_posteriors[sample].cbegin(), genotype_posteriors[sample].cend(), [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    
    cout << it->first << " " << it->second << endl;
    
    //cout << genotype_posteriors[sample].cbegin()->second << endl;
}

int main(int argc, const char **argv)
{
    test();
    
//    auto options = Octopus::parse_options(argc, argv);
//    
//    if (options.second) {
//        Octopus::run_octopus(options.first);
//        std::cout << "finished running Octopus" << std::endl;
//    } else {
//        std::cout << "did not run Octopus" << std::endl;
//    }
    
    return EXIT_SUCCESS;
}
