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

#include <iostream>
#include <cstdlib>

#include "program_options.h"
#include "octopus.h"

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
#include "candidate_generators.h"
#include "mappable_set.h"
#include "mappable_map.h"
#include "haplotype_tree.h"
#include "genotype_model.h"
#include "population_genotype_model.h"
#include "vcf.h"

#include "sequence_utils.h"

using std::cout;
using std::endl;

void test()
{
    auto reference = make_reference(human_reference_fasta);
    
    ReadManager read_manager {human_1000g_bam1, human_1000g_bam2, human_1000g_bam3};
    
    //auto region = parse_region("2:104142870-104142884", reference);
    auto region = parse_region("11:67503118-67503253", reference);
    
    CandidateVariantGenerator generator {};
    generator.register_generator(std::make_unique<AlignmentCandidateVariantGenerator>(reference, 10));
    
    auto reads = make_mappable_map(read_manager.fetch_reads(region));
    
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
    
    auto genotype_posteriors = genotype_model->evaluate(haplotypes, reads);
    
    auto sample = read_manager.get_sample_ids()[2];
    auto it = std::max_element(genotype_posteriors[sample].cbegin(), genotype_posteriors[sample].cend(), [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    
    it->first.at(0).print_explicit_alleles();
    cout << "\n";
    it->first.at(1).print_explicit_alleles();
    cout << "\n" << it->second << endl;
    
    //cout << genotype_posteriors[sample].cbegin()->second << endl;
}

int main(int argc, const char **argv)
{
    test();
    exit(0);
    
    try {
        auto options = Octopus::parse_options(argc, argv);
        
        if (options.second) {
            Octopus::run_octopus(options.first);
            std::cout << "finished running Octopus" << std::endl;
        } else {
            std::cout << "did not run Octopus" << std::endl;
        }
        
    } catch (...) {
        std::cerr << "Encountered unknown error. Quiting now" << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
