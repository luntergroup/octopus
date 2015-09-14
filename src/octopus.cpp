//
//  octopus.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "octopus.h"

#include <iostream>
#include <thread>
#include <future>
#include <algorithm>

#include "common.h"
#include "program_options.h"
#include "mappable_map.h"
#include "reference_genome.h"
#include "read_manager.h"
#include "read_utils.h"
#include "candidate_generators.h"
#include "vcf.h"

#include "haplotype_tree.h"
#include "genotype_model.h"
#include "population_genotype_model.h"

using std::cout;
using std::endl;

namespace Octopus
{
    void run_octopus(po::variables_map& options)
    {
        cout << "starting Octopus" << endl;
        
        //auto num_system_threads = std::thread::hardware_concurrency(); // just a hint
        //if (num_system_threads == 0) num_system_threads = 1;
        
        //auto max_threads  = Octopus::get_num_threads(options);
        
        auto reference    = Octopus::get_reference(options);
        auto read_manager = Octopus::get_read_manager(options);
        auto regions      = Octopus::get_search_regions(options, reference);
        auto vcf          = Octopus::get_output_vcf(options);
        auto samples = read_manager.get_sample_ids();
        
        auto read_filter = Octopus::get_read_filter(options);
        
        auto candidate_generator = Octopus::get_candidate_generator(options, reference);
        
        auto region = *regions.cbegin()->second.cbegin();
        
        cout << "looking at region " << region << endl;
        
        auto reads = make_mappable_map(read_manager.fetch_reads(samples, region));
        
        auto filtered_reads = filter_reads(reads, read_filter);
        auto& good_reads    = filtered_reads.first;
        auto& bad_reads     = filtered_reads.second;
        
        cout << "found " << good_reads[samples.front()].size() << " good reads" << endl;
        cout << "found " << bad_reads[samples.front()].size() << " bad reads" << endl;
        
        add_reads(good_reads, candidate_generator);
        
        auto candidates = candidate_generator.get_candidates(region);
        
        cout << "found " << candidates.size() << " candidates" << endl;
        
        Octopus::HaplotypeTree tree {reference};
        extend_tree(candidates, tree);
        auto haplotypes = tree.get_haplotypes(region);
        
        cout << "there are " << haplotypes.size() << " haplotypes" << endl;
        
        auto genotype_model = std::make_unique<Octopus::PopulationGenotypeModel>(1, 2);
        
        auto genotype_posteriors = genotype_model->evaluate(haplotypes, reads);
        
        auto sample = samples.front();
        auto it = std::max_element(genotype_posteriors[sample].cbegin(), genotype_posteriors[sample].cend(), [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
        
        it->first.at(0).print_explicit_alleles();
        cout << "\n";
        it->first.at(1).print_explicit_alleles();
        cout << "\n" << it->second << endl;
    }
    
} // end namespace Octopus


