//
//  octopus.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "octopus.hpp"

#include <iostream>
#include <thread>
#include <future>
#include <memory>
#include <algorithm>

#include "common.hpp"
#include "program_options.hpp"
#include "mappable_map.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "read_filter.hpp"
#include "read_transform.hpp"
#include "read_utils.hpp"
#include "candidate_generators.hpp"
#include "vcf.hpp"
#include "variant_caller.hpp"
#include "basic_caller.hpp"
#include "cancer_caller.hpp"

using std::cout;
using std::endl;

namespace Octopus
{
    auto get_contigs(const SearchRegions& regions)
    {
        std::vector<GenomicRegion::StringType> result {};
        result.reserve(regions.size());
        std::transform(std::cbegin(regions), std::cend(regions), std::back_inserter(result),
                       [] (const auto& p) { return p.first; });
        return result;
    }
    
    void run_octopus(po::variables_map& options)
    {
        cout << "starting Octopus" << endl;
        
        //auto num_system_threads = std::thread::hardware_concurrency(); // just a hint
        //if (num_system_threads == 0) num_system_threads = 1;
        
        //auto max_threads  = Octopus::get_num_threads(options);
        
        auto reference           = Octopus::get_reference(options);
        auto read_manager        = Octopus::get_read_manager(options);
        auto regions             = Octopus::get_search_regions(options, reference);
        auto read_filter         = Octopus::get_read_filter(options);
        auto read_transform      = Octopus::get_read_transformer(options);
        auto candidate_generator = Octopus::get_candidate_generator(options, reference);
        auto vcf                 = Octopus::get_output_vcf(options);
        
        cout << "writing results to " << vcf.path().string() << endl;
        
        auto samples = read_manager.get_samples();
        
        auto contigs = get_contigs(regions);
        
        auto vcf_header_builder = get_default_header_builder().set_samples(samples);
        for (const auto& contig : contigs) {
            vcf_header_builder.add_contig(contig);
        }
        
        auto vcf_header = vcf_header_builder.build_once();
        
        vcf.write(vcf_header);
        
        for (const auto& contig_region : regions) {
            auto region = *contig_region.second.cbegin();
            
            cout << "looking at region " << region << endl;
            
            std::unique_ptr<VariantCaller> caller = std::make_unique<CancerVariantCaller>(reference, read_manager, read_filter, read_transform, candidate_generator);
            
            auto calls = caller->call_variants(region);
            
            for (auto& call : calls) {
                vcf.write(call);
            }
        }
    }
    
} // end namespace Octopus


