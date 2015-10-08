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
#include <algorithm> // std::transform, std::find
#include <cstddef>   // size_t
#include <stdexcept>

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
#include "population_caller.hpp"
#include "cancer_caller.hpp"

#include "test_common.hpp"

#include "genotype.hpp"
#include "cancer_genotype.hpp"

namespace Octopus
{
    size_t count_reads(ReadManager& read_manager, ReferenceGenome& reference)
    {
        size_t result {};
        for (const auto& contig : reference.get_contig_names()) {
            result += read_manager.count_reads(reference.get_contig_region(contig));
        }
        return result;
    }
    
    std::vector<GenomicRegion> split_region_by_coverage(const GenomicRegion& region, size_t target_region_coverage,
                                                        ReadManager& read_manager)
    {
        std::vector<GenomicRegion> result {};
        
        
        
        return result;
    }
    
    auto get_contigs(const SearchRegions& regions)
    {
        std::vector<GenomicRegion::StringType> result {};
        result.reserve(regions.size());
        std::transform(std::cbegin(regions), std::cend(regions), std::back_inserter(result),
                       [] (const auto& p) { return p.first; });
        return result;
    }
    
    std::vector<SampleIdType> get_samples(const po::variables_map& options, const ReadManager& read_manager)
    {
        auto user_samples = Options::get_samples(options);
        auto file_samples = read_manager.get_samples();
        
        if (!user_samples.empty()) {
            std::vector<SampleIdType> bad_samples {};
            
            for (const auto& user_sample : user_samples) {
                if (std::find(std::cbegin(file_samples), std::cend(file_samples), user_sample) == std::cend(file_samples)) {
                    bad_samples.push_back(user_sample);
                }
            }
            
            if (!bad_samples.empty()) {
                std::string error {"input samples not in read files: "};
                for (const auto& sample : bad_samples) error += sample + ',';
                error.erase(--error.end()); // removes last ','
                throw std::runtime_error {error};
            }
            
            return user_samples;
        } else {
            return file_samples;
        }
    }
    
    size_t approx_num_reads(size_t bytes_available)
    {
        return bytes_available / sizeof(AlignedRead);
    }
    
    void run_octopus(po::variables_map& options)
    {
        using std::cout; using std::endl;
        
        //auto num_system_threads = std::thread::hardware_concurrency(); // just a hint
        //if (num_system_threads == 0) num_system_threads = 1;
        
        //auto max_threads  = Octopus::get_num_threads(options);
        
        //auto memory_quota = Options::get_memory_quota(options);
        
        auto reference           = Options::get_reference(options);
        auto read_manager        = Options::get_read_manager(options);
        auto regions             = Options::get_search_regions(options, reference);
        auto read_filter         = Options::get_read_filter(options);
        auto read_transform      = Options::get_read_transformer(options);
        auto candidate_generator = Options::get_candidate_generator(options, reference);
        auto caller              = Options::get_variant_caller(options, reference, candidate_generator);
        auto vcf                 = Options::get_output_vcf(options);
        
        //std::cout << "there are " << read_filter.num_filters() << " read filters" << std::endl;
        //std::cout << "there are " << read_transform.num_transforms() << " read transforms" << std::endl;
        
        const auto samples = get_samples(options, read_manager);
        
        cout << "model details: " << caller->get_details() << endl;
        cout << "writing results to " << vcf.path().string() << endl;
        
        auto contigs = get_contigs(regions);
        
        auto vcf_header_builder = get_default_header_builder().set_samples(samples);
        for (const auto& contig : contigs) vcf_header_builder.add_contig(contig);
        vcf_header_builder.add_basic_field("reference", reference.get_name());
        vcf_header_builder.add_structured_field("Octopus", {{"some", "option"}});
        
        const auto vcf_header = vcf_header_builder.build_once();
        
        vcf.write(vcf_header);
        
        for (const auto& contig_region : regions) {
            auto region = *contig_region.second.cbegin();
            
            cout << "processing region " << region << endl;
            
            auto good_reads = filter_reads(make_mappable_map(read_manager.fetch_reads(region)), read_filter).first;
            
            transform_reads(good_reads, read_transform);
            
            auto calls = caller->call_variants(region, std::move(good_reads));
            
            cout << "writing " << calls.size() << " calls to VCF" << endl;
            
            for (auto& call : calls) {
                cout << call << endl;
                vcf.write(call);
            }
        }
    }
    
} // namespace Octopus


