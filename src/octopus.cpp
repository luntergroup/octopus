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
#include "downsampler.hpp"
#include "read_transform.hpp"
#include "read_pipe.hpp"
#include "read_utils.hpp"
#include "candidate_generators.hpp"
#include "vcf.hpp"
#include "vcf_utils.hpp"
#include "variant_caller.hpp"

#include "test_common.hpp"

#include "haplotype_tree.hpp"

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
    
    VcfHeader make_header(const std::vector<SampleIdType>& samples,
                          const std::vector<GenomicRegion::StringType>& contigs,
                          const ReferenceGenome& reference)
    {
        auto vcf_header_builder = get_default_header_builder().set_samples(samples);
        for (const auto& contig : contigs) {
            vcf_header_builder.add_contig(contig, {{"length", std::to_string(reference.get_contig_size(contig))}});
        }
        vcf_header_builder.add_basic_field("reference", reference.get_name());
        vcf_header_builder.add_structured_field("Octopus", {{"some", "option"}});
        
        return vcf_header_builder.build_once();
    }
    
    void run_octopus(po::variables_map& options)
    {
        using std::cout; using std::endl;
        
        //auto num_system_threads = std::thread::hardware_concurrency(); // just a hint
        //if (num_system_threads == 0) num_system_threads = 1;
        
        //auto max_threads  = Octopus::get_num_threads(options);
        
        //auto memory_quota = Options::get_memory_quota(options);
        
        const size_t max_reads = 1'000'000;
        
        auto reference           = Options::get_reference(options);
        auto read_manager        = Options::get_read_manager(options);
        auto regions             = Options::get_search_regions(options, reference);
        auto read_filter         = Options::get_read_filter(options);
        auto downsampler         = Options::get_downsampler(options);
        auto read_transform      = Options::get_read_transformer(options);
        auto candidate_generator = Options::get_candidate_generator(options, reference);
        auto output              = Options::get_output_vcf(options);
        
        ReadPipe read_pipe {read_manager, read_filter, downsampler, read_transform};
        
        const auto samples = get_samples(options, read_manager);
        
        cout << "there are " << samples.size() << " samples" << endl;
        
        cout << "writing results to " << output.path().string() << endl;
        
        const auto contigs = get_contigs(regions);
        
        auto vcf_header = make_header(samples, contigs, reference);
        
        output.write(vcf_header);
        
        for (const auto& contig_regions : regions) {
            const auto& contig = contig_regions.first;
            
            size_t num_buffered_reads {};
            
            auto caller = Options::get_variant_caller(options, reference, candidate_generator, contig);
            
            for (const auto& region : contig_regions.second) {
                cout << "processing input region " << region << endl;
                
                auto subregion = read_manager.find_covered_subregion(samples, region, max_reads);
                
                while (get_begin(subregion) != get_end(region)) {
                    cout << "processing subregion " << subregion << endl;
                    
                    auto reads = read_pipe.fetch_reads(samples, subregion);
                    
                    //return; // uncomment for ReadPipe performance benchmarking
                    
                    auto calls = caller->call_variants(subregion, std::move(reads));
                    
                    cout << "writing " << calls.size() << " calls to VCF" << endl;
                    
                    for (auto&& call : calls) {
                        cout << call << endl;
                        output.write(std::move(call));
                    }
                    
                    num_buffered_reads = caller->num_buffered_reads();
                    
                    subregion = get_right_overhang(region, subregion);
                    subregion = read_manager.find_covered_subregion(samples, subregion, max_reads);
                }
            }
        }
    }
    
} // namespace Octopus


