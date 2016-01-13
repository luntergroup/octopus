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
#include <cstddef>
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
#include "candidate_generator_builder.hpp"
#include "vcf.hpp"
#include "vcf_utils.hpp"
#include "variant_caller.hpp"

#include "test_common.hpp"
#include "mappable_debug.hpp"

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
        std::vector<GenomicRegion::ContigNameType> result {};
        result.reserve(regions.size());
        std::transform(std::cbegin(regions), std::cend(regions), std::back_inserter(result),
                       [] (const auto& p) { return p.first; });
        return result;
    }
    
    bool check_search_regions(const SearchRegions& regions, const ReferenceGenome& reference)
    {
        bool result {true};
        
        for (const auto& p : regions) {
            if (!reference.has_contig(p.first)) {
                std::cout << "Bad input: contig " << p.first <<
                    " does not exist in reference " << reference.get_name() << std::endl;
                result = false;
            }
        }
        
        return result;
    }
    
    std::vector<SampleIdType> get_samples(const po::variables_map& options,
                                          const ReadManager& read_manager)
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
                          const std::vector<GenomicRegion::ContigNameType>& contigs,
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
    
    GenomicRegion::SizeType search_regions_size(const SearchRegions& regions)
    {
        return std::accumulate(std::cbegin(regions), std::cend(regions), GenomicRegion::SizeType {},
                               [] (const auto curr, const auto& p) {
                                   return curr + sum_sizes(p.second);
                               });
    }
    
    void run_octopus(const po::variables_map& options)
    {
        using std::cout; using std::endl;
        
        const auto reference = Options::make_reference(options);
        
        if (!reference.is_good()) {
            cout << "Octopus: quiting as got bad reference genome" << endl;
            return;
        }
        
        const auto regions = Options::get_search_regions(options, reference);
        
        if (!check_search_regions(regions, reference)) {
            cout << "Octopus: quiting as got bad input regions" << endl;
            return;
        }
        
        auto read_manager = Options::make_read_manager(options);
        
        const auto samples = get_samples(options, read_manager);
        
        auto read_filter    = Options::make_read_filter(options);
        auto downsampler    = Options::make_downsampler(options);
        auto read_transform = Options::make_read_transform(options);
        
        auto output = Options::make_output_vcf_writer(options);
        
        if (!output.is_open()) {
            cout << "Octopus: quiting as could not make output file" << endl;
            return;
        }
        
        auto candidate_generator_builder = Options::make_candidate_generator_builder(options, reference);
        
        if (candidate_generator_builder.num_generators() == 0) {
            std::cout << "Octopus: quiting as there are no candidate generators" << std::endl;
            return;
        }
        
        ReadPipe read_pipe {read_manager, read_filter, downsampler, read_transform};
        
        cout << "writing results to " << output.path().string() << endl;
        cout << "there are " << samples.size() << " samples" << endl;
        
        const auto contigs = get_contigs(regions);
        
        auto vcf_header = make_header(samples, contigs, reference);
        
        output.write(vcf_header);
        
        const size_t max_reads = 1'000'000;
        
        for (const auto& contig_regions : regions) {
            const auto& contig = contig_regions.first;
            
            size_t num_buffered_reads {};
            
            auto caller = Options::make_variant_caller(options, reference, candidate_generator_builder, contig);
            
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
                        //cout << call << endl;
                        output.write(std::move(call));
                    }
                    
                    num_buffered_reads = caller->num_buffered_reads();
                    
                    subregion = get_right_overhang(region, subregion);
                    subregion = read_manager.find_covered_subregion(samples, subregion, max_reads);
                }
            }
        }
        
        cout << "processed " << search_regions_size(regions) << "bp" << std::endl;
    }
    
} // namespace Octopus


