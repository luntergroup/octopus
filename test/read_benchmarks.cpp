//
//  read_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <memory>
#include <boost/filesystem/path.hpp>

#include "test_common.h"
#include "benchmark_utils.h"
#include "htslib_read_facade.h"
#include "read_manager.h"
#include "aligned_read.h"
#include "read_utils.h"
#include "mappable_algorithms.h"

using std::cout;
using std::endl;

//TEST_CASE("read_benchmark", "[read,benchmark]")
//{
//    HtslibReadFacade a_reader {human_1000g_bam1};
//    
//    auto f_read = [&a_reader] () {
//        a_reader.fetch_reads(GenomicRegion("10", 1000000, 1010000));
//    };
//    
//    auto without_vptr = benchmark<std::chrono::microseconds>(f_read, 10).count();
//    
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
//    
//    auto sample_ids = a_read_manager.get_sample_ids();
//    
//    auto the_sample_id = sample_ids.at(0);
//    
//    auto f_factory = [&a_read_manager, &the_sample_id] () {
//        a_read_manager.fetch_reads(the_sample_id, GenomicRegion("10", 1000000, 1010000));
//    };
//    
//    auto with_vptr = benchmark<std::chrono::microseconds>(f_factory, 10).count();
//    
//    std::cout << "Without vptr: " << without_vptr << "us" << std::endl;
//    std::cout << "With vptr: " << with_vptr << "us" << std::endl;
//}

//// Best time:
//// 93 milliseconds (MacBook Pro)
//TEST_CASE("reader_construct_destory_benchmark", "[read,benchmark]")
//{
//    boost::filesystem::path the_path {human_1000g_bam};
//    unsigned num_reference_contigs {};
//    
//    auto f_construct = [&the_path, &num_reference_contigs] () {
//        ReadReader the_reader {the_path, std::make_unique<HtslibReadFacade>(the_path)};
//        num_reference_contigs = the_reader.get_num_reference_contigs();
//    };
//    
//    auto time = benchmark<std::chrono::milliseconds>(f_construct, 100).count();
//    
//    std::cout << "Read construct time: " << time << "ms" << std::endl;
//}

//TEST_CASE("aligned_read_copy_benchmark", "[read,benchmark]")
//{
//    ReadFactory a_read_factory(std::vector<std::string> {human_1000g_bam});
//    auto sample_ids = a_read_factory.get_sample_ids();
//    auto the_sample_id = sample_ids.at(0);
//    
//    GenomicRegion a_big_region {"1", 0, 10000000};
//    GenomicRegion another_big_region {"2", 0, 10000000};
//    
//    auto reads_in_first_region  = a_read_factory.fetch_reads(the_sample_id, a_big_region);
//    auto reads_in_second_region = a_read_factory.fetch_reads(the_sample_id, another_big_region);
//    
//    auto f_lots_of_copies = [&reads_in_first_region, &reads_in_second_region] () {
//        auto a_copy = reads_in_first_region;
//        a_copy.insert(a_copy.end(), reads_in_second_region.begin(), reads_in_second_region.end());
//    };
//    
//    auto copy_time = benchmark<std::chrono::microseconds>(f_lots_of_copies, 100).count();
//    
//    //std::cout << "Copy: " << copy_time << "us" << std::endl;
//}

//TEST_CASE("reader coverage calculation benchmarks", "[read,benchmark]")
//{
//    ReadManager a_read_manager {std::vector<std::string> {human_1000g_bam1}};
//    
//    GenomicRegion a_region {"1", 1000000, 2000000};
//    
//    auto sample = a_read_manager.get_sample_ids().front();
//    
//    auto reads = a_read_manager.fetch_reads(sample, a_region);
//    
//    auto f_min_coverage = [&reads, &a_region] () {
//        unsigned coverage = min_coverage(reads, a_region);
//    };
//    
//    auto min_coverage_time = benchmark<std::chrono::nanoseconds>(f_min_coverage, 100).count();
//    auto time_per_position = min_coverage_time / size(a_region);
//    auto time_per_read     = min_coverage_time / reads.size();
//    
//    cout << "total time: " << min_coverage_time << "ns" << endl;
//    cout << "time per position: " << time_per_position << "ns" <<  endl;
//    cout << "time per read: " << time_per_read << "ns" <<  endl;
//}

