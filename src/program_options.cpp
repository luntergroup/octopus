//
//  program_options.cpp
//  Octopus
//
//  Created by Daniel Cooke on 27/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "program_options.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <iterator>
#include <algorithm>  // std::transform, std::min
#include <functional> // std::function
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include "genomic_region.h"
#include "reference_genome.h"
#include "string_utils.h"

namespace fs = boost::filesystem;

std::pair<po::variables_map, bool> parse_options(int argc, const char** argv)
{
    try {
        po::positional_options_description p;
        p.add("command", -1);
        
        po::options_description general("General options");
        general.add_options()
        ("help,h", "produce help message")
        ("version", "output the version number")
        ("verbosity", po::value<unsigned>()->default_value(0), "level of logging. Verbosity 0 switches off logging")
        ;
        
        po::options_description backend("Backend options");
        backend.add_options()
        ("num-threads,t", po::value<unsigned>(), "the number of threads")
        ("compress-reads", po::value<bool>()->default_value(false), "compress the reads (slower)")
        ("max-open-files", po::value<unsigned>()->default_value(20), "the maximum number of files that can be open at one time")
        ;
        
        po::options_description input("Input/output options");
        input.add_options()
        ("reference-file,R", po::value<std::string>(), "the reference genome file")
        ("read-file,I", po::value<std::vector<std::string>>()->multitoken(), "space-seperated list of aligned read file paths, or a path to a text file containing read file paths")
        ("regions,R", po::value<std::vector<std::string>>()->multitoken(), "space-seperated list of one-indexed variant search regions (chrom:begin-end), or a path to a file containing such regions")
        ("skip-regions", po::value<std::vector<std::string>>()->multitoken(), "space-seperated list of one-indexed regions (chrom:begin-end) to skip, or a path to a file containing such regions")
        ("known-variants", po::value<std::string>(), "a variant file path containing known variants. These variants will automatically become candidates")
        ("output,o", po::value<std::string>(), "the path of the output variant file")
        ("log-file", po::value<std::string>(), "the path of the output log file")
        ;
        
        po::options_description filters("Read filter options");
        filters.add_options()
        ("min-map-qual", po::value<unsigned>()->default_value(20), "reads with smaller mapping quality are ignored")
        ("no-duplicates", po::value<bool>()->default_value(false), "removes duplicate reads")
        ("trim-soft-clipped", po::value<bool>()->default_value(false), "trims soft clipped parts of the read")
        ("trim-flanks", po::value<bool>()->default_value(false), "trims the flanks of all reads")
        ("trim-adapters", po::value<bool>()->default_value(true), "trims any overlapping regions that pass the fragment size")
        ;
        
        po::options_description candidates("Candidate generation options");
        candidates.add_options()
        ("candidates-from-alignments", po::value<bool>()->default_value(true), "generate candidate variants from the aligned reads")
        ("candidates-from-assembler", po::value<bool>()->default_value(true), "generate candidate variants with the assembler")
        ("min-base-quality", po::value<unsigned>()->default_value(15), "only base changes with quality above this value are considered for snp generation")
        ("k", po::value<unsigned>()->default_value(15), "k-mer size to use")
        ("no-cycles", po::value<bool>()->default_value(false), "dissalow cycles in assembly graph")
        ;
        
        po::options_description model("Model options");
        model.add_options()
        ("ploidy", po::value<unsigned>()->default_value(2), "the organism ploidy")
        ("snp-prior", po::value<double>()->default_value(0.003), "the prior probability of a snp")
        ("insertion-prior", po::value<double>()->default_value(0.003), "the prior probability of an insertion into the reference")
        ("deletion-prior", po::value<double>()->default_value(0.003), "the prior probability of a deletion from the reference")
        ;
        
        po::options_description calling("Caller options");
        calling.add_options()
        ("min-posterior", po::value<unsigned>()->default_value(15), "the minimum variant posterior probability")
        ;
        
        po::options_description all("Allowed options");
        all.add(general).add(backend).add(candidates);
        
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(all).positional(p).run(), vm);
        po::notify(vm);
        
        if (vm.count("help")) {
            std::cout << "Usage: octopus <command> [options]" << std::endl;
            std::cout << all << std::endl;
            return {vm, false};
        }
        
        return {vm, true};
    }
    catch(std::exception& e) {
        std::cout << e.what() << std::endl;
        return {po::variables_map {}, false};
    }
}

namespace detail
{
    bool is_region_file_path(const std::string& region_option)
    {
        return fs::native(region_option);
    }
    
    struct Line
    {
        std::string line_data;
        
        operator std::string() const
        {
            return line_data;
        }
    };
    
    std::istream& operator>>(std::istream& str, Line& data)
    {
        std::getline(str, data.line_data);
        return str;
    }
    
    std::string to_region_format(const std::string& bed_line)
    {
        auto tokens = split(bed_line, '\t');
        
        switch (tokens.size()) {
            case 0:
                throw std::runtime_error {"Empty line in input region bed file"};
            case 1:
                return std::string {tokens[0]};
            case 2:
                // Assume this represents a half range rather than a point
                return std::string {tokens[0] + ':' + tokens[1] + '-'};
            default:
                return std::string {tokens[0] + ':' + tokens[1] + '-' + tokens[2]};
        }
    }
    
    std::function<GenomicRegion(std::string)> get_line_parser(const fs::path& the_region_path,
                                                              const ReferenceGenome& the_reference)
    {
        if (the_region_path.extension().string() == ".bed") {
            return [&the_reference] (const std::string& line) {
                return parse_region(detail::to_region_format(line), the_reference);
            };
        } else {
            return [&the_reference] (const std::string& line) {
                return parse_region(line, the_reference);
            };
        }
    }
    
    std::vector<GenomicRegion> get_regions_from_file(const std::string& region_path, const ReferenceGenome& the_reference)
    {
        std::vector<GenomicRegion> result {};
        
        fs::path the_path {region_path};
        
        if (!fs::exists(the_path)) {
            throw std::runtime_error {"cannot find given region file " + the_path.string()};
        }
        
        std::ifstream the_file {the_path.string()};
        
        std::transform(std::istream_iterator<Line>(the_file), std::istream_iterator<Line>(),
                       std::back_inserter(result), get_line_parser(the_path, the_reference));
        
        return result;
    }
    
} // end namespace detail

std::vector<GenomicRegion> parse_region_option(const po::variables_map& options, const std::string& region_option,
                                               const ReferenceGenome& the_reference)
{
    if (options.count(region_option) == 0) {
        return get_all_contig_regions(the_reference);
    } else {
        std::vector<GenomicRegion> result {};
        
        const auto& given_regions = options.at(region_option).as<std::vector<std::string>>();
        
        for (const auto& region : given_regions) {
            if (detail::is_region_file_path(region)) {
                auto regions_from_file = detail::get_regions_from_file(region, the_reference);
                result.insert(result.end(), std::make_move_iterator(regions_from_file.begin()),
                              std::make_move_iterator(regions_from_file.end()));
            } else {
                result.emplace_back(parse_region(region, the_reference));
            }
        }
        
        return result;
    }
}
