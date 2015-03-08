//
//  program_options.cpp
//  Octopus
//
//  Created by Daniel Cooke on 27/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "program_options.h"

std::pair<po::variables_map, bool> parse_options(int argc, char** argv)
{
    try {
        po::positional_options_description p;
        p.add("command", -1);
        
        po::options_description general("General options");
        general.add_options()
        ("help", "produce help message")
        ("version", "output the version number")
        ("verbosity", po::value<unsigned>()->default_value(0), "level of logging. Verbosity 0 switches off logging")
        ;
        
        po::options_description backend("Backend options");
        backend.add_options()
        ("num-threads", po::value<unsigned>(), "the number of threads")
        ("compress-reads", po::value<bool>()->default_value(false), "compress the reads (slower)")
        ("max-open-files", po::value<unsigned>()->default_value(20), "the maximum number of files that can be open at one time")
        ;
        
        po::options_description input("Input options");
        input.add_options()
        ("ref-file", po::value<std::string>(), "the reference genome")
        ("read-file", po::value<std::string>(), "comma-delimited list of aligned read file paths, or a path to a text file containing read file paths")
        ("regions", po::value<std::string>(), "comma-delimited list of one-indexed regions (chrom:begin-end) to consider, or a path to a file containing such regions")
        ("skip-regions", po::value<std::string>(), "comma-delimited list of one-indexed regions (chrom:begin-end) to skip, or a path to a file containing such regions")
        ("known-variants", po::value<std::string>(), "a variant file path containing known variants. These variants will automatically become candidates")
        ("output", po::value<std::string>(), "the path of the output variant file")
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
        ("from-reads", po::value<bool>()->default_value(true), "generate candidate variants from the aligned reads")
        ("assemble", po::value<bool>()->default_value(true), "generate candidate variants with the assembler")
        ("k", po::value<unsigned>()->default_value(15), "k-mer size to use")
        ("no-cycles", po::value<bool>()->default_value(false), "dissalow cycles in assembly graph")
        ;
        
        po::options_description model("Model options");
        model.add_options()
        ("ploidy", po::value<unsigned>()->default_value(2), "the organism ploidy")
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
