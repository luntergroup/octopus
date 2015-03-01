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
        ;
        
        po::options_description backend("Backend options");
        backend.add_options()
        ("t", po::value<unsigned>(), "the number of threads")
        ("c", po::value<bool>()->default_value(false), "compress reads")
        ("verbosity", po::value<unsigned>()->default_value(0), "Level of logging")
        ;
        
        po::options_description calling("Caller options");
        calling.add_options()
        ("ref-calls", po::value<bool>()->default_value(false), "Output ref-calls")
        ;
        
        po::options_description assembler("Assembler options");
        assembler.add_options()
        ("assemble", po::value<bool>()->default_value(true), "Assemble reads")
        ("k", po::value<unsigned>()->default_value(15), "K-mer length")
        ("no-cycles", po::value<bool>()->default_value(false), "Dissalow cycles in assembly graph")
        ;
        
        po::options_description all("Allowed options");
        all.add(general).add(backend).add(assembler);
        
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
