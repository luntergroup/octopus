//
//  reference_genome_implementor_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "reference_genome_implementor_factory.h"

#include "fasta.h"

#include <iostream>

std::unique_ptr<IReferenceGenomeImplementor>
ReferenceGenomeImplementorFactory::make(std::string genome_file_path) const
{
    return std::make_unique<Fasta>(genome_file_path);
}
