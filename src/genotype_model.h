//
//  genotype_model.h
//  Octopus
//
//  Created by Daniel Cooke on 20/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_genotype_model_h
#define Octopus_genotype_model_h

#include <string>
#include <vector>
#include <unordered_map>

#include "mappable_map.h"
#include "genotype.h"

class AlignedRead;
class Haplotype;

namespace Octopus
{
    class GenotypeModel
    {
    public:
        using SampleIdType                = std::string;
        using SampleGenotypeProbabilities = std::unordered_map<Genotype<Haplotype>, double>;
        using GenotypeProbabilities       = std::unordered_map<SampleIdType, SampleGenotypeProbabilities>;
        using ReadMap                     = MappableMap<SampleIdType, AlignedRead>;
        
        GenotypeProbabilities evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads);
        
        virtual ~GenotypeModel() = default;
        
    private:
        virtual GenotypeProbabilities do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads);
    };
    
    inline GenotypeModel::GenotypeProbabilities
    GenotypeModel::evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads)
    {
        return do_evaluate(haplotypes, reads);
    }
    
} // end namespace Octopus

#endif
