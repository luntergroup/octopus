//
//  genotype_reader.hpp
//  Octopus
//
//  Created by Daniel Cooke on 23/05/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef genotype_reader_hpp
#define genotype_reader_hpp

#include <functional>
#include <unordered_map>
#include <vector>

#include "common.hpp"
#include "reference_genome.hpp"
#include "vcf_reader.hpp"
#include "genomic_region.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "mappable_flat_set.hpp"

namespace Octopus
{
    class GenotypeReader
    {
    public:
        using GenotypeMap = std::unordered_map<SampleIdType, MappableFlatSet<Genotype<Haplotype>>>;
        
        explicit GenotypeReader(const ReferenceGenome& reference, VcfReader&& variant_reader);
        
        GenotypeMap extract_genotype(const GenomicRegion& region);
        
    private:
        std::reference_wrapper<const ReferenceGenome> reference_;
        
        VcfReader variant_reader_;
        std::vector<SampleIdType> samples_;
    };
} // namespace Octopus

#endif /* genotype_reader_hpp */
