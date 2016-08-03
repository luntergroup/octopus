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

#include <boost/optional.hpp>

#include <config/common.hpp>
#include <core/types/haplotype.hpp>
#include <core/types/genotype.hpp>
#include <containers/mappable_flat_set.hpp>
#include <core/types/variant.hpp>
#include <io/variant/vcf_record.hpp>

namespace octopus {

class VcfHeader;
class ReferenceGenome;
class GenomicRegion;

using GenotypeMap = std::unordered_map<SampleName, MappableFlatSet<Genotype<Haplotype>>>;

GenotypeMap extract_genotypes(const std::vector<VcfRecord>& calls, const VcfHeader& header,
                              const ReferenceGenome& reference,
                              boost::optional<GenomicRegion> call_region = boost::none);

}

#endif /* genotype_reader_hpp */
