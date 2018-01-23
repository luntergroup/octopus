// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef genotype_reader_hpp
#define genotype_reader_hpp

#include <functional>
#include <unordered_map>
#include <vector>
#include <utility>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "core/types/allele.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "containers/mappable_flat_set.hpp"
#include "core/types/variant.hpp"
#include "io/variant/vcf_record.hpp"

namespace octopus {

class ReferenceGenome;
class GenomicRegion;

using GenotypeMap = std::unordered_map<SampleName, MappableFlatSet<Genotype<Haplotype>>>;

std::pair<std::vector<Allele>, bool>
get_called_alleles(const VcfRecord& call, const VcfRecord::SampleName& sample, const bool trim_padding = true);

GenotypeMap extract_genotypes(const std::vector<VcfRecord>& calls,
                              const std::vector<VcfRecord::SampleName>& samples,
                              const ReferenceGenome& reference,
                              boost::optional<GenomicRegion> call_region = boost::none);

}

#endif
