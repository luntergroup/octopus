// Copyright (c) 2015-2020 Daniel Cooke
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

enum class ReferencePadPolicy { leave, trim_alt_alleles, trim_all_alleles };

std::vector<boost::optional<Allele>>
get_resolved_alleles(const std::vector<VcfRecord>& calls,
                     std::size_t call_idx,
                     const VcfRecord::SampleName& sample,
                     ReferencePadPolicy ref_pad_policy = ReferencePadPolicy::trim_all_alleles);

GenotypeMap
extract_genotypes(const std::vector<VcfRecord>& calls,
                  const std::vector<VcfRecord::SampleName>& samples,
                  const ReferenceGenome& reference,
                  boost::optional<GenomicRegion> call_region = boost::none);

}

#endif
