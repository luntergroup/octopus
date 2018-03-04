// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_assigner_hpp
#define read_assigner_hpp

#include <unordered_map>
#include <vector>
#include <deque>
#include <functional>
#include <iosfwd>

#include "core/types/haplotype.hpp"
#include "basics/aligned_read.hpp"
#include "basics/cigar_string.hpp"
#include "core/types/genotype.hpp"
#include "core/types/allele.hpp"

namespace octopus {

class HaplotypeLikelihoodModel;

using ReadSupportSet = std::vector<AlignedRead>;
using HaplotypeSupportMap = std::unordered_map<Haplotype, ReadSupportSet>;
using ReadRefSupportSet = std::vector<std::reference_wrapper<const AlignedRead>>;
using AlleleSupportMap = std::unordered_map<Allele, ReadRefSupportSet>;

HaplotypeSupportMap compute_haplotype_support(const Genotype<Haplotype>& genotype,
                                              const std::vector<AlignedRead>& reads);

HaplotypeSupportMap compute_haplotype_support(const Genotype<Haplotype>& genotype,
                                              const std::vector<AlignedRead>& reads,
                                              std::deque<AlignedRead>& unassigned);

HaplotypeSupportMap compute_haplotype_support(const Genotype<Haplotype>& genotype,
                                              const std::vector<AlignedRead>& reads,
                                              HaplotypeLikelihoodModel model);

HaplotypeSupportMap compute_haplotype_support(const Genotype<Haplotype>& genotype,
                                              const std::vector<AlignedRead>& reads,
                                              std::deque<AlignedRead>& unassigned,
                                              HaplotypeLikelihoodModel model);

AlleleSupportMap compute_allele_support(const std::vector<Allele>& alleles,
                                        const HaplotypeSupportMap& haplotype_support);

} // namespace octopus

#endif
