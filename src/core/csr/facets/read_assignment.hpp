// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_assignment_hpp
#define read_assignment_hpp

#include <unordered_map>
#include <vector>
#include <iosfwd>

#include <core/types/haplotype.hpp>
#include <basics/aligned_read.hpp>
#include <core/types/genotype.hpp>
#include <core/types/allele.hpp>

namespace octopus {

class HaplotypeLikelihoodModel;

using SupportSet = std::vector<AlignedRead>;

using HaplotypeSupportMap = std::unordered_map<Haplotype, SupportSet>;

using AlleleSupportMap = std::unordered_map<Allele, SupportSet>;

HaplotypeSupportMap compute_haplotype_support(const Genotype<Haplotype>& genotype,
                                              const std::vector<AlignedRead>& reads);

HaplotypeSupportMap compute_haplotype_support(const Genotype<Haplotype>& genotype,
                                              const std::vector<AlignedRead>& reads,
                                              HaplotypeLikelihoodModel model);

AlleleSupportMap compute_allele_support(const std::vector<Allele>& alleles,
                                        const HaplotypeSupportMap& haplotype_support);

} // namespace octopus

#endif /* read_assignment_hpp */
