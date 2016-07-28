//
//  read_assignment.hpp
//  Octopus
//
//  Created by Daniel Cooke on 15/07/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef read_assignment_hpp
#define read_assignment_hpp

#include <unordered_map>
#include <vector>
#include <iosfwd>

#include "haplotype.hpp"
#include "aligned_read.hpp"
#include "genotype.hpp"
#include "allele.hpp"

namespace octopus
{
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
