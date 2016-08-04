// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef haplotype_filter_hpp
#define haplotype_filter_hpp

#include <vector>
#include <cstddef>
#include <functional>
#include <unordered_map>

#include <config/common.hpp>
#include <core/types/haplotype.hpp>

namespace octopus {

class HaplotypeLikelihoodCache;

std::vector<Haplotype>
filter_to_n(std::vector<Haplotype>& haplotypes, const std::vector<SampleName>& samples,
            const HaplotypeLikelihoodCache& haplotype_likelihoods, const std::size_t n);

using HaplotypeReference    = std::reference_wrapper<const Haplotype>;
using HaplotypePosteriorMap = std::unordered_map<HaplotypeReference, double>;

std::vector<HaplotypeReference>
extract_removable(const std::vector<Haplotype>& haplotypes,
                  const HaplotypePosteriorMap& haplotype_posteriors,
                  const std::vector<SampleName>& samples,
                  const HaplotypeLikelihoodCache& haplotype_likelihoods,
                  std::size_t max_to_remove, double min_posterior);

} // namespace octopus

#endif /* haplotype_filter_hpp */
