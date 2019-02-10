// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_assigner_hpp
#define read_assigner_hpp

#include <unordered_map>
#include <vector>
#include <deque>
#include <functional>
#include <iosfwd>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/mappable_reference_wrapper.hpp"
#include "basics/aligned_read.hpp"
#include "basics/cigar_string.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/allele.hpp"

namespace octopus {

class HaplotypeLikelihoodModel;

using HaplotypeProbabilityMap = std::unordered_map<Haplotype, double>;
using ReadSupportSet = std::vector<AlignedRead>;
using HaplotypeSupportMap = std::unordered_map<Haplotype, ReadSupportSet>;
using AlignedReadConstReference = MappableReferenceWrapper<const AlignedRead>;
using ReadRefSupportSet = std::vector<AlignedReadConstReference>;
using AlleleSupportMap = std::unordered_map<Allele, ReadRefSupportSet>;

struct AmbiguousRead : public Mappable<AmbiguousRead>
{
    AlignedRead read;
    boost::optional<std::vector<Haplotype>> haplotypes = boost::none;
    AmbiguousRead(AlignedRead read) : read {std::move(read)} {}
    const auto& mapped_region() const noexcept { return octopus::mapped_region(read); }
};

using AmbiguousReadList = std::deque<AmbiguousRead>;

struct AssignmentConfig
{
    enum class AmbiguousAction { drop, first, random, all } ambiguous_action = AmbiguousAction::drop;
    enum class AmbiguousRecord { read_only, haplotypes, haplotypes_if_three_or_more_options };
    AmbiguousRecord ambiguous_record = AmbiguousRecord::haplotypes_if_three_or_more_options;
};

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          const HaplotypeProbabilityMap& log_priors,
                          AssignmentConfig config = AssignmentConfig {});

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AssignmentConfig config = AssignmentConfig {});

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          const HaplotypeProbabilityMap& log_priors,
                          AssignmentConfig config = AssignmentConfig {});

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          AssignmentConfig config = AssignmentConfig {});

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config = AssignmentConfig {});

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          const HaplotypeProbabilityMap& log_priors,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config = AssignmentConfig {});

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config = AssignmentConfig {});

template <typename BinaryPredicate>
AlleleSupportMap
compute_allele_support(const std::vector<Allele>& alleles,
                       const HaplotypeSupportMap& haplotype_support,
                       BinaryPredicate inclusion_pred)
{
    AlleleSupportMap result {};
    result.reserve(alleles.size());
    for (const auto& allele : alleles) {
        ReadRefSupportSet allele_support {};
        for (const auto& p : haplotype_support) {
            if (inclusion_pred(p.first, allele)) {
                allele_support.insert(std::cend(allele_support), std::cbegin(p.second), std::cend(p.second));
            }
        }
        std::sort(std::begin(allele_support), std::end(allele_support));
        result.emplace(allele, std::move(allele_support));
    }
    return result;
}

AlleleSupportMap
compute_allele_support(const std::vector<Allele>& alleles,
                       const HaplotypeSupportMap& haplotype_support);

std::size_t
try_assign_ambiguous_reads_to_alleles(const std::vector<Allele>& alleles,
                                      const AmbiguousReadList& ambiguous_reads,
                                      AlleleSupportMap& allele_support);

AlleleSupportMap
compute_allele_support(const std::vector<Allele>& alleles,
                       const HaplotypeSupportMap& haplotype_support,
                       const AmbiguousReadList& ambiguous_reads);

} // namespace octopus

#endif
