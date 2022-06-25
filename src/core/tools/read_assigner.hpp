// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_assigner_hpp
#define read_assigner_hpp

#include <unordered_map>
#include <vector>
#include <deque>
#include <functional>
#include <memory>
#include <iosfwd>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/mappable_reference_wrapper.hpp"
#include "basics/aligned_read.hpp"
#include "basics/aligned_template.hpp"
#include "basics/cigar_string.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/allele.hpp"
#include "utils/thread_pool.hpp"

namespace octopus {

class HaplotypeLikelihoodModel;

using HaplotypeProbabilityMap = std::unordered_map<Haplotype, double>;
using ReadSupportSet = std::vector<AlignedRead>;
using HaplotypeSupportMap = std::unordered_map<Haplotype, ReadSupportSet>;
using TemplateSupportSet = std::vector<AlignedTemplate>;
using HaplotypeTemplateSupportMap = std::unordered_map<Haplotype, TemplateSupportSet>;
using AlignedReadConstReference = MappableReferenceWrapper<const AlignedRead>;
using ReadRefSupportSet = std::vector<AlignedReadConstReference>;
using AlleleSupportMap = std::unordered_map<Allele, ReadRefSupportSet>;

struct AmbiguousRead : public Mappable<AmbiguousRead>
{
    AlignedRead read;
    boost::optional<std::vector<std::shared_ptr<Haplotype>>> haplotypes = boost::none;
    AmbiguousRead(AlignedRead read) : read {std::move(read)} {}
    AmbiguousRead(AlignedRead read, std::vector<std::shared_ptr<Haplotype>> haplotypes)
    : read {std::move(read)}, haplotypes {std::move(haplotypes)} {}
    const auto& mapped_region() const noexcept { return octopus::mapped_region(read); }
};

using AmbiguousReadList = std::deque<AmbiguousRead>;

struct AmbiguousTemplate : public Mappable<AmbiguousTemplate>
{
    AlignedTemplate read_template;
    boost::optional<std::vector<std::shared_ptr<Haplotype>>> haplotypes = boost::none;
    AmbiguousTemplate(AlignedTemplate read_template) : read_template {std::move(read_template)} {}
    AmbiguousTemplate(AlignedTemplate read_template, std::vector<std::shared_ptr<Haplotype>> haplotypes)
    : read_template {std::move(read_template)}, haplotypes {std::move(haplotypes)} {}
    const auto& mapped_region() const noexcept { return octopus::mapped_region(read_template); }
};

using AmbiguousTemplateList = std::deque<AmbiguousTemplate>;

struct AssignmentConfig
{
    enum class AmbiguousAction { drop, first, random, all } ambiguous_action = AmbiguousAction::drop;
    enum class AmbiguousRecord { read_only, haplotypes, haplotypes_if_three_or_more_options } ambiguous_record = AmbiguousRecord::haplotypes;
};

using OptionalThreadPool = boost::optional<ThreadPool&>;

// AlignedRead

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          const HaplotypeProbabilityMap& log_priors,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          const HaplotypeProbabilityMap& log_priors,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          const HaplotypeProbabilityMap& log_priors,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

// AlignedTemplate

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          const HaplotypeProbabilityMap& log_priors,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          AmbiguousTemplateList& ambiguous,
                          const HaplotypeProbabilityMap& log_priors,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          AmbiguousTemplateList& ambiguous,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          AmbiguousTemplateList& ambiguous,
                          const HaplotypeProbabilityMap& log_priors,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          AmbiguousTemplateList& ambiguous,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config = AssignmentConfig {},
                          OptionalThreadPool workers = boost::none);

const static auto default_inclusion_pred = [] (const Haplotype& haplotype, const Allele& allele) { return haplotype.includes(allele); };

AlleleSupportMap
compute_allele_support(const std::vector<Allele>& alleles,
                       const HaplotypeSupportMap& haplotype_support,
                       std::function<bool(const Haplotype&, const Allele&)> inclusion_pred = default_inclusion_pred);

std::size_t
try_assign_ambiguous_reads_to_alleles(const std::vector<Allele>& alleles,
                                      const AmbiguousReadList& ambiguous_reads,
                                      AlleleSupportMap& allele_support,
                                      std::function<bool(const Haplotype&, const Allele&)> inclusion_pred = default_inclusion_pred);

AlleleSupportMap
compute_allele_support(const std::vector<Allele>& alleles,
                       const HaplotypeSupportMap& haplotype_support,
                       const AmbiguousReadList& ambiguous_reads,
                       std::function<bool(const Haplotype&, const Allele&)> inclusion_pred = default_inclusion_pred);

} // namespace octopus

#endif
