// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_assigner.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cassert>

#include <boost/optional.hpp>

#include "core/models/haplotype_likelihood_model.hpp"
#include "core/models/error/error_model_factory.hpp"
#include "utils/maths.hpp"
#include "utils/kmer_mapper.hpp"
#include "utils/random_select.hpp"
#include "utils/parallel_transform.hpp"

namespace octopus {

namespace {

using HaplotypeLikelihoods = std::vector<std::vector<double>>;

auto vectorise(const Genotype<Haplotype>& genotype, const HaplotypeProbabilityMap& priors)
{
    std::vector<double> result(genotype.ploidy());
    std::transform(std::cbegin(genotype), std::cend(genotype), std::begin(result),
                   [&] (const auto& haplotype) { return priors.at(haplotype); });
    return result;
}

auto get_priors(const Genotype<Haplotype>& genotype, const HaplotypeProbabilityMap& log_priors)
{
    if (log_priors.empty()) {
        return std::vector<double>(genotype.ploidy());
    } else {
        return vectorise(genotype, log_priors);
    }
}

void find_map_haplotypes(const Genotype<Haplotype>& genotype, const unsigned read,
                         const HaplotypeLikelihoods& likelihoods, const std::vector<double>& log_priors,
                         std::vector<unsigned>& result)
{
    assert(result.empty());
    auto max_likelihood = std::numeric_limits<double>::lowest();
    for (unsigned k {0}; k < genotype.ploidy(); ++k) {
        const auto curr = likelihoods[k][read] + log_priors[k];
        if (maths::almost_equal(curr, max_likelihood)) {
            result.push_back(k);
        } else if (curr > max_likelihood) {
            result.assign({k});
            max_likelihood = curr;
        }
    }
    if (result.empty()) {
        result.resize(genotype.ploidy());
        std::iota(std::begin(result), std::end(result), 0);
    }
}

template <typename Map, typename Aligned, typename Ambiguous>
void calculate_support(Map& result,
                       const Genotype<Haplotype>& genotype,
                       const std::vector<Aligned>& reads,
                       const std::vector<double>& log_priors,
                       const HaplotypeLikelihoods& likelihoods,
                       boost::optional<Ambiguous&> ambiguous,
                       const AssignmentConfig& config)
{
    std::vector<unsigned> top {};
    top.reserve(genotype.ploidy());
    std::vector<std::shared_ptr<Haplotype>> haplotype_ptrs {};
    if (config.ambiguous_record != AssignmentConfig::AmbiguousRecord::read_only) {
        haplotype_ptrs.resize(genotype.ploidy());
    }
    for (unsigned i {0}; i < reads.size(); ++i) {
        const auto& read = reads[i];
        find_map_haplotypes(genotype, i, likelihoods, log_priors, top);
        if (top.size() == 1) {
            result[genotype[top.front()]].push_back(read);
        } else {
            using UA = AssignmentConfig::AmbiguousAction;
            switch (config.ambiguous_action) {
                case UA::first:
                    result[genotype[top.front()]].push_back(read);
                    break;
                case UA::all: {
                    for (auto idx : top) result[genotype[idx]].push_back(read);
                    break;
                }
                case UA::random: {
                    result[genotype[random_select(top)]].push_back(read);
                    break;
                }
                case UA::drop:
                default:
                    break;
            }
            if (ambiguous) {
                ambiguous->emplace_back(read);
                if (config.ambiguous_record == AssignmentConfig::AmbiguousRecord::haplotypes
                    || (config.ambiguous_record == AssignmentConfig::AmbiguousRecord::haplotypes_if_three_or_more_options && top.size() >= 3)) {
                    ambiguous->back().haplotypes.emplace();
                    ambiguous->back().haplotypes->reserve(top.size());
                    for (auto idx : top) {
                        if (!haplotype_ptrs[idx]) haplotype_ptrs[idx] = std::make_shared<Haplotype>(genotype[idx]);
                        ambiguous->back().haplotypes->push_back(haplotype_ptrs[idx]);
                    }
                }
            }
        }
        top.clear();
    }
}

auto calculate_support(const Genotype<Haplotype>& genotypes,
                       const std::vector<AlignedRead>& reads,
                       const std::vector<double>& log_priors,
                       const HaplotypeLikelihoods& likelihoods,
                       boost::optional<AmbiguousReadList&> ambiguous,
                       const AssignmentConfig& config)
{
    HaplotypeSupportMap result {};
    calculate_support(result, genotypes, reads, log_priors, likelihoods, ambiguous, config);
    return result;
}

auto calculate_support(const Genotype<Haplotype>& genotype,
                       const std::vector<AlignedTemplate>& reads,
                       const std::vector<double>& log_priors,
                       const HaplotypeLikelihoods& likelihoods,
                       boost::optional<AmbiguousTemplateList&> ambiguous,
                       const AssignmentConfig& config)
{
    HaplotypeTemplateSupportMap result {};
    calculate_support(result, genotype, reads, log_priors, likelihoods, ambiguous, config);
    return result;
}

template <typename MappableTp>
GenomicRegion::Size estimate_max_indel_size_helper(const MappableTp& mappable)
{
    const auto p = std::minmax({region_size(mappable), static_cast<GenomicRegion::Size>(sequence_size(mappable))});
    return p.second - p.first;
}

GenomicRegion::Size estimate_max_indel_size_helper(const AlignedTemplate& reads)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), GenomicRegion::Size {0},
                           [] (auto curr, const auto& read) { return curr + estimate_max_indel_size_helper(read); });
}

template <typename Range>
auto estimate_max_indel_size(const Range& mappables)
{
    GenomicRegion::Size result {0};
    for (const auto& mappable : mappables) {
        result = std::max(result, estimate_max_indel_size_helper(mappable));
    }
    return result;
}

auto compute_read_hashes(const std::vector<AlignedRead>& reads, OptionalThreadPool workers = boost::none)
{
    static constexpr unsigned char mapperKmerSize {6};
    std::vector<KmerPerfectHashes> result(reads.size());
    using octopus::transform;
    transform(std::cbegin(reads), std::cend(reads), std::begin(result),
              [=] (const AlignedRead& read) { return compute_kmer_hashes<mapperKmerSize>(read.sequence()); },
              workers);
    return result;
}

auto compute_read_hashes(const std::vector<AlignedTemplate>& templates, OptionalThreadPool workers = boost::none)
{
    static constexpr unsigned char mapperKmerSize {6};
    std::vector<std::vector<KmerPerfectHashes>> result(templates.size());
    using octopus::transform;
    transform(std::cbegin(templates), std::cend(templates), std::begin(result),
              [&] (const auto& reads) { 
                  std::vector<KmerPerfectHashes> hashes {};
                  hashes.reserve(reads.size());
                  std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(hashes),
                                [=] (const AlignedRead& read) { return compute_kmer_hashes<mapperKmerSize>(read.sequence()); });
                  return hashes;
              },
              workers);
    return result;
}

auto expand_for_alignment(const Haplotype& haplotype, const GenomicRegion& reads_region,
                          const GenomicRegion::Size indel_factor, const HaplotypeLikelihoodModel& model)
{
    const auto min_flank_pad = 2 * model.pad_requirement();
    const auto& haplotype_region = mapped_region(haplotype);
    unsigned min_lhs_expansion {min_flank_pad}, min_rhs_expansion {min_flank_pad};
    if (begins_before(reads_region, haplotype_region)) {
        min_lhs_expansion += begin_distance(reads_region, haplotype_region);
    }
    if (ends_before(haplotype_region, reads_region)) {
        min_rhs_expansion += end_distance(haplotype_region, reads_region);
    }
    const auto min_expansion = std::max(min_lhs_expansion, min_rhs_expansion) + indel_factor;
    return expand(haplotype, min_expansion);
}

auto map_query_to_target_helper(const KmerPerfectHashes& query, const KmerHashTable& target,
                                MappedIndexCounts& mapping_counts)
{
    return map_query_to_target(query, target, mapping_counts);
}

std::vector<std::vector<std::size_t>>
map_query_to_target_helper(const std::vector<KmerPerfectHashes>& queries, const KmerHashTable& target,
                           MappedIndexCounts& mapping_counts)
{
    std::vector<std::vector<std::size_t>> result {};
    result.reserve(queries.size());
    for (const auto& query : queries) {
        result.push_back(map_query_to_target(query, target, mapping_counts));
        reset_mapping_counts(mapping_counts);
    }
    return result;
}

template <typename Container, typename ReadHashKmerHashTable>
auto calculate_likelihoods(const Haplotype& haplotype,
                           const Container& reads,
                           const KmerHashTable& haplotype_hashes,
                           const ReadHashKmerHashTable& read_hashes,
                           MappedIndexCounts& haplotype_mapping_counts,
                           const HaplotypeLikelihoodModel& model)
{
    std::vector<double> result(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::cbegin(read_hashes), std::begin(result),
                    [&] (const auto& read, const auto& read_hash) {
                        auto mapping_positions = map_query_to_target_helper(read_hash, haplotype_hashes, haplotype_mapping_counts);
                        reset_mapping_counts(haplotype_mapping_counts);
                        return model.evaluate(read, mapping_positions);
                    });
    return result;
}

template <typename Container>
auto calculate_likelihoods(const Genotype<Haplotype>& genotype,
                           const Container& reads,
                           HaplotypeLikelihoodModel& model,
                           OptionalThreadPool workers)
{
    const auto reads_region = encompassing_region(reads);
    const auto read_hashes = compute_read_hashes(reads, workers);
    static constexpr unsigned char mapperKmerSize {6};
    HaplotypeLikelihoods result {};
    const auto indel_factor = estimate_max_indel_size(genotype) + estimate_max_indel_size(reads);
    if (workers && workers->size() > 1 && genotype.ploidy() > 1 && reads.size() > 3) {
        result.resize(genotype.ploidy());
        using octopus::transform;
        transform(std::cbegin(genotype), std::cend(genotype), std::begin(result),
                  [&reads, &reads_region, &read_hashes, indel_factor, model] (const auto& haplotype) mutable {
                      const auto expanded_haplotype = expand_for_alignment(haplotype, reads_region, indel_factor, model);
                      const auto haplotype_hashes = make_kmer_hash_table<mapperKmerSize>(expanded_haplotype.sequence());
                      auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
                      model.reset(expanded_haplotype);
                      return calculate_likelihoods(haplotype, reads, haplotype_hashes, read_hashes, haplotype_mapping_counts, model);
                  }, workers);
    } else {
        result.reserve(genotype.ploidy());
        auto haplotype_hashes = init_kmer_hash_table<mapperKmerSize>();
        for (const auto& haplotype : genotype) {
            const auto expanded_haplotype = expand_for_alignment(haplotype, reads_region, indel_factor, model);
            populate_kmer_hash_table<mapperKmerSize>(expanded_haplotype.sequence(), haplotype_hashes);
            auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
            model.reset(expanded_haplotype);
            auto likelihoods = calculate_likelihoods(haplotype, reads, haplotype_hashes, read_hashes, haplotype_mapping_counts, model);
            clear_kmer_hash_table(haplotype_hashes);
            result.push_back(std::move(likelihoods));
        }
    }
    return result;
}

template <typename ReadType, typename AmbiguousReadListType>
auto
compute_haplotype_support_helper2(const Genotype<Haplotype>& genotype,
                                  const std::vector<ReadType>& reads,
                                  const HaplotypeProbabilityMap& log_priors,
                                  HaplotypeLikelihoodModel model,
                                  boost::optional<AmbiguousReadListType&> ambiguous,
                                  AssignmentConfig config,
                                  OptionalThreadPool workers)
{
    assert(genotype.ploidy() > 1);
    const auto priors = get_priors(genotype, log_priors);
    const auto likelihoods = calculate_likelihoods(genotype, reads, model, workers);
    return calculate_support(genotype, reads, priors, likelihoods, ambiguous, config);
}

template <typename ReadType, typename AmbiguousReadListType>
auto
compute_haplotype_support_helper(const Genotype<Haplotype>& genotype,
                                 const std::vector<ReadType>& reads,
                                 const HaplotypeProbabilityMap& log_priors,
                                 HaplotypeLikelihoodModel model,
                                 boost::optional<AmbiguousReadListType&> ambiguous,
                                 AssignmentConfig config,
                                 OptionalThreadPool workers)
{
    if (is_max_zygosity(genotype)) {
        return compute_haplotype_support_helper2(genotype, reads, log_priors, std::move(model), ambiguous, std::move(config), workers);
    } else {
        return compute_haplotype_support_helper2(collapse(genotype), reads, log_priors, std::move(model), ambiguous, std::move(config), workers);
    }
}

} // namespace

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          const HaplotypeProbabilityMap& log_priors,
                          HaplotypeLikelihoodModel model,
                          boost::optional<AmbiguousReadList&> ambiguous,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    if (!reads.empty()) {
        if (is_heterozygous(genotype)) {
            return compute_haplotype_support_helper(genotype, reads, log_priors, std::move(model), ambiguous, std::move(config), workers);
        } else if (config.ambiguous_action != AssignmentConfig::AmbiguousAction::drop) {
            HaplotypeSupportMap result {};
            result.emplace(genotype[0], reads);
            return result;
        }
    }
    return {};
}

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          const HaplotypeProbabilityMap& log_priors,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    return compute_haplotype_support(genotype, reads, log_priors, model, ambiguous, config, workers);
}

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          HaplotypeLikelihoodModel model,
                          AmbiguousReadList& ambiguous,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    return compute_haplotype_support(genotype, reads, ambiguous, {}, model, config, workers);
}

static HaplotypeLikelihoodModel make_default_haplotype_likelihood_model()
{
    HaplotypeLikelihoodModel::Config config {};
    config.max_indel_error = 8;
    config.use_flank_state = false;
    config.use_mapping_quality = false;
    return {config};
}

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          const HaplotypeProbabilityMap& log_priors,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    auto model = make_default_haplotype_likelihood_model();
    return compute_haplotype_support(genotype, reads, log_priors, model, ambiguous, config, workers);
}

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    auto model = make_default_haplotype_likelihood_model();
    return compute_haplotype_support(genotype, reads, model, ambiguous, config, workers);
}

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          const HaplotypeProbabilityMap& log_priors,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    auto model = make_default_haplotype_likelihood_model();
    return compute_haplotype_support(genotype, reads, log_priors, model, boost::none, config, workers);
}

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    auto model = make_default_haplotype_likelihood_model();
    return compute_haplotype_support(genotype, reads, model, config, workers);
}

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    return compute_haplotype_support(genotype, reads, {}, std::move(model), boost::none, config, workers);
}

HaplotypeSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedRead>& reads,
                          AmbiguousReadList& ambiguous,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    return compute_haplotype_support(genotype, reads, std::move(model), ambiguous, config, workers);
}

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          const HaplotypeProbabilityMap& log_priors,
                          HaplotypeLikelihoodModel model,
                          boost::optional<AmbiguousTemplateList&> ambiguous,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    if (!reads.empty()) {
        if (is_heterozygous(genotype)) {
            return compute_haplotype_support_helper(genotype, reads, log_priors, std::move(model), ambiguous, std::move(config), workers);
        } else if (config.ambiguous_action != AssignmentConfig::AmbiguousAction::drop) {
            HaplotypeTemplateSupportMap result {};
            result.emplace(genotype[0], reads);
            return result;
        }
    }
    return {};
}

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          AmbiguousTemplateList& ambiguous,
                          const HaplotypeProbabilityMap& log_priors,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    return compute_haplotype_support(genotype, reads, log_priors, model, ambiguous, config, workers);
}

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          HaplotypeLikelihoodModel model,
                          AmbiguousTemplateList& ambiguous,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    return compute_haplotype_support(genotype, reads, ambiguous, {}, model, config, workers);
}

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          AmbiguousTemplateList& ambiguous,
                          const HaplotypeProbabilityMap& log_priors,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    auto model = make_default_haplotype_likelihood_model();
    return compute_haplotype_support(genotype, reads, log_priors, model, ambiguous, config, workers);
}

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          AmbiguousTemplateList& ambiguous,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    auto model = make_default_haplotype_likelihood_model();
    return compute_haplotype_support(genotype, reads, model, ambiguous, config, workers);
}

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          const HaplotypeProbabilityMap& log_priors,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    auto model = make_default_haplotype_likelihood_model();
    return compute_haplotype_support(genotype, reads, log_priors, model, boost::none, config, workers);
}

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    auto model = make_default_haplotype_likelihood_model();
    return compute_haplotype_support(genotype, reads, model, config, workers);
}

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    return compute_haplotype_support(genotype, reads, {}, std::move(model), boost::none, config, workers);
}

HaplotypeTemplateSupportMap
compute_haplotype_support(const Genotype<Haplotype>& genotype,
                          const std::vector<AlignedTemplate>& reads,
                          AmbiguousTemplateList& ambiguous,
                          HaplotypeLikelihoodModel model,
                          AssignmentConfig config,
                          OptionalThreadPool workers)
{
    return compute_haplotype_support(genotype, reads, std::move(model), ambiguous, config, workers);
}

bool is_redundant_allele(const Allele& allele) noexcept
{
    return is_empty_region(allele) && is_sequence_empty(allele);
}

bool default_includes(const Haplotype& haplotype, const Allele& allele)
{
    return haplotype.includes(allele) || (is_redundant_allele(allele) && (haplotype.contains(Allele {expand_rhs(mapped_region(allele), 1), ""})));
}

AlleleSupportMap
compute_allele_support(const std::vector<Allele>& alleles,
                       const HaplotypeSupportMap& haplotype_support,
                       std::function<bool(const Haplotype&, const Allele&)> inclusion_pred)
{
    AlleleSupportMap result {};
    result.reserve(alleles.size());
    for (const auto& allele : alleles) {
        ReadRefSupportSet allele_support {};
        for (const auto& p : haplotype_support) {
            if (inclusion_pred(p.first, allele)) {
                allele_support.reserve(p.second.size());
                std::copy_if(std::cbegin(p.second), std::cend(p.second), std::back_inserter(allele_support),
                             [&allele] (const auto& read) { return overlaps(read, allele); });
            }
        }
        std::sort(std::begin(allele_support), std::end(allele_support));
        result.emplace(allele, std::move(allele_support));
    }
    return result;
}

template <typename BinaryPredicate>
auto copy_included(const std::vector<Allele>& alleles, 
                   const Haplotype& haplotype, 
                   const BinaryPredicate& includes_pred)
{
    std::vector<Allele> result {};
    result.reserve(alleles.size());
    std::copy_if(std::cbegin(alleles), std::cend(alleles), std::back_inserter(result),
                 [&] (const auto& allele) { return includes_pred(haplotype, allele); });
    return result;
}

template <typename BinaryPredicate>
struct HaveDifferentAlleles
{
    bool operator()(const std::shared_ptr<Haplotype>& lhs, const std::shared_ptr<Haplotype>& rhs) const
    {
        const auto lhs_includes = copy_included(alleles, *lhs, includes_pred);
        const auto rhs_includes = copy_included(alleles, *rhs, includes_pred);
        return lhs_includes != rhs_includes;
    }
    HaveDifferentAlleles(const std::vector<Allele>& alleles, BinaryPredicate includes_pred) 
    : alleles {alleles}, includes_pred {includes_pred} {}
    const std::vector<Allele>& alleles;
    BinaryPredicate includes_pred;
};

template <typename BinaryPredicate>
bool have_common_alleles(const std::vector<std::shared_ptr<Haplotype>>& haplotypes, 
                         const std::vector<Allele>& alleles, 
                         const BinaryPredicate& includes_pred)
{
    return std::adjacent_find(std::cbegin(haplotypes), std::cend(haplotypes), HaveDifferentAlleles<BinaryPredicate> {alleles, includes_pred}) == std::cend(haplotypes);
}

void sort_and_merge(std::deque<AlignedReadConstReference>& src, ReadRefSupportSet& dst)
{
    std::sort(std::begin(src), std::end(src));
    auto itr = dst.insert(std::end(dst), std::begin(src), std::end(src));
    std::inplace_merge(std::begin(dst), itr, std::end(dst));
}

std::size_t
try_assign_ambiguous_reads_to_alleles(const std::vector<Allele>& alleles,
                                      const AmbiguousReadList& ambiguous_reads,
                                      AlleleSupportMap& allele_support,
                                      std::function<bool(const Haplotype&, const Allele&)> inclusion_pred)
{
    std::size_t num_assigned {0};
    std::unordered_map<Allele, std::deque<AlignedReadConstReference>> assigned {};
    assigned.reserve(alleles.size());
    for (const auto& ambiguous_read : ambiguous_reads) {
        if (ambiguous_read.haplotypes && have_common_alleles(*ambiguous_read.haplotypes, alleles, inclusion_pred)) {
            const auto supported_alleles = copy_included(alleles, *ambiguous_read.haplotypes->front(), inclusion_pred);
            for (const auto& allele : supported_alleles) {
                if (overlaps(ambiguous_read, allele)) {
                    assigned[allele].emplace_back(ambiguous_read.read);
                }
            }
        }
    }
    for (auto& p : assigned) sort_and_merge(p.second, allele_support[p.first]);
    return num_assigned;
}

AlleleSupportMap
compute_allele_support(const std::vector<Allele>& alleles,
                       const HaplotypeSupportMap& haplotype_support,
                       const AmbiguousReadList& ambiguous_reads,
                       std::function<bool(const Haplotype&, const Allele&)> inclusion_pred)
{
    auto result = compute_allele_support(alleles, haplotype_support, inclusion_pred);
    try_assign_ambiguous_reads_to_alleles(alleles, ambiguous_reads, result, inclusion_pred);
    return result;
}

} // namespace octopus
