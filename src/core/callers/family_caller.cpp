// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "family_caller.hpp"

#include <typeinfo>
#include <functional>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <map>
#include <utility>
#include <limits>

#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "containers/probability_matrix.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/calls/germline_variant_call.hpp"
#include "core/types/calls/denovo_call.hpp"
#include "core/types/calls/denovo_reference_reversion_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "exceptions/unimplemented_feature_error.hpp"
#include "core/models/genotype/population_prior_model.hpp"
#include "core/models/genotype/coalescent_population_prior_model.hpp"
#include "core/models/genotype/inheritance_model.hpp"
#include "utils/append.hpp"

namespace octopus {

namespace {

class BadPloidy : public UnimplementedFeatureError
{
    std::string do_help() const override
    {
        return "Use the population caller and/or submit a feature request";
    }
public:
    BadPloidy(unsigned max_ploidy)
    : UnimplementedFeatureError {"family calling with ploidies greater than " + std::to_string(max_ploidy), "FamilyCaller"}
    {}
};

template <typename T>
std::size_t index_of(const T& value, const std::vector<T>& values)
{
    auto itr = std::find(std::cbegin(values), std::cend(values), value);
    return std::distance(std::cbegin(values), itr);
}

} // namespace

FamilyCaller::FamilyCaller(Caller::Components&& components,
                           Caller::Parameters general_parameters,
                           Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {std::move(specific_parameters)}
{
    unique_ploidies_ = parameters_.ploidies;
    std::sort(std::begin(unique_ploidies_), std::end(unique_ploidies_));
    unique_ploidies_.erase(std::unique(std::begin(unique_ploidies_), std::end(unique_ploidies_)), std::end(unique_ploidies_));
    if (unique_ploidies_.back() > 2) {
        throw BadPloidy {2};
    }
    if (unique_ploidies_.back() == 0) {
        throw std::logic_error {"At least one sample must have positive ploidy"};
    }
    offspring_with_two_parents_.reserve(samples_.size());
    offspring_with_one_parent_.reserve(samples_.size());
    for (std::size_t s {0}; s < samples_.size(); ++s) {
        const auto num_parents = parameters_.family.num_parents(samples_[s]);
        if (num_parents == 2) {
            auto mother_idx = index_of(*parameters_.family.mother_of(samples_[s]), samples_);
            auto father_idx = index_of(*parameters_.family.father_of(samples_[s]), samples_);
            offspring_with_two_parents_.emplace_back(s, std::make_pair(mother_idx, father_idx));
        } else if (num_parents == 1) {
            auto parent = parameters_.family.mother_of(samples_[s]);
            if (!parent) parent = parameters_.family.father_of(samples_[s]);
            assert(parent);
            offspring_with_one_parent_.emplace_back(s, index_of(*parent, samples_));
        }
    }
    offspring_with_two_parents_.shrink_to_fit();
    offspring_with_one_parent_.shrink_to_fit();
}

std::string FamilyCaller::do_name() const
{
    return "family";
}

Caller::CallTypeSet FamilyCaller::do_call_types() const
{
    return {std::type_index(typeid(GermlineVariantCall)),
            std::type_index(typeid(DenovoCall)),
            std::type_index(typeid(DenovoReferenceReversionCall))};
}

unsigned FamilyCaller::do_min_callable_ploidy() const
{
    return unique_ploidies_.front();
}

unsigned FamilyCaller::do_max_callable_ploidy() const
{
    return unique_ploidies_.back();
}

std::size_t FamilyCaller::do_remove_duplicates(HaplotypeBlock& haplotypes) const
{
    if (parameters_.deduplicate_haplotypes_with_germline_model) {
        if (haplotypes.size() < 2) return 0;
        CoalescentModel::Parameters model_params {};
        if (parameters_.population_prior_model_params) model_params = *parameters_.population_prior_model_params;
        Haplotype reference {mapped_region(haplotypes), reference_.get()};
        CoalescentModel model {std::move(reference), model_params, haplotypes.size(), CoalescentModel::CachingStrategy::none};
        const CoalescentProbabilityGreater cmp {std::move(model)};
        return octopus::remove_duplicates(haplotypes, cmp);
    } else {
        return Caller::do_remove_duplicates(haplotypes);
    }
}

std::unique_ptr<Caller::Latents>
FamilyCaller::infer_latents(const HaplotypeBlock& haplotypes,
                            const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto indexed_haplotypes = index(haplotypes);
    const auto prior_model = make_prior_model(haplotypes);
    const model::FamilyModel model {parameters_.family, *prior_model, {parameters_.max_genotype_combinations}, debug_log_};
    prior_model->prime(haplotypes);
    if (unique_ploidies_.size() == 1) {
        auto genotypes = generate_all_genotypes(indexed_haplotypes, parameters_.ploidies.front());
        if (debug_log_) stream(*debug_log_) << "There are " << genotypes.size() << " candidate genotypes";
        auto inferences = model.evaluate(samples_, haplotypes, genotypes, haplotype_likelihoods);
        return std::make_unique<Latents>(samples_, parameters_.family, indexed_haplotypes, std::move(genotypes), std::move(inferences));
    } else {
        model::FamilyModel::GenotypeVector genotypes {};
        for (const auto ploidy : unique_ploidies_) {
            if (ploidy > 0) {
                utils::append(generate_all_genotypes(indexed_haplotypes, ploidy), genotypes);
            } else {
                genotypes.push_back(Genotype<IndexedHaplotype<>> {});
            }
        }
        auto inferences = model.evaluate(samples_, parameters_.ploidies, haplotypes, genotypes, haplotype_likelihoods);
        return nullptr;
        // return std::make_unique<Latents>(samples_, parameters_.family, indexed_haplotypes, std::move(genotypes), std::move(inferences));
    }
}

std::vector<std::unique_ptr<VariantCall>>
FamilyCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

namespace {

using JointProbabilityVector = model::FamilyModel::Latents::JointProbabilityVector;
using GenotypeCombination = model::FamilyModel::Latents::GenotypeCombination;
using GenotypeBlock = MappableBlock<Genotype<IndexedHaplotype<>>>;

struct MarginalisationInfo
{
    const GenotypeBlock& genotypes;
    const std::vector<IndexedHaplotype<>>& haplotypes;
};

bool contains_helper(const Haplotype& haplotype, const Allele& allele)
{
    if (!is_indel(allele)) {
        return haplotype.contains(allele);
    } else {
        return haplotype.includes(allele);
    }
}

bool contains_helper(const Genotype<IndexedHaplotype<>>& genotype, const Allele& allele)
{
    if (!is_indel(allele)) {
        return contains(genotype, allele);
    } else {
        return includes(genotype, allele);
    }
}

bool contains(const GenotypeCombination& combination, const Allele& allele, const MarginalisationInfo& info)
{
    return std::any_of(std::cbegin(combination), std::cend(combination), [&] (const auto genotype_index) { 
        return contains_helper(info.genotypes[genotype_index], allele); });
}

using HaplotypeBoolCache = std::vector<boost::optional<bool>>;
using GenotypeCountCache = std::vector<boost::optional<unsigned>>;

// allele posterior calculation

bool contains(const IndexedHaplotype<>& haplotype, const Allele& allele, HaplotypeBoolCache& cache)
{
    if (index_of(haplotype) >= cache.size()) {
        cache.resize(2 * (index_of(haplotype) + 1));
    }
    if (!cache[index_of(haplotype)]) {
        cache[index_of(haplotype)] = contains_helper(haplotype, allele);
    }
    return *cache[index_of(haplotype)];
}

unsigned count_occurrences(const Allele& allele, const Genotype<IndexedHaplotype<>>& genotype,
                           HaplotypeBoolCache& cache)
{
	return std::count_if(std::cbegin(genotype), std::cend(genotype),
                        [&] (const auto& haplotype) { return contains(haplotype, allele, cache); });
}

bool contains(const Genotype<IndexedHaplotype<>>& genotype, const Allele& allele, HaplotypeBoolCache& cache)
{
    return std::any_of(std::cbegin(genotype), std::cend(genotype),
                      [&] (const auto& haplotype) { return contains(haplotype, allele, cache); });
}

auto count_occurrences(const Allele& allele, 
                       const std::size_t genotype_index,
                       HaplotypeBoolCache& haplotype_cache,
                       GenotypeCountCache& genotype_cache,
                       const MarginalisationInfo& info)
{
    if (genotype_index >= genotype_cache.size()) {
        genotype_cache.resize(2 * (genotype_index + 1));
    }
    if (!genotype_cache[genotype_index]) {
        genotype_cache[genotype_index] = count_occurrences(allele, info.genotypes[genotype_index], haplotype_cache);
    }
    return *genotype_cache[genotype_index];
}

bool contains(const std::size_t genotype_index,
              const Allele& allele,
              HaplotypeBoolCache& haplotype_cache,
              GenotypeCountCache& genotype_cache,
              const MarginalisationInfo& info)
{
    return count_occurrences(allele, genotype_index, haplotype_cache, genotype_cache, info) > 0;
}

bool contains(const GenotypeCombination& combination, 
              const Allele& allele,
              HaplotypeBoolCache& haplotype_cache,
              GenotypeCountCache& genotype_cache,
              const MarginalisationInfo& info)
{
    return std::any_of(std::cbegin(combination), std::cend(combination), [&] (const auto genotype_index) { 
        return contains(genotype_index, allele, haplotype_cache, genotype_cache, info); });
}

template <typename UnaryPredicate>
Phred<double>
marginalise_condition(const JointProbabilityVector& genotype_posteriors, UnaryPredicate&& pred)
{
    thread_local std::vector<double> buffer {};
    buffer.clear();
    for (const auto& p : genotype_posteriors) {
        if (!pred(p.combination)) {
            buffer.push_back(p.log_probability);
        }
    }
    if (!buffer.empty()) {
        return log_probability_false_to_phred(std::min(maths::log_sum_exp(buffer), 0.0));
    } else {
        return Phred<double> {std::numeric_limits<double>::infinity()};
    }
}

auto compute_segregation_posterior_uncached(const Allele& allele, const JointProbabilityVector& genotype_posteriors, const MarginalisationInfo& info)
{
    return marginalise_condition(genotype_posteriors, [&] (const auto& combination) { return contains(combination, allele, info); });
}

auto compute_segregation_posterior_cached(const Allele& allele, const JointProbabilityVector& genotype_posteriors, 
                                          const MarginalisationInfo& info)
{
    HaplotypeBoolCache haplotype_cache(info.haplotypes.size());
    GenotypeCountCache genotype_cache(info.genotypes.size());
    return marginalise_condition(genotype_posteriors, [&] (const auto& combination) { return contains(combination, allele, haplotype_cache, genotype_cache, info); });
}

auto compute_segregation_posterior(const Allele& allele, const JointProbabilityVector& genotype_posteriors,
                                   const MarginalisationInfo& info)
{
    if (genotype_posteriors.size() >= 500) {
        return compute_segregation_posterior_cached(allele, genotype_posteriors, info);
    } else {
        return compute_segregation_posterior_uncached(allele, genotype_posteriors, info);
    }
}

auto compute_segregation_posteriors(const std::vector<Allele>& alleles, const JointProbabilityVector& genotype_posteriors,
                                    const MarginalisationInfo& info)
{
    std::vector<Phred<double>> result(alleles.size());
    std::transform(std::cbegin(alleles), std::cend(alleles), std::begin(result),
                   [&] (const auto& allele) { return compute_segregation_posterior(allele, genotype_posteriors, info); });
    return result;
}

auto call_indices(const std::vector<Phred<double>>& posteriors, const Phred<double> min_posterior)
{
    std::vector<std::size_t> result {};
    result.reserve(posteriors.size());
    for (std::size_t idx {0}; idx < posteriors.size(); ++idx) {
        if (posteriors[idx] >= min_posterior) {
            result.push_back(idx);
        }
    }
    return result;
}

struct TrioIndices
{
    std::size_t offspring, mother, father;
};
struct DuoIndices
{
    std::size_t offspring, parent;
};

unsigned count_occurrences(const Allele& allele, const Genotype<IndexedHaplotype<>>& genotype)
{
	return std::count_if(std::cbegin(genotype), std::cend(genotype),
                        [&] (const auto& haplotype) { return contains_helper(haplotype, allele); });
}

struct TrioGenotypeIndices
{
    std::size_t offspring, mother, father;
};
struct DuoGenotypeIndices
{
    std::size_t offspring, parent;
};

TrioGenotypeIndices get_genotype_indices(const TrioIndices& trio, const GenotypeCombination& combination) noexcept
{
    return {combination[trio.offspring], combination[trio.mother], combination[trio.father]};
}
DuoGenotypeIndices get_genotype_indices(const DuoIndices& trio, const GenotypeCombination& combination) noexcept
{
    return {combination[trio.offspring], combination[trio.parent]};
}

bool is_denovo(const Allele& allele, const TrioGenotypeIndices& trio, const MarginalisationInfo& info)
{
	const auto offspring_occurrences = count_occurrences(allele, info.genotypes[trio.offspring]);
	switch(offspring_occurrences) {
		case 0: return false;
		case 1: return !(contains_helper(info.genotypes[trio.mother], allele)
                		|| contains_helper(info.genotypes[trio.father], allele));
		case 2: return !(contains_helper(info.genotypes[trio.mother], allele)
						&& contains_helper(info.genotypes[trio.father], allele));
		default: {
			auto maternal_occurrences = count_occurrences(allele, info.genotypes[trio.mother]);
			auto paternal_occurrences = count_occurrences(allele, info.genotypes[trio.father]);
			return maternal_occurrences > 0 && paternal_occurrences > 0 && (maternal_occurrences + paternal_occurrences) >= offspring_occurrences;
		}
	}
}

bool is_denovo(const Allele& allele, 
               const TrioGenotypeIndices& trio,
               HaplotypeBoolCache& haplotype_cache,
               GenotypeCountCache& genotype_cache,
               const MarginalisationInfo& info)
{
	const auto offspring_occurrences = count_occurrences(allele, trio.offspring, haplotype_cache, genotype_cache, info);
	switch(offspring_occurrences) {
		case 0: return false;
		case 1: return !(contains(trio.mother, allele, haplotype_cache, genotype_cache, info)
                		|| contains(trio.father, allele, haplotype_cache, genotype_cache, info));
		case 2: return !(contains(trio.mother, allele, haplotype_cache, genotype_cache, info)
						&& contains(trio.father, allele, haplotype_cache, genotype_cache, info));
		default: {
			auto maternal_occurrences = count_occurrences(allele, trio.mother, haplotype_cache, genotype_cache, info);
			auto paternal_occurrences = count_occurrences(allele, trio.father, haplotype_cache, genotype_cache, info);
			return maternal_occurrences > 0 && paternal_occurrences > 0 && (maternal_occurrences + paternal_occurrences) >= offspring_occurrences;
		}
	}
}

template <typename TrioOrDuo>
auto compute_denovo_posterior_uncached(const Allele& allele, 
                                       const TrioOrDuo& samples,
                                       const JointProbabilityVector& genotype_posteriors, 
                                       const MarginalisationInfo& info)
{
    return marginalise_condition(genotype_posteriors, [&] (const auto& combination) { 
        return is_denovo(allele, get_genotype_indices(samples, combination), info); });
}

template <typename TrioOrDuo>
auto compute_denovo_posterior_cached(const Allele& allele, 
                                     const TrioOrDuo& samples,
                                     const JointProbabilityVector& genotype_posteriors,
                                     const MarginalisationInfo& info)
{
    HaplotypeBoolCache haplotype_cache(info.haplotypes.size());
    GenotypeCountCache genotype_cache(info.genotypes.size());
    return marginalise_condition(genotype_posteriors, [&] (const auto& combination) { 
        return is_denovo(allele, get_genotype_indices(samples, combination), haplotype_cache, genotype_cache, info); });
}

template <typename TrioOrDuo>
auto compute_denovo_posterior(const Allele& allele, 
                              const TrioOrDuo& samples,
                              const JointProbabilityVector& genotype_posteriors,
                              const MarginalisationInfo& info)
{
    if (genotype_posteriors.size() >= 500) {
        return compute_denovo_posterior_cached(allele, samples, genotype_posteriors, info);
    } else {
        return compute_denovo_posterior_uncached(allele, samples, genotype_posteriors, info);
    }
}

boost::optional<std::size_t> find_variant_index(const Allele& allele, const std::vector<Variant>& variants)
{
    const auto er = std::equal_range(std::cbegin(variants), std::cend(variants), allele,
                                     [] (const auto& lhs, const auto& rhs) { return mapped_region(lhs) < mapped_region(rhs); });
    const auto itr = std::find_if(er.first, er.second, [&allele] (const Variant& v) { return v.alt_allele() == allele; });
    if (itr != er.second) {
        return std::distance(std::cbegin(variants), itr);
    } else {
        return boost::none;
    }
}

auto call_genotypes(const JointProbabilityVector& joint_genotype_posteriors)
{
    assert(!joint_genotype_posteriors.empty());
    const static auto log_probability_less = [] (const auto& lhs, const auto& rhs) { return lhs.log_probability < rhs.log_probability; };
    const auto map_itr = std::max_element(std::cbegin(joint_genotype_posteriors), std::cend(joint_genotype_posteriors), log_probability_less);
    return map_itr->combination;
}

using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<IndexedHaplotype<>>>;

auto call_genotype(const GenotypeProbabilityMap::InnerMap& genotype_posteriors)
{
    return std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                            [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; })->first;
}

using GenotypeVector = std::vector<Genotype<IndexedHaplotype<>>>;

auto call_genotypes(const std::vector<SampleName>& samples, const GenotypeProbabilityMap& genotype_posteriors)
{
    GenotypeVector result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.push_back(call_genotype(genotype_posteriors[sample]));
    }
    return result;
}

auto get_genotypes(const GenotypeCombination& combination, const GenotypeBlock& genotypes)
{
    GenotypeVector result {};
    result.reserve(combination.size());
    for (const auto idx : combination) {
        result.push_back(genotypes[idx]);
    }
    return result;
}

auto call_genotypes(const std::vector<SampleName>& samples, 
                    const JointProbabilityVector& joint_genotype_posteriors, 
                    const GenotypeProbabilityMap& marginal_genotype_posteriors, 
                    const GenotypeBlock& genotypes,
                    const bool use_marginals = false)
{
    if (use_marginals) {
        return call_genotypes(samples, marginal_genotype_posteriors);
    } else {
        return get_genotypes(call_genotypes(joint_genotype_posteriors), genotypes);
    }
}

bool includes(const GenotypeVector& genotypes, const Allele& allele)
{
    const auto includes_allele = [&] (const auto& genotype) { return includes(genotype, allele); };
    return std::any_of(std::cbegin(genotypes), std::cend(genotypes), includes_allele);
}

auto compute_posterior(const Genotype<Allele>& genotype, const GenotypeProbabilityMap::InnerMap& posteriors)
{
    auto p = std::accumulate(std::cbegin(posteriors), std::cend(posteriors), 0.0,
                             [&] (const double curr, const auto& p) {
                                 return curr + (contains(p.first, genotype) ? 0.0 : p.second);
                             });
    return probability_false_to_phred(p);
}

struct GenotypePosterior
{
    Genotype<Allele> genotype;
    Phred<double> posterior;
};

auto marginalise(const Genotype<Allele>& genotype, const GenotypeProbabilityMap::InnerMap& genotype_posteriors)
{
    double mass_not_contained {0};
    if (genotype.ploidy() > 0) {
        for_each_contains(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), genotype,
                          [&] (const auto& p, bool contain) { if (!contain) mass_not_contained += p.second; },
                          [] (const auto& p) { return p.first; });
    }
    return probability_false_to_phred(mass_not_contained);
}

auto marginalise(const GenomicRegion& region,
                 const GenotypeProbabilityMap::InnerMap& posteriors, 
                 const Genotype<IndexedHaplotype<>>& genotype)
{
    GenotypePosterior result {};
    result.genotype = copy<Allele>(genotype, region);
    result.posterior = marginalise(result.genotype, posteriors);
    return result;
}

octopus::VariantCall::GenotypeCall convert(GenotypePosterior&& call)
{
    return octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
}

auto convert_to_call_genotypes(const std::vector<SampleName>& samples, std::vector<GenotypePosterior> genotypes_calls)
{
    std::vector<std::pair<SampleName, Call::GenotypeCall>> result {};
    result.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples),
                   std::make_move_iterator(std::begin(genotypes_calls)),
                   std::back_inserter(result),
                   [] (const auto& sample, auto&& genotype) {
                       return std::make_pair(sample, convert(std::move(genotype)));
                   });
    return result;
}

std::unique_ptr<octopus::VariantCall>
make_germline_call(Variant variant, Phred<double> posterior, 
                  const std::vector<SampleName>& samples, 
                  std::vector<GenotypePosterior> genotypes_calls)
{
    return std::make_unique<GermlineVariantCall>(std::move(variant), convert_to_call_genotypes(samples, genotypes_calls), posterior);
}

auto max(const std::vector<boost::optional<Phred<double>>>& phreds)
{
    boost::optional<Phred<double>> result {};
    for (const auto& phred : phreds) {
        if (phred) {
            if (!result || *result < *phred) {
                result = phred;
            }
        }
    }
    return result;
}

std::unique_ptr<octopus::VariantCall>
make_denovo_call(Variant variant, Phred<double> variant_posterior, Phred<double> denovo_posterior, 
                 const std::vector<SampleName>& samples, 
                 std::vector<GenotypePosterior> genotypes_calls)
{
    return std::make_unique<DenovoCall>(std::move(variant), convert_to_call_genotypes(samples, genotypes_calls), variant_posterior, denovo_posterior);
}

} // namespace

std::vector<std::unique_ptr<VariantCall>>
FamilyCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    const auto alleles = decompose(candidates);
    const auto& joint_genotype_posteriors = latents.inferences_.posteriors.joint_genotype_probabilities;
    const MarginalisationInfo info {latents.genotypes_, latents.haplotypes_};
    const auto allele_posteriors = compute_segregation_posteriors(alleles, joint_genotype_posteriors, info);
    const auto called_allele_indices = call_indices(allele_posteriors, parameters_.min_variant_posterior);
    std::vector<std::vector<boost::optional<Phred<double>>>> denovo_posteriors(called_allele_indices.size(), std::vector<boost::optional<Phred<double>>>(samples_.size()));
    std::vector<std::size_t> called_denovo_allele_indices {};
    called_denovo_allele_indices.reserve(called_allele_indices.size());
    for (const auto& p : offspring_with_two_parents_) {
        const TrioIndices trio {p.first, p.second.first, p.second.second};
        for (const auto allele_idx : called_allele_indices) {
            denovo_posteriors[allele_idx][p.first] = compute_denovo_posterior(alleles[allele_idx], trio, joint_genotype_posteriors, info);
            if (*denovo_posteriors[allele_idx][p.first] >= parameters_.min_denovo_posterior) {
                called_denovo_allele_indices.push_back(allele_idx);
            }
        }
    }
    std::sort(std::begin(called_denovo_allele_indices), std::end(called_denovo_allele_indices));
    called_denovo_allele_indices.erase(std::unique(std::begin(called_denovo_allele_indices), std::end(called_denovo_allele_indices)), 
                                       std::end(called_denovo_allele_indices));
    const auto& marginal_genotype_posteriors = *latents.genotype_posteriors();
    const auto genotype_calls = call_genotypes(samples_, joint_genotype_posteriors, marginal_genotype_posteriors, latents.genotypes_);
    std::vector<std::unique_ptr<VariantCall>> result {};
    for (const auto allele_idx : called_allele_indices) {
        if (includes(genotype_calls, alleles[allele_idx])) {
            auto variant_idx = find_variant_index(alleles[allele_idx], candidates);
            if (variant_idx) {
                const Variant& variant {candidates[*variant_idx]};
                std::vector<GenotypePosterior> variant_genotype_calls {};
                variant_genotype_calls.reserve(samples_.size());
                for (std::size_t sample_idx {0}; sample_idx < samples_.size(); ++sample_idx) {
                    auto variant_genotype = marginalise(mapped_region(variant), marginal_genotype_posteriors[samples_[sample_idx]], genotype_calls[sample_idx]);
                    variant_genotype_calls.push_back(std::move(variant_genotype));
                }
                const auto variant_posterior = allele_posteriors[allele_idx];
                const bool denovo = std::find(std::cbegin(called_denovo_allele_indices), std::cend(called_denovo_allele_indices), allele_idx) != std::cend(called_denovo_allele_indices);
                if (denovo) {
                    auto denovo_posterior = max(denovo_posteriors[allele_idx]);
                    assert(denovo_posterior);
                    result.push_back(make_denovo_call(variant, variant_posterior, *denovo_posterior, samples_, variant_genotype_calls));
                } else {
                    result.push_back(make_germline_call(variant, variant_posterior, samples_, variant_genotype_calls));
                }
            }
        }
    }
    return result;
}

std::vector<std::unique_ptr<ReferenceCall>>
FamilyCaller::call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents,
                             const ReadPileupMap& pileups) const
{
    return call_reference(alleles, dynamic_cast<const Latents&>(latents), pileups);
}

std::vector<std::unique_ptr<ReferenceCall>>
FamilyCaller::call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                             const ReadPileupMap& pileups) const
{
    return {};
}

std::unique_ptr<PedigreePriorModel> 
FamilyCaller::make_prior_model(const HaplotypeBlock& haplotypes) const
{
    std::unique_ptr<PopulationPriorModel> population_prior = std::make_unique<CoalescentPopulationPriorModel>(CoalescentModel{Haplotype {mapped_region(haplotypes), reference_}, *parameters_.population_prior_model_params});
    DeNovoModel denovo_model {parameters_.denovo_model_params, haplotypes.size(), DeNovoModel::CachingStrategy::none};
    InheritanceModel inheritance_prior {std::move(denovo_model)};
    return std::make_unique<PedigreePriorModel>(
        samples_,
        parameters_.family,
        std::move(population_prior),
        std::move(inheritance_prior)
        );
}

namespace {

using InverseGenotypeTable  = std::vector<std::vector<std::size_t>>;

auto make_inverse_genotype_table(const MappableBlock<Genotype<IndexedHaplotype<>>>& genotypes,
                                 const std::size_t num_haplotypes)
{
    InverseGenotypeTable result(num_haplotypes);
    for (auto& indices : result) indices.reserve(genotypes.size() / num_haplotypes);
    for (std::size_t genotype_idx {0}; genotype_idx < genotypes.size(); ++genotype_idx) {
        for (const auto& haplotype : genotypes[genotype_idx]) {
            result[index_of(haplotype)].push_back(genotype_idx);
        }
    }
    for (auto& indices : result) {
        std::sort(std::begin(indices), std::end(indices));
        indices.erase(std::unique(std::begin(indices), std::end(indices)), std::end(indices));
        indices.shrink_to_fit();
    }
    return result;
}

} // namespace

FamilyCaller::Latents::Latents(const std::vector<SampleName>& samples,
                               const Pedigree& family,
                               IndexedHaplotypeBlock haplotypes,
                               Genotypeblock genotypes,
                               ModelInferences latents)
: genotypes_ {std::move(genotypes)}
, haplotypes_ {std::move(haplotypes)}
, inferences_ {std::move(latents)}
{
    std::vector<std::vector<double>> marginal_genotype_posteriors(samples.size(), std::vector<double>(genotypes_.size(), 0.0));
    for (const auto& p : inferences_.posteriors.joint_genotype_probabilities) {
        assert(p.combination.size() == samples.size());
        const auto probability = std::exp(p.log_probability);
        for (std::size_t s {0}; s < samples.size(); ++s) {
            marginal_genotype_posteriors[s][p.combination[s]] += probability;
        }
    }
    GenotypeProbabilityMap genotype_posteriors {std::begin(genotypes_), std::end(genotypes_)};
    for (std::size_t s {0}; s < samples.size(); ++s) {
        insert_sample(samples[s], marginal_genotype_posteriors[s], genotype_posteriors);
    }
    genotype_posteriors_ = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));

    // haplotype posteriors
    const auto inverse_genotypes = make_inverse_genotype_table(genotypes, haplotypes_.size());
    std::vector<double> haplotype_posteriors(haplotypes_.size());
    auto itr = std::cbegin(inverse_genotypes);
    std::vector<std::size_t> genotype_indices(genotypes_.size());
    std::iota(std::begin(genotype_indices), std::end(genotype_indices), 0);
    // noncontaining genotypes are genotypes that do not contain a particular haplotype.
    const auto num_noncontaining_genotypes = genotypes_.size() - itr->size();
    std::vector<std::size_t> noncontaining_genotype_indices(num_noncontaining_genotypes);
    for (const auto& haplotype : haplotypes_) {
        std::set_difference(std::cbegin(genotype_indices), std::cend(genotype_indices),
                            std::cbegin(*itr), std::cend(*itr),
                            std::begin(noncontaining_genotype_indices));
        double prob_not_observed {1};
        for (const auto& sample_genotype_posteriors : marginal_genotype_posteriors) {
            prob_not_observed *= std::accumulate(std::cbegin(noncontaining_genotype_indices), std::cend(noncontaining_genotype_indices),
                                                 0.0, [&] (const auto total, const auto i) {
                return total + sample_genotype_posteriors[i];
            });
        }
        haplotype_posteriors[index_of(haplotype)] = 1.0 - prob_not_observed;
        ++itr;
    }
    haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>();
    haplotype_posteriors_->reserve(haplotypes_.size());
    for (const auto& haplotype : haplotypes_) {
        haplotype_posteriors_->emplace(haplotype, haplotype_posteriors[index_of(haplotype)]);
    }
}

std::shared_ptr<FamilyCaller::Latents::HaplotypeProbabilityMap>
FamilyCaller::Latents::haplotype_posteriors() const noexcept
{
    return haplotype_posteriors_;
}

std::shared_ptr<FamilyCaller::Latents::GenotypeProbabilityMap>
FamilyCaller::Latents::genotype_posteriors() const noexcept
{
    return genotype_posteriors_;
}

} // namespace octopus
