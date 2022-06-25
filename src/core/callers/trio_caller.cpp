// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "trio_caller.hpp"

#include <typeinfo>
#include <functional>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <map>
#include <utility>
#include <limits>

#include <boost/multiprecision/cpp_dec_float.hpp>

#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "containers/probability_matrix.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/calls/germline_variant_call.hpp"
#include "core/types/calls/denovo_call.hpp"
#include "core/types/calls/denovo_reference_reversion_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "core/models/genotype/uniform_population_prior_model.hpp"
#include "core/models/genotype/coalescent_population_prior_model.hpp"
#include "core/models/genotype/population_model.hpp"
#include "core/models/genotype/coalescent_genotype_prior_model.hpp"
#include "core/models/genotype/uniform_genotype_prior_model.hpp"
#include "core/models/genotype/individual_model.hpp"
#include "utils/map_utils.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"
#include "utils/concat.hpp"
#include "utils/select_top_k.hpp"
#include "exceptions/unimplemented_feature_error.hpp"

namespace octopus {

class BadPloidy : public UnimplementedFeatureError
{
    std::string do_help() const override
    {
        return "Use the population caller and/or submit a feature request";
    }
public:
    BadPloidy(unsigned max_ploidy)
    : UnimplementedFeatureError {"trio calling with ploidies greater than " + std::to_string(max_ploidy), "TrioCaller"}
    {}
};

TrioCaller::TrioCaller(Caller::Components&& components,
                       Caller::Parameters general_parameters,
                       Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {std::move(specific_parameters)}
{
    if (parameters_.maternal_ploidy == 0 && parameters_.paternal_ploidy == 0 && parameters_.child_ploidy == 0) {
        throw std::logic_error {"At least one sample must have positive ploidy"};
    }
    if (parameters_.child_ploidy == 0 && parameters_.maternal_ploidy > 0 && parameters_.paternal_ploidy > 0) {
        throw std::logic_error {"There must be at least one inherited haplotype if both parents have zygosity"};
    }
    const auto max_ploidy = model::TrioModel::max_ploidy();
    if (parameters_.maternal_ploidy > max_ploidy || parameters_.paternal_ploidy > max_ploidy || parameters_.child_ploidy > max_ploidy) {
        throw BadPloidy {max_ploidy};
    }
}

std::string TrioCaller::do_name() const
{
    return "trio";
}

Caller::CallTypeSet TrioCaller::do_call_types() const
{
    return {std::type_index(typeid(GermlineVariantCall)),
            std::type_index(typeid(DenovoCall)),
            std::type_index(typeid(DenovoReferenceReversionCall))};
}

unsigned TrioCaller::do_min_callable_ploidy() const
{
    return std::max({parameters_.maternal_ploidy, parameters_.paternal_ploidy, parameters_.child_ploidy});
}

unsigned TrioCaller::do_max_callable_ploidy() const
{
    return std::max({parameters_.maternal_ploidy, parameters_.paternal_ploidy, parameters_.child_ploidy});
}

std::size_t TrioCaller::do_remove_duplicates(HaplotypeBlock& haplotypes) const
{
    if (parameters_.deduplicate_haplotypes_with_germline_model) {
        if (haplotypes.size() < 2) return 0;
        CoalescentModel::Parameters model_params {};
        if (parameters_.germline_prior_model_params) model_params = *parameters_.germline_prior_model_params;
        Haplotype reference {mapped_region(haplotypes), reference_.get()};
        CoalescentModel model {std::move(reference), model_params, haplotypes.size(), CoalescentModel::CachingStrategy::none};
        const CoalescentProbabilityGreater cmp {std::move(model)};
        return octopus::remove_duplicates(haplotypes, cmp);
    } else {
        return Caller::do_remove_duplicates(haplotypes);
    }
}

// TrioCaller::Latents

TrioCaller::Latents::Latents(IndexedHaplotypeBlock haplotypes,
                             GenotypeBlock genotypes,
                             model::TrioModel::InferredLatents latents,
                             const Trio& trio)
: trio {trio}
, haplotypes {std::move(haplotypes)}
, maternal_genotypes {std::move(genotypes)}
, model_latents {std::move(latents)}
, concatenated_genotypes_ {}
, padded_marginal_maternal_posteriors_ {}
, padded_marginal_paternal_posteriors_ {}
, padded_marginal_child_posteriors_ {}
, child_ploidy_ {maternal_genotypes.front().ploidy()}
{
    set_genotype_posteriors_shared_genotypes(trio);
    set_haplotype_posteriors_shared_genotypes();
}

TrioCaller::Latents::Latents(IndexedHaplotypeBlock haplotypes,
                             GenotypeBlock maternal_genotypes,
                             GenotypeBlock paternal_genotypes,
                             const unsigned child_ploidy,
                             ModelInferences latents,
                             const Trio& trio)
: trio {trio}
, haplotypes {std::move(haplotypes)}
, maternal_genotypes {std::move(maternal_genotypes)}
, paternal_genotypes {std::move(paternal_genotypes)}
, model_latents {std::move(latents)}
, concatenated_genotypes_ {}
, padded_marginal_maternal_posteriors_ {}
, padded_marginal_paternal_posteriors_ {}
, padded_marginal_child_posteriors_ {}
, child_ploidy_ {child_ploidy}
{
    set_genotype_posteriors_unique_genotypes(trio);
    set_haplotype_posteriors_unique_genotypes();
    padded_marginal_maternal_posteriors_.clear();
    padded_marginal_maternal_posteriors_.shrink_to_fit();
    padded_marginal_paternal_posteriors_.clear();
    padded_marginal_paternal_posteriors_.shrink_to_fit();
    padded_marginal_child_posteriors_.clear();
    padded_marginal_child_posteriors_.shrink_to_fit();
}

std::shared_ptr<TrioCaller::Latents::HaplotypeProbabilityMap>
TrioCaller::Latents::haplotype_posteriors() const noexcept
{
    return marginal_haplotype_posteriors;
}

std::shared_ptr<TrioCaller::Latents::GenotypeProbabilityMap>
TrioCaller::Latents::genotype_posteriors() const noexcept
{
    return marginal_genotype_posteriors;
}

namespace {

using model::TrioModel;
using JointProbability = TrioModel::Latents::JointProbability;

using GenotypeReference = std::reference_wrapper<const Genotype<IndexedHaplotype<>>>;

template <typename Function>
auto marginalise(const std::vector<Genotype<IndexedHaplotype<>>>& genotypes,
                 const std::vector<JointProbability>& joint_posteriors,
                 Function&& who)
{
    std::vector<double> result(genotypes.size(), 0.0);
    if (genotypes.empty()) return result;
    for (const auto& jp : joint_posteriors) {
        result[who(jp)] += jp.probability;
    }
    return result;
}

auto marginalise_mother(const std::vector<Genotype<IndexedHaplotype<>>>& genotypes,
                        const std::vector<JointProbability>& joint_posteriors)
{
    return marginalise(genotypes, joint_posteriors, [] (const JointProbability& p) noexcept { return p.maternal; });
}

auto marginalise_father(const std::vector<Genotype<IndexedHaplotype<>>>& genotypes,
                        const std::vector<JointProbability>& joint_posteriors)
{
    return marginalise(genotypes, joint_posteriors, [] (const JointProbability& p) noexcept { return p.paternal; });
}

auto marginalise_child(const std::vector<Genotype<IndexedHaplotype<>>>& genotypes,
                       const std::vector<JointProbability>& joint_posteriors)
{
    return marginalise(genotypes, joint_posteriors, [] (const JointProbability& p) noexcept { return p.child; });
}

} // namespace

void TrioCaller::Latents::set_genotype_posteriors(const Trio& trio)
{
    if (paternal_genotypes) {
        set_genotype_posteriors_unique_genotypes(trio);
    } else {
        set_genotype_posteriors_shared_genotypes(trio);
    }
}

void TrioCaller::Latents::set_genotype_posteriors_shared_genotypes(const Trio& trio)
{
    auto& trio_posteriors = model_latents.posteriors.joint_genotype_probabilities;
    marginal_maternal_posteriors = marginalise_mother(maternal_genotypes, trio_posteriors);
    marginal_paternal_posteriors = marginalise_father(maternal_genotypes, trio_posteriors);
    marginal_child_posteriors    = marginalise_child(maternal_genotypes, trio_posteriors);
    GenotypeProbabilityMap genotype_posteriors {std::begin(maternal_genotypes), std::end(maternal_genotypes)};
    insert_sample(trio.mother(), marginal_maternal_posteriors, genotype_posteriors);
    insert_sample(trio.father(), marginal_paternal_posteriors, genotype_posteriors);
    insert_sample(trio.child(), marginal_child_posteriors, genotype_posteriors);
    marginal_genotype_posteriors = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
}

void TrioCaller::Latents::set_genotype_posteriors_unique_genotypes(const Trio& trio)
{
    auto& trio_posteriors = model_latents.posteriors.joint_genotype_probabilities;
    assert(paternal_genotypes);
    marginal_paternal_posteriors = marginalise_father(*paternal_genotypes, trio_posteriors);
    if (!maternal_genotypes.empty()) {
        marginal_maternal_posteriors = marginalise_mother(maternal_genotypes, trio_posteriors);
        for (auto& t : trio_posteriors) t.paternal += maternal_genotypes.size();
    } else {
        maternal_genotypes.assign({Genotype<IndexedHaplotype<>> {}});
        marginal_maternal_posteriors.assign({1.0});
        for (auto& t : trio_posteriors) ++t.paternal;
    }
    const bool child_shares_paternal_genotypes {child_ploidy_ == paternal_genotypes->front().ploidy()};
    if (child_shares_paternal_genotypes) {
        marginal_child_posteriors = marginalise_child(*paternal_genotypes, trio_posteriors);
        for (auto& t : trio_posteriors) t.child += maternal_genotypes.size();
    } else {
        if (maternal_genotypes.size() > 1) {
            marginal_child_posteriors = marginalise_child(maternal_genotypes, trio_posteriors);
        } else {
            marginal_child_posteriors.assign({1.0});
        }
    }
    concatenated_genotypes_ = concat(maternal_genotypes, *paternal_genotypes);
    GenotypeProbabilityMap genotype_posteriors {std::cbegin(concatenated_genotypes_), std::cend(concatenated_genotypes_)};
    const auto num_unique_genotypes = concatenated_genotypes_.size();
    padded_marginal_maternal_posteriors_.resize(num_unique_genotypes);
    std::copy(std::cbegin(marginal_maternal_posteriors), std::cend(marginal_maternal_posteriors),
              std::begin(padded_marginal_maternal_posteriors_));
    padded_marginal_paternal_posteriors_.resize(num_unique_genotypes);
    std::copy(std::crbegin(marginal_paternal_posteriors), std::crend(marginal_paternal_posteriors),
              std::rbegin(padded_marginal_paternal_posteriors_));
    padded_marginal_child_posteriors_.resize(num_unique_genotypes);
    if (child_shares_paternal_genotypes) {
        std::copy(std::crbegin(marginal_child_posteriors), std::crend(marginal_child_posteriors),
                  std::rbegin(padded_marginal_child_posteriors_));
    } else {
        std::copy(std::cbegin(marginal_child_posteriors), std::cend(marginal_child_posteriors),
                  std::begin(padded_marginal_child_posteriors_));
    }
    insert_sample(trio.mother(), padded_marginal_maternal_posteriors_, genotype_posteriors);
    insert_sample(trio.father(), padded_marginal_paternal_posteriors_, genotype_posteriors);
    insert_sample(trio.child(), padded_marginal_child_posteriors_, genotype_posteriors);
    marginal_genotype_posteriors = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
}

namespace {

using JointProbability      = TrioModel::Latents::JointProbability;
using TrioProbabilityVector = std::vector<JointProbability>;
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

using GenotypeMarginalPosteriorMatrix = std::vector<std::vector<double>>;

auto calculate_haplotype_posteriors(const MappableBlock<IndexedHaplotype<>>& haplotypes,
                                    const MappableBlock<Genotype<IndexedHaplotype<>>>& genotypes,
                                    const GenotypeMarginalPosteriorMatrix& genotype_posteriors)
{
    const auto inverse_genotypes = make_inverse_genotype_table(genotypes, haplotypes.size());
    std::vector<double> result(haplotypes.size());
    auto itr = std::cbegin(inverse_genotypes);
    std::vector<std::size_t> genotype_indices(genotypes.size());
    std::iota(std::begin(genotype_indices), std::end(genotype_indices), 0);
    // noncontaining genotypes are genotypes that do not contain a particular haplotype.
    const auto num_noncontaining_genotypes = genotypes.size() - itr->size();
    std::vector<std::size_t> noncontaining_genotype_indices(num_noncontaining_genotypes);
    for (const auto& haplotype : haplotypes) {
        std::set_difference(std::cbegin(genotype_indices), std::cend(genotype_indices),
                            std::cbegin(*itr), std::cend(*itr),
                            std::begin(noncontaining_genotype_indices));
        double prob_not_observed {1};
        for (const auto& sample_genotype_posteriors : genotype_posteriors) {
            prob_not_observed *= std::accumulate(std::cbegin(noncontaining_genotype_indices), std::cend(noncontaining_genotype_indices),
                                                 0.0, [&] (const auto total, const auto i) {
                return total + sample_genotype_posteriors[i];
            });
        }
        result[index_of(haplotype)] = 1.0 - prob_not_observed;
        ++itr;
    }
    return result;
}

} // namespace

void TrioCaller::Latents::set_haplotype_posteriors()
{
    if (paternal_genotypes) {
        set_haplotype_posteriors_unique_genotypes();
    } else {
        set_haplotype_posteriors_shared_genotypes();
    }
}

void TrioCaller::Latents::set_haplotype_posteriors_shared_genotypes()
{
    const GenotypeMarginalPosteriorMatrix genotype_posteriors {marginal_maternal_posteriors,
                                                               marginal_paternal_posteriors,
                                                               marginal_child_posteriors};
    auto haplotype_posteriors = calculate_haplotype_posteriors(haplotypes, maternal_genotypes, genotype_posteriors);
    marginal_haplotype_posteriors = std::make_shared<HaplotypeProbabilityMap>();
    marginal_haplotype_posteriors->reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) {
        marginal_haplotype_posteriors->emplace(haplotype, haplotype_posteriors[index_of(haplotype)]);
    }
}

void TrioCaller::Latents::set_haplotype_posteriors_unique_genotypes()
{
    const GenotypeMarginalPosteriorMatrix genotype_posteriors {padded_marginal_maternal_posteriors_,
                                                               padded_marginal_paternal_posteriors_,
                                                               padded_marginal_child_posteriors_};
    auto haplotype_posteriors = calculate_haplotype_posteriors(haplotypes, concatenated_genotypes_, genotype_posteriors);
    marginal_haplotype_posteriors = std::make_shared<HaplotypeProbabilityMap>();
    marginal_haplotype_posteriors->reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) {
        marginal_haplotype_posteriors->emplace(haplotype, haplotype_posteriors[index_of(haplotype)]);
    }
}

// TrioCaller

std::unique_ptr<Caller::Latents>
TrioCaller::infer_latents(const HaplotypeBlock& haplotypes,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods,
                          OptionalThreadPool workers) const
{
    const auto indexed_haplotypes = index(haplotypes);
    if (parameters_.child_ploidy == 0) {
        assert(parameters_.maternal_ploidy == 0 || parameters_.paternal_ploidy == 0);
        const auto prior_model = make_single_sample_prior_model(haplotypes);
        prior_model->prime(haplotypes);
        const model::IndividualModel sample_model {*prior_model};
        GenotypeBlock parent_genotypes {mapped_region(haplotypes)}, empty_genotypes {mapped_region(haplotypes)};
        if (parameters_.maternal_ploidy > 0) {
            parent_genotypes = generate_all_genotypes(indexed_haplotypes, parameters_.maternal_ploidy);
            haplotype_likelihoods.prime(parameters_.trio.mother());
        } else {
            parent_genotypes = generate_all_genotypes(indexed_haplotypes, parameters_.paternal_ploidy);
            haplotype_likelihoods.prime(parameters_.trio.father());
        }
        auto sample_latents = sample_model.evaluate(parent_genotypes, haplotype_likelihoods);
        model::TrioModel::InferredLatents trio_latents {};
        trio_latents.log_evidence = sample_latents.log_evidence;
        trio_latents.posteriors.joint_genotype_probabilities.resize(parent_genotypes.size());
        using GI = TrioModel::Latents::JointProbability::GenotypeIndex;
        for (GI idx {0}; idx < parent_genotypes.size(); ++idx) {
            const auto log_posterior = sample_latents.posteriors.genotype_log_probabilities[idx];
            auto maternal_idx = parameters_.maternal_ploidy > 0 ? idx : GI {0};
            auto paternal_idx = parameters_.paternal_ploidy > 0 ? idx : GI {0};
            trio_latents.posteriors.joint_genotype_probabilities[idx] = {log_posterior, std::exp(log_posterior), maternal_idx, paternal_idx, 0};
        }
        if (parameters_.maternal_ploidy == 0) std::swap(parent_genotypes, empty_genotypes);
        return std::make_unique<Latents>(std::move(indexed_haplotypes), std::move(parent_genotypes), std::move(empty_genotypes),
                                         parameters_.child_ploidy, std::move(trio_latents), parameters_.trio);
    }
    auto germline_prior_model = make_prior_model(haplotypes);
    germline_prior_model->prime(haplotypes);
    DeNovoModel denovo_model {parameters_.denovo_model_params, haplotypes.size(), DeNovoModel::CachingStrategy::none};
    denovo_model.prime(haplotypes);
    const model::TrioModel model {
        parameters_.trio, *germline_prior_model, denovo_model,
        TrioModel::Options {parameters_.max_genotype_combinations},
        debug_log_
    };
    auto maternal_genotypes = generate_all_genotypes(indexed_haplotypes, parameters_.maternal_ploidy);
    if (parameters_.maternal_ploidy == parameters_.paternal_ploidy) {
        auto latents = model.evaluate(maternal_genotypes, haplotype_likelihoods);
        return std::make_unique<Latents>(std::move(indexed_haplotypes), std::move(maternal_genotypes),
                                         std::move(latents), parameters_.trio);
    } else {
        auto paternal_genotypes = generate_all_genotypes(indexed_haplotypes, parameters_.paternal_ploidy);
        if (parameters_.maternal_ploidy == parameters_.child_ploidy) {
            auto latents = model.evaluate(maternal_genotypes, paternal_genotypes,
                                          maternal_genotypes, haplotype_likelihoods);
            return std::make_unique<Latents>(std::move(indexed_haplotypes), std::move(maternal_genotypes), std::move(paternal_genotypes),
                                             parameters_.child_ploidy, std::move(latents), parameters_.trio);
        } else {
            auto latents = model.evaluate(maternal_genotypes, paternal_genotypes,
                                          paternal_genotypes, haplotype_likelihoods);
            return std::make_unique<Latents>(std::move(indexed_haplotypes), std::move(maternal_genotypes), std::move(paternal_genotypes),
                                             parameters_.child_ploidy, std::move(latents), parameters_.trio);
        }
    }
}

boost::optional<Caller::ModelPosterior>
TrioCaller::calculate_model_posterior(const HaplotypeBlock& haplotypes,
                                      const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                      const Caller::Latents& latents) const
{
    return calculate_model_posterior(haplotypes, haplotype_likelihoods, dynamic_cast<const Latents&>(latents));
}

namespace {

static auto calculate_model_posterior(const double normal_model_log_evidence,
                                      const double dummy_model_log_evidence)
{
    constexpr double normalModelPrior {0.9999999};
    constexpr double dummyModelPrior {1.0 - normalModelPrior};
    const auto normal_model_ljp = std::log(normalModelPrior) + normal_model_log_evidence;
    const auto dummy_model_ljp  = std::log(dummyModelPrior) + dummy_model_log_evidence;
    const auto norm = maths::log_sum_exp(normal_model_ljp, dummy_model_ljp);
    return std::exp(normal_model_ljp - norm);
}

} // namespace

boost::optional<Caller::ModelPosterior>
TrioCaller::calculate_model_posterior(const HaplotypeBlock& haplotypes,
                                      const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                      const Latents& latents) const
{
    const auto indexed_haplotypes = index(haplotypes);
    const auto prior_model = make_single_sample_prior_model(haplotypes);
    prior_model->prime(haplotypes);
    const model::IndividualModel model {*prior_model, debug_log_};
    ModelPosterior result {};
    result.samples.resize(3);
    const auto propose_model_check_genotypes_helper = [&] (const std::size_t sample_idx) {
        if (samples_[sample_idx] == parameters_.trio.mother()) {
            return propose_model_check_genotypes(haplotypes, indexed_haplotypes, latents.maternal_genotypes, latents.marginal_maternal_posteriors);
        } else if (samples_[sample_idx] == parameters_.trio.father()) {
            if (latents.paternal_genotypes) {
                return propose_model_check_genotypes(haplotypes, indexed_haplotypes, *latents.paternal_genotypes, latents.marginal_paternal_posteriors);
            } else {
                return propose_model_check_genotypes(haplotypes, indexed_haplotypes, latents.maternal_genotypes, latents.marginal_paternal_posteriors);
            }
        } else {
            if (parameters_.maternal_ploidy == parameters_.child_ploidy) {
                return propose_model_check_genotypes(haplotypes, indexed_haplotypes, latents.maternal_genotypes, latents.marginal_child_posteriors);
            } else {
                return propose_model_check_genotypes(haplotypes, indexed_haplotypes, *latents.paternal_genotypes, latents.marginal_child_posteriors);
            }
        }
    };
    for (std::size_t sample_idx {0}; sample_idx < samples_.size(); ++sample_idx) {
        haplotype_likelihoods.prime(samples_[sample_idx]);
        const auto genotypes = propose_model_check_genotypes_helper(sample_idx);
        const auto inferences1 = model.evaluate(genotypes.first, haplotype_likelihoods);
        const auto inferences2 = model.evaluate(genotypes.second, haplotype_likelihoods);
        result.samples[sample_idx] = octopus::calculate_model_posterior(inferences1.log_evidence, inferences2.log_evidence); 
    }
    // const auto max_ploidy = std::max({parameters_.maternal_ploidy, parameters_.paternal_ploidy, parameters_.child_ploidy});
    // if (max_ploidy + 1 <= model::TrioModel::max_ploidy()) {
    //     const auto genotypes = generate_all_genotypes(indexed_haplotypes, max_ploidy + 1);
    //     const auto germline_prior_model = make_prior_model(haplotypes);
    //     DeNovoModel denovo_model {parameters_.denovo_model_params};
    //     germline_prior_model->prime(haplotypes);
    //     denovo_model.prime(haplotypes);
    //     if (debug_log_) *debug_log_ << "Calculating model posterior";
    //     const model::TrioModel model {parameters_.trio, *germline_prior_model, denovo_model,
    //                                   TrioModel::Options {parameters_.max_genotype_combinations},
    //                                   debug_log_};
    //     const auto inferences = model.evaluate(genotypes, haplotype_likelihoods);
    //     ModelPosterior result {};
    //     result.joint = octopus::calculate_model_posterior(latents.model_latents.log_evidence, inferences.log_evidence);
    //     return result;
    // }
    return result;
}

std::vector<std::unique_ptr<VariantCall>>
TrioCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents, OptionalThreadPool workers) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

namespace {

struct TrioGenotypeInfo
{
    const MappableBlock<Genotype<IndexedHaplotype<>>>& genotypes;
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

bool contains(const JointProbability& trio, const Allele& allele, const TrioGenotypeInfo& info)
{
    return contains_helper(info.genotypes[trio.maternal], allele)
           || contains_helper(info.genotypes[trio.paternal], allele)
           || contains_helper(info.genotypes[trio.child], allele);
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
                       const TrioGenotypeInfo& info)
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
              const TrioGenotypeInfo& info)
{
    return count_occurrences(allele, genotype_index, haplotype_cache, genotype_cache, info) > 0;
}

bool contains(const JointProbability& trio, 
              const Allele& allele,
              HaplotypeBoolCache& haplotype_cache,
              GenotypeCountCache& genotype_cache,
              const TrioGenotypeInfo& info)
{
    return contains(trio.maternal, allele, haplotype_cache, genotype_cache, info)
        || contains(trio.paternal, allele, haplotype_cache, genotype_cache, info)
        || contains(trio.child, allele, haplotype_cache, genotype_cache, info);
}

template <typename UnaryPredicate>
Phred<double>
marginalise_condition(const TrioProbabilityVector& trio_posteriors, UnaryPredicate&& pred)
{
    thread_local std::vector<double> buffer {};
    buffer.clear();
    for (const auto& trio : trio_posteriors) {
        if (!pred(trio)) {
            buffer.push_back(trio.log_probability);
        }
    }
    if (!buffer.empty()) {
        return log_probability_false_to_phred(std::min(maths::log_sum_exp(buffer), 0.0));
    } else {
        return Phred<double> {std::numeric_limits<double>::infinity()};
    }
}

auto compute_segregation_posterior_uncached(const Allele& allele, const TrioProbabilityVector& trio_posteriors, const TrioGenotypeInfo& info)
{
    return marginalise_condition(trio_posteriors, [&] (const auto& trio) { return contains(trio, allele, info); });
}

auto compute_segregation_posterior_cached(const Allele& allele, const TrioProbabilityVector& trio_posteriors, 
                                          const TrioGenotypeInfo& info)
{
    HaplotypeBoolCache haplotype_cache(info.haplotypes.size());
    GenotypeCountCache genotype_cache(info.genotypes.size());
    return marginalise_condition(trio_posteriors, [&] (const auto& trio) { return contains(trio, allele, haplotype_cache, genotype_cache, info); });
}

auto compute_segregation_posterior(const Allele& allele, const TrioProbabilityVector& trio_posteriors,
                                   const TrioGenotypeInfo& info)
{
    if (trio_posteriors.size() >= 500) {
        return compute_segregation_posterior_cached(allele, trio_posteriors, info);
    } else {
        return compute_segregation_posterior_uncached(allele, trio_posteriors, info);
    }
}

using AllelePosteriorMap = std::map<Allele, Phred<double>>;

auto compute_segregation_posteriors(const std::vector<Allele>& alleles, const TrioProbabilityVector& trio_posteriors,
                                    const TrioGenotypeInfo& info)
{
    AllelePosteriorMap result {};
    for (const auto& allele : alleles) {
        result.emplace(allele, compute_segregation_posterior(allele, trio_posteriors, info));
    }
    return result;
}

auto call_alleles(const AllelePosteriorMap& allele_posteriors, const Phred<double> min_posterior)
{
    AllelePosteriorMap result {};
    std::copy_if(std::cbegin(allele_posteriors), std::cend(allele_posteriors),
                 std::inserter(result, std::begin(result)),
                 [=] (const auto& p) { return p.second >= min_posterior; });
    return result;
}

// de novo posterior calculation

unsigned count_occurrences(const Allele& allele, const Genotype<IndexedHaplotype<>>& genotype)
{
	return std::count_if(std::cbegin(genotype), std::cend(genotype),
                        [&] (const auto& haplotype) { return contains_helper(haplotype, allele); });
}

bool is_denovo(const Allele& allele, const JointProbability& trio, const TrioGenotypeInfo& info)
{
	const auto child_occurrences = count_occurrences(allele, info.genotypes[trio.child]);
	switch(child_occurrences) {
		case 0: return false;
		case 1: return !(contains_helper(info.genotypes[trio.maternal], allele)
                		|| contains_helper(info.genotypes[trio.paternal], allele));
		case 2: return !(contains_helper(info.genotypes[trio.maternal], allele)
						&& contains_helper(info.genotypes[trio.paternal], allele));
		default: {
			auto maternal_occurrences = count_occurrences(allele, info.genotypes[trio.maternal]);
			auto paternal_occurrences = count_occurrences(allele, info.genotypes[trio.paternal]);
			return maternal_occurrences > 0 && paternal_occurrences > 0 && (maternal_occurrences + paternal_occurrences) >= child_occurrences;
		}
	}
}

bool is_denovo(const Allele& allele, 
               const JointProbability& trio,
               HaplotypeBoolCache& haplotype_cache,
               GenotypeCountCache& genotype_cache,
               const TrioGenotypeInfo& info)
{
	const auto child_occurrences = count_occurrences(allele, trio.child, haplotype_cache, genotype_cache, info);
	switch(child_occurrences) {
		case 0: return false;
		case 1: return !(contains(trio.maternal, allele, haplotype_cache, genotype_cache, info)
                		|| contains(trio.paternal, allele, haplotype_cache, genotype_cache, info));
		case 2: return !(contains(trio.maternal, allele, haplotype_cache, genotype_cache, info)
						&& contains(trio.paternal, allele, haplotype_cache, genotype_cache, info));
		default: {
			auto maternal_occurrences = count_occurrences(allele, trio.maternal, haplotype_cache, genotype_cache, info);
			auto paternal_occurrences = count_occurrences(allele, trio.paternal, haplotype_cache, genotype_cache, info);
			return maternal_occurrences > 0 && paternal_occurrences > 0 && (maternal_occurrences + paternal_occurrences) >= child_occurrences;
		}
	}
}

auto compute_denovo_posterior_uncached(const Allele& allele, const TrioProbabilityVector& trio_posteriors, 
                                       const TrioGenotypeInfo& info)
{
    return marginalise_condition(trio_posteriors, [&] (const auto& trio) { return is_denovo(allele, trio, info); });
}

auto compute_denovo_posterior_cached(const Allele& allele, const TrioProbabilityVector& trio_posteriors,
                                     const TrioGenotypeInfo& info)
{
    HaplotypeBoolCache haplotype_cache(info.haplotypes.size());
    GenotypeCountCache genotype_cache(info.genotypes.size());
    return marginalise_condition(trio_posteriors, [&] (const auto& trio) { return is_denovo(allele, trio, haplotype_cache, genotype_cache, info); });
}

auto compute_denovo_posterior(const Allele& allele, const TrioProbabilityVector& trio_posteriors,
                              const TrioGenotypeInfo& info)
{
    if (trio_posteriors.size() >= 500) {
        return compute_denovo_posterior_cached(allele, trio_posteriors, info);
    } else {
        return compute_denovo_posterior_uncached(allele, trio_posteriors, info);
    }
}

auto compute_denovo_posteriors(const AllelePosteriorMap& called_alleles,
                               const TrioProbabilityVector& trio_posteriors,
                               const TrioGenotypeInfo& info)
{
    AllelePosteriorMap result {};
    for (const auto& p : called_alleles) {
        result.emplace(p.first, compute_denovo_posterior(p.first, trio_posteriors, info));
    }
    return result;
}

struct CalledDenovo : public Mappable<CalledDenovo>
{
    Allele allele;
    Phred<double> allele_posterior, denovo_posterior;
    CalledDenovo(Allele allele, Phred<double> allele_posterior, Phred<double> denovo_posterior)
    : allele {std::move(allele)}
    , allele_posterior {allele_posterior}
    , denovo_posterior {denovo_posterior}
    {}
    const GenomicRegion& mapped_region() const noexcept { return allele.mapped_region(); }
};

auto call_denovos(const AllelePosteriorMap& denovo_posteriors, const AllelePosteriorMap& segregating_posteriors,
                  const Phred<double> min_denovo_posterior)
{
    std::vector<CalledDenovo> result {};
    result.reserve(denovo_posteriors.size());
    for (const auto& p : denovo_posteriors) {
        if (p.second >= min_denovo_posterior) {
            result.emplace_back(p.first, segregating_posteriors.at(p.first), p.second);
        }
    }
    return result;
}

struct CallCompare
{
    bool operator()(const AllelePosteriorMap::value_type& lhs, const CalledDenovo& rhs) const
    {
        return lhs.first < rhs.allele;
    }
    bool operator()(const CalledDenovo& lhs, const AllelePosteriorMap::value_type& rhs) const
    {
        return lhs.allele < rhs.first;
    }
};

auto get_germline_alleles(const AllelePosteriorMap& called_alleles,
                          const std::vector<CalledDenovo>& denovos)
{
    std::vector<AllelePosteriorMap::value_type> result {};
    result.reserve(called_alleles.size() - denovos.size());
    std::set_difference(std::cbegin(called_alleles), std::cend(called_alleles),
                        std::cbegin(denovos), std::cend(denovos),
                        std::back_inserter(result), CallCompare {});
    return result;
}

struct CalledGermlineVariant : public Mappable<CalledGermlineVariant>
{
    Variant variant;
    Phred<double> posterior;
    CalledGermlineVariant(Variant variant, Phred<double> posterior)
    : variant {std::move(variant)}
    , posterior {posterior}
    {}
    const GenomicRegion& mapped_region() const noexcept { return variant.mapped_region(); }
};

boost::optional<Variant> find_variant(const Allele& allele, const std::vector<Variant>& variants)
{
    const auto er = std::equal_range(std::cbegin(variants), std::cend(variants), allele,
                                     [] (const auto& lhs, const auto& rhs) { return mapped_region(lhs) < mapped_region(rhs); });
    const auto itr = std::find_if(er.first, er.second, [&allele] (const Variant& v) { return v.alt_allele() == allele; });
    if (itr != er.second) {
        return *itr;
    } else {
        return boost::none;
    }
}

auto call_germline_variants(const std::vector<AllelePosteriorMap::value_type>& germline_allele_posteriors,
                            const std::vector<Variant>& variants, const Phred<double> min_posterior)
{
    std::vector<CalledGermlineVariant> result {};
    result.reserve(germline_allele_posteriors.size());
    for (const auto& p : germline_allele_posteriors) {
        if (p.second >= min_posterior) {
            const auto variant = find_variant(p.first, variants);
            if (variant) result.emplace_back(*variant, p.second);
        }
    }
    return result;
}

struct TrioCall
{
    Genotype<IndexedHaplotype<>> mother, father, child;
};

bool includes(const TrioCall& trio, const Allele& allele)
{
    return includes(trio.mother, allele) || includes(trio.father, allele) || includes(trio.child, allele);
}

bool none_mendilian_errors(const JointProbability& call, const std::vector<CalledGermlineVariant>& germline_calls, const TrioGenotypeInfo& info)
{
    return std::none_of(std::cbegin(germline_calls), std::cend(germline_calls),
                        [&] (const auto& germline) { return is_denovo(germline.variant.alt_allele(), call, info); });
}

bool all_mendilian_errors(const JointProbability& call, const std::vector<CalledDenovo>& denovo_calls, const TrioGenotypeInfo& info)
{
    return std::all_of(std::cbegin(denovo_calls), std::cend(denovo_calls),
                       [&] (const auto& denovo) { return is_denovo(denovo.allele, call, info); });
}

bool is_viable_genotype_call(const JointProbability& call,
                             const std::vector<CalledGermlineVariant>& germline_calls,
                             const std::vector<CalledDenovo>& denovo_calls,
                             const TrioGenotypeInfo& info)
{
    return none_mendilian_errors(call, germline_calls, info) && all_mendilian_errors(call, denovo_calls, info);
}

TrioCall to_call(const JointProbability& p, const TrioGenotypeInfo& info) noexcept
{
    return TrioCall {info.genotypes[p.maternal], info.genotypes[p.paternal], info.genotypes[p.child]};
}

auto call_trio(const TrioProbabilityVector& trio_posteriors,
               const std::vector<CalledGermlineVariant>& germline_calls,
               const std::vector<CalledDenovo>& denovo_calls,
               const TrioGenotypeInfo& info)
{
    assert(!trio_posteriors.empty());
    const auto map_itr = std::max_element(std::cbegin(trio_posteriors), std::cend(trio_posteriors),
                                          [] (const auto& lhs, const auto& rhs) {
                                              return lhs.probability < rhs.probability;
                                          });
    if (trio_posteriors.size() == 1 || is_viable_genotype_call(*map_itr, germline_calls, denovo_calls, info)) {
        return to_call(*map_itr, info);
    } else {
        std::vector<std::reference_wrapper<const JointProbability>> trio_posterior_refs {};
        trio_posterior_refs.reserve(trio_posteriors.size());
        std::copy(std::cbegin(trio_posteriors), std::cend(trio_posteriors), std::back_inserter(trio_posterior_refs));
        std::sort(std::begin(trio_posterior_refs), std::end(trio_posterior_refs),
                  [] (const auto& lhs, const auto& rhs) { return lhs.get().probability > rhs.get().probability; });
        auto viable_map_itr = std::find_if(std::next(std::cbegin(trio_posterior_refs)), std::cend(trio_posterior_refs),
                                           [&] (const auto& p) {
                                               return is_viable_genotype_call(p,  germline_calls, denovo_calls, info);
                                           });
        if (viable_map_itr != std::cend(trio_posterior_refs)) {
            return to_call(*viable_map_itr, info);
        } else {
            return to_call(*map_itr, info);
        }
    }
}

bool includes(const TrioCall& trio, const CalledGermlineVariant& call)
{
    return includes(trio, call.variant.alt_allele());
}

bool includes(const TrioCall& trio, const CalledDenovo& call)
{
    return includes(trio, call.allele);
}

template <typename T>
void remove_ungenotyped_allele(std::vector<T>& calls, const TrioCall& trio)
{
    auto itr = std::remove_if(std::begin(calls), std::end(calls),
                              [&] (const auto& call) { return !includes(trio, call); });
    calls.erase(itr, std::end(calls));
}

void remove_ungenotyped_allele(std::vector<CalledGermlineVariant>& germline_calls,
                               std::vector<CalledDenovo>& denovo_calls,
                               const TrioCall& trio)
{
    remove_ungenotyped_allele(germline_calls, trio);
    remove_ungenotyped_allele(denovo_calls, trio);
}

using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<IndexedHaplotype<>>>;

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

struct GenotypedTrio
{
    GenotypePosterior mother, father, child;
};

auto call_genotypes(const Trio& trio, const TrioCall& called_trio,
                    const GenotypeProbabilityMap& trio_posteriors,
                    const std::vector<GenomicRegion>& regions)
{
    std::vector<GenotypedTrio> result {};
    result.reserve(regions.size());
    for (const auto& region : regions) {
        auto mother_genotype = copy<Allele>(called_trio.mother, region);
        auto mother_posterior = compute_posterior(mother_genotype, trio_posteriors[trio.mother()]);
        auto father_genotype = copy<Allele>(called_trio.father, region);
        auto father_posterior = compute_posterior(father_genotype, trio_posteriors[trio.father()]);
        auto child_genotype = copy<Allele>(called_trio.child, region);
        auto child_posterior = compute_posterior(child_genotype, trio_posteriors[trio.child()]);
        result.push_back({{std::move(mother_genotype), mother_posterior},
                          {std::move(father_genotype), father_posterior},
                          {std::move(child_genotype), child_posterior}});
    }
    return result;
}

bool is_reference_reversion(const Allele& denovo, const std::map<Allele, Allele>& reference_alleles)
{
    return reference_alleles.at(denovo) == denovo;
}

auto make_variant(Allele&& denovo, const std::map<Allele, Allele>& reference_alleles)
{
    return Variant {reference_alleles.at(denovo), std::move(denovo)};
}

auto make_genotype_calls(GenotypedTrio&& call, const Trio& trio)
{
    return std::vector<std::pair<SampleName, Call::GenotypeCall>> {
        {trio.mother(), {std::move(call.mother.genotype), call.mother.posterior}},
        {trio.father(), {std::move(call.father.genotype), call.father.posterior}},
        {trio.child(), {std::move(call.child.genotype), call.child.posterior}}
    };
}

auto make_calls(std::vector<CalledDenovo>&& alleles,
                std::vector<GenotypedTrio>&& genotypes,
                const Trio& trio,
                const std::vector<Variant>& candidates,
                const boost::optional<Phred<double>> max_quality = boost::none)
{
    std::map<Allele, Allele> reference_alleles {};
    for (const auto& denovo : alleles) {
        auto itr = std::find_if(std::cbegin(candidates), std::cend(candidates),
                                [denovo] (const auto& c) { return is_same_region(c, denovo); });
        reference_alleles.emplace(denovo.allele, itr->ref_allele());
    }
    std::vector<std::unique_ptr<VariantCall>> result {};
    result.reserve(alleles.size());
    std::transform(std::make_move_iterator(std::begin(alleles)), std::make_move_iterator(std::end(alleles)),
                   std::make_move_iterator(std::begin(genotypes)), std::back_inserter(result),
                   [&trio, &reference_alleles, max_quality] (auto&& denovo, auto&& genotype) -> std::unique_ptr<DenovoCall> {
                       if (max_quality) denovo.allele_posterior = std::min(denovo.allele_posterior, *max_quality);
                       if (is_reference_reversion(denovo.allele, reference_alleles)) {
                           return std::make_unique<DenovoReferenceReversionCall>(std::move(denovo.allele),
                                                                                 make_genotype_calls(std::move(genotype), trio),
                                                                                 denovo.allele_posterior, denovo.denovo_posterior);
                       } else {
                           return std::make_unique<DenovoCall>(make_variant(std::move(denovo.allele), reference_alleles),
                                                               make_genotype_calls(std::move(genotype), trio),
                                                               denovo.allele_posterior, denovo.denovo_posterior);
                       }
                   });
    return result;
}

auto make_calls(std::vector<CalledGermlineVariant>&& variants,
                std::vector<GenotypedTrio>&& genotypes,
                const Trio& trio,
                const boost::optional<Phred<double>> max_quality = boost::none)
{
    std::vector<std::unique_ptr<VariantCall>> result {};
    result.reserve(variants.size());
    std::transform(std::make_move_iterator(std::begin(variants)), std::make_move_iterator(std::end(variants)),
                   std::make_move_iterator(std::begin(genotypes)), std::back_inserter(result),
                   [&trio, max_quality] (auto&& variant, auto&& genotype) {
                       if (max_quality) variant.posterior = std::min(variant.posterior, *max_quality);
                       return std::make_unique<GermlineVariantCall>(std::move(variant.variant),
                                                                    make_genotype_calls(std::move(genotype), trio),
                                                                    variant.posterior);
                   });
    return result;
}

auto make_calls(std::vector<CalledGermlineVariant>&& variants,
                std::vector<GenotypedTrio>&& germline_genotypes,
                std::vector<CalledDenovo>&& alleles,
                std::vector<GenotypedTrio>&& denovo_genotypes,
                const Trio& trio,
                const std::vector<Variant>& candidates,
                const boost::optional<Phred<double>> max_quality = boost::none)
{
    auto germline_calls = make_calls(std::move(variants), std::move(germline_genotypes), trio, max_quality);
    auto denovo_calls = make_calls(std::move(alleles), std::move(denovo_genotypes), trio, candidates, max_quality);
    std::vector<std::unique_ptr<VariantCall>> result {};
    result.reserve(germline_calls.size() + denovo_calls.size());
    std::merge(std::make_move_iterator(std::begin(germline_calls)), std::make_move_iterator(std::end(germline_calls)),
               std::make_move_iterator(std::begin(denovo_calls)), std::make_move_iterator(std::end(denovo_calls)),
               std::back_inserter(result), [] (const auto& lhs, const auto& rhs) {
                   return lhs->mapped_region() < rhs->mapped_region();
               });
    return result;
}

} // namespace

namespace debug {

template <typename Container>
void log(const Container& maternal_genotypes,
         const boost::optional<Container>& paternal_genotypes,
         const TrioProbabilityVector& posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log);
void log(const AllelePosteriorMap& posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log,
         Phred<double> min_posterior, bool denovo = false);

} // namespace debug

std::vector<std::unique_ptr<VariantCall>>
TrioCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    const auto alleles = decompose(candidates);
    const auto& trio_posteriors = latents.model_latents.posteriors.joint_genotype_probabilities;
    if (latents.concatenated_genotypes_.empty()) {
        debug::log(latents.maternal_genotypes, latents.paternal_genotypes, trio_posteriors, debug_log_, trace_log_);
    }
    const TrioGenotypeInfo info {latents.concatenated_genotypes_.empty() ? latents.maternal_genotypes : latents.concatenated_genotypes_, latents.haplotypes};
    const auto allele_posteriors = compute_segregation_posteriors(alleles, trio_posteriors, info);
    debug::log(allele_posteriors, debug_log_, trace_log_, parameters_.min_variant_posterior);
    const auto called_alleles = call_alleles(allele_posteriors, parameters_.min_variant_posterior);
    const auto denovo_posteriors = compute_denovo_posteriors(called_alleles, trio_posteriors, info);
    debug::log(denovo_posteriors, debug_log_, trace_log_, parameters_.min_denovo_posterior, true);
    auto denovos = call_denovos(denovo_posteriors, allele_posteriors, parameters_.min_denovo_posterior);
    const auto germline_alleles = get_germline_alleles(called_alleles, denovos);
    auto germline_variants = call_germline_variants(germline_alleles, candidates, parameters_.min_variant_posterior);
    auto called_trio = call_trio(trio_posteriors, germline_variants, denovos, info);
    remove_ungenotyped_allele(germline_variants, denovos, called_trio);
    if (parameters_.maternal_ploidy == 0) called_trio.mother = Genotype<IndexedHaplotype<>> {};
    if (parameters_.paternal_ploidy == 0) called_trio.father = Genotype<IndexedHaplotype<>> {};
    if (parameters_.child_ploidy == 0) called_trio.child = Genotype<IndexedHaplotype<>> {};
    auto denovo_genotypes = call_genotypes(parameters_.trio, called_trio, *latents.genotype_posteriors(), extract_regions(denovos));
    auto germline_genotypes = call_genotypes(parameters_.trio, called_trio, *latents.genotype_posteriors(), extract_regions(germline_variants));
    boost::optional<Phred<double>> max_quality {};
    if (latents.model_latents.estimated_lost_log_posterior_mass) {
        max_quality = log_probability_false_to_phred(*latents.model_latents.estimated_lost_log_posterior_mass);
    }
    return make_calls(std::move(germline_variants), std::move(germline_genotypes),
                      std::move(denovos), std::move(denovo_genotypes),
                      parameters_.trio, candidates,
                      max_quality);
}

std::vector<std::unique_ptr<ReferenceCall>>
TrioCaller::call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents,
                           const ReadPileupMap& pileups) const
{
    return call_reference(alleles, dynamic_cast<const Latents&>(latents), pileups);
}

std::vector<std::unique_ptr<ReferenceCall>>
TrioCaller::call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                           const ReadPileupMap& pileups) const
{
    return {};
}

std::unique_ptr<PopulationPriorModel> TrioCaller::make_prior_model(const HaplotypeBlock& haplotypes) const
{
    if (parameters_.germline_prior_model_params) {
        return std::make_unique<CoalescentPopulationPriorModel>(CoalescentModel {Haplotype {mapped_region(haplotypes), reference_},
                                                                                 *parameters_.germline_prior_model_params,
                                                                                 haplotypes.size(), CoalescentModel::CachingStrategy::address});
    } else {
        return std::make_unique<UniformPopulationPriorModel>();
    }
}

std::unique_ptr<GenotypePriorModel> TrioCaller::make_single_sample_prior_model(const HaplotypeBlock& haplotypes) const
{
    if (parameters_.germline_prior_model_params) {
        return std::make_unique<CoalescentGenotypePriorModel>(CoalescentModel {Haplotype {mapped_region(haplotypes), reference_},
                                                                               *parameters_.germline_prior_model_params,
                                                                               haplotypes.size(), CoalescentModel::CachingStrategy::address});
    } else {
        return std::make_unique<UniformGenotypePriorModel>();
    }
}

namespace {

template <typename IndexType>
void erase_duplicates(MappableBlock<Genotype<IndexedHaplotype<IndexType>>>& genotypes)
{
    using std::begin; using std::end;
    std::sort(begin(genotypes), end(genotypes), GenotypeLess {});
    genotypes.erase(std::unique(begin(genotypes), end(genotypes)), end(genotypes));
}

} // namespace

std::pair<TrioCaller::GenotypeBlock, TrioCaller::GenotypeBlock>
TrioCaller::propose_model_check_genotypes(const HaplotypeBlock& haplotypes,
                                          const IndexedHaplotypeBlock& indexed_haplotypes,
                                          const GenotypeBlock& genotypes,
                                          const std::vector<double>& genotype_posteriors) const
{
    constexpr std::size_t max_model_check_genotypes {5};
    const auto num_model_check_genotypes = std::min(max_model_check_genotypes, genotypes.size());
    const auto best_indices = select_top_k_indices(genotype_posteriors, num_model_check_genotypes, false);
    GenotypeBlock assumed {mapped_region(haplotypes)};
    for (const auto genotype_idx : best_indices) {
        assumed.push_back(genotypes[genotype_idx]);
    }
    auto augmented = extend(assumed, indexed_haplotypes);
    erase_duplicates(augmented);
    return {std::move(assumed), std::move(augmented)};
}

namespace debug {

template <typename S, typename Container>
void print(S&& stream, 
           const Container& maternal_genotypes,
           const Container& paternal_genotypes,
           const Container& child_genotypes,
           const TrioProbabilityVector& posteriors,
           const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    const auto m = std::min(n, posteriors.size());
    if (m == posteriors.size()) {
        stream << "Printing all trio joint genotype posteriors (maternal | paternal | child)" << '\n';
    } else {
        stream << "Printing top " << m << " trio joint genotype posteriors (maternal | paternal | child)" << '\n';
    }
    std::vector<JointProbability> v {};
    v.reserve(posteriors.size());
    std::copy(std::cbegin(posteriors), std::cend(posteriors), std::back_inserter(v));
    const auto mth = std::next(std::begin(v), m);
    std::partial_sort(std::begin(v), mth, std::end(v),
                      [] (const auto& lhs, const auto& rhs) { return lhs.probability > rhs.probability; });
    std::for_each(std::begin(v), mth,
                  [&] (const auto& p) {
                      using octopus::debug::print_variant_alleles;
                      print_variant_alleles(stream, maternal_genotypes[p.maternal]);
                      stream << " | ";
                      print_variant_alleles(stream, paternal_genotypes[p.paternal]);
                      stream << " | ";
                      print_variant_alleles(stream, child_genotypes[p.child]);
                      stream << " " << p.probability << "\n";
                  });
}

template <typename Container>
void log(const Container& maternal_genotypes,
         const Container& paternal_genotypes,
         const Container& child_genotypes,
         const TrioProbabilityVector& posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log)
{
    if (trace_log) {
        print(stream(*trace_log), maternal_genotypes, paternal_genotypes, child_genotypes, posteriors);
    }
    if (debug_log) {
        print(stream(*debug_log), maternal_genotypes, paternal_genotypes, child_genotypes, posteriors, 10);
    }
}

template <typename Container>
void log(const Container& maternal_genotypes,
         const boost::optional<Container>& paternal_genotypes,
         const TrioProbabilityVector& posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log)
{
    if (paternal_genotypes) {
        log(maternal_genotypes, *paternal_genotypes, maternal_genotypes, posteriors, debug_log, trace_log);
    } else {
        log(maternal_genotypes, maternal_genotypes, maternal_genotypes, posteriors, debug_log, trace_log);
    }
}

template <typename S>
void print(S&& stream, const AllelePosteriorMap& posteriors,
           const std::string& type = "allele",
           const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    const auto m = std::min(n, posteriors.size());
    if (m == posteriors.size()) {
        stream << "Printing all " << type << " posteriors" << '\n';
    } else {
        stream << "Printing top " << m << " " << type << " posteriors" << '\n';
    }
    std::vector<std::pair<Allele, Phred<double>>> v {};
    v.reserve(posteriors.size());
    std::copy(std::cbegin(posteriors), std::cend(posteriors), std::back_inserter(v));
    const auto mth = std::next(std::begin(v), m);
    std::partial_sort(std::begin(v), mth, std::end(v),
                      [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
    std::for_each(std::begin(v), mth,
                  [&] (const auto& p) {
                      stream << p.first << " " << p.second.probability_true() << '\n';
                  });
}

void log(const AllelePosteriorMap& posteriors,
         boost::optional<logging::DebugLogger>& debug_log,
         boost::optional<logging::TraceLogger>& trace_log,
         Phred<double> min_posterior, const bool denovo)
{
    if (!denovo || !posteriors.empty()) {
        const std::string type {denovo ? "denovo allele" : "allele"};
        if (trace_log) {
            print(stream(*trace_log), posteriors, type);
        }
        if (debug_log) {
            const auto n = std::count_if(std::cbegin(posteriors), std::cend(posteriors),
                                         [=] (const auto& p) { return p.second >= min_posterior; });
            print(stream(*debug_log), posteriors, type, std::max(n, decltype(n) {10}));
        }
    }
}
    
} // namespace debug

} // namespace octopus
