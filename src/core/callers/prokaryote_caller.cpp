// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "prokaryote_caller.hpp"

#include <typeinfo>
#include <unordered_map>
#include <deque>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <utility>
#include <stdexcept>
#include <iostream>

#include "basics/genomic_region.hpp"
#include "containers/probability_matrix.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/calls/germline_variant_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "core/models/genotype/uniform_genotype_prior_model.hpp"
#include "core/models/genotype/coalescent_genotype_prior_model.hpp"
#include "utils/maths.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "utils/concat.hpp"
#include "logging/logging.hpp"

namespace octopus {

ProkaryoteCaller::ProkaryoteCaller(Caller::Components&& components,
                                   Caller::Parameters general_parameters,
                                   Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {std::move(specific_parameters)}
{
    if (parameters_.max_clones < 1) {
        throw std::logic_error {"ProkaryoteCaller: max_clones must be > 1"};
    }
}

std::string ProkaryoteCaller::do_name() const
{
    return "prokaryote";
}

ProkaryoteCaller::CallTypeSet ProkaryoteCaller::do_call_types() const
{
    return {std::type_index(typeid(GermlineVariantCall))};
}

std::size_t ProkaryoteCaller::do_remove_duplicates(std::vector<Haplotype>& haplotypes) const
{
    if (parameters_.deduplicate_haplotypes_with_germline_model) {
        if (haplotypes.size() < 2) return 0;
        CoalescentModel::Parameters model_params {};
        if (parameters_.prior_model_params) model_params = *parameters_.prior_model_params;
        Haplotype reference {mapped_region(haplotypes.front()), reference_.get()};
        CoalescentModel model {std::move(reference), model_params, haplotypes.size(), CoalescentModel::CachingStrategy::none};
        const CoalescentProbabilityGreater cmp {std::move(model)};
        return octopus::remove_duplicates(haplotypes, cmp);
    } else {
        return Caller::do_remove_duplicates(haplotypes);
    }
}

// ProkaryoteCaller::Latents public methods

ProkaryoteCaller::Latents::Latents(std::vector<Genotype<Haplotype>> haploid_genotypes, std::vector<Genotype<Haplotype>> polyploid_genotypes,
                                   HaploidModelInferences haploid_model_inferences, SubloneModelInferences subclone_model_inferences,
                                   const SampleName& sample)
: haploid_genotypes_ {std::move(haploid_genotypes)}
, polyploid_genotypes_ {std::move(polyploid_genotypes)}
, haploid_model_inferences_ {std::move(haploid_model_inferences)}
, subclone_model_inferences_ {std::move(subclone_model_inferences)}
, model_posteriors_ {}
, sample_ {sample}
{
    // Bayes Factor
    const auto norm = maths::log_sum_exp({haploid_model_inferences_.log_evidence, subclone_model_inferences_.approx_log_evidence});
    model_posteriors_.clonal = std::exp(haploid_model_inferences_.log_evidence - norm);
    model_posteriors_.subclonal = std::exp(subclone_model_inferences_.approx_log_evidence - norm);
}

std::shared_ptr<ProkaryoteCaller::Latents::HaplotypeProbabilityMap>
ProkaryoteCaller::Latents::haplotype_posteriors() const noexcept
{
    if (haplotype_posteriors_ == nullptr) {
        haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>();
        for (const auto& p : (*(this->genotype_posteriors()))[sample_]) {
            for (const auto& haplotype : p.first.copy_unique_ref()) {
                (*haplotype_posteriors_)[haplotype] += p.second;
            }
        }
    }
    return haplotype_posteriors_;
}

std::shared_ptr<ProkaryoteCaller::Latents::GenotypeProbabilityMap>
ProkaryoteCaller::Latents::genotype_posteriors() const noexcept
{
    if (genotype_posteriors_ == nullptr) {
        const auto genotypes = concat(haploid_genotypes_, polyploid_genotypes_);
        auto posteriors = concat(haploid_model_inferences_.posteriors.genotype_probabilities,
                                 subclone_model_inferences_.posteriors.genotype_probabilities);
        std::for_each(std::begin(posteriors), std::next(std::begin(posteriors), haploid_genotypes_.size()),
                      [this] (auto& p) { p *= model_posteriors_.clonal; });
        std::for_each(std::next(std::begin(posteriors), haploid_genotypes_.size()), std::end(posteriors),
                      [this] (auto& p) { p *= model_posteriors_.subclonal; });
        genotype_posteriors_ = std::make_shared<GenotypeProbabilityMap>(std::make_move_iterator(std::begin(genotypes)),
                                                                        std::make_move_iterator(std::end(genotypes)));
        insert_sample(sample_, posteriors, *genotype_posteriors_);
    }
    return genotype_posteriors_;
}

// ProkaryoteCaller::Latents private methods

namespace {

auto make_sublone_model_mixture_prior_map(const SampleName& sample, const unsigned num_clones, const double alpha = 0.5)
{
    model::SubcloneModel::Priors::GenotypeMixturesDirichletAlphaMap result {};
    model::SubcloneModel::Priors::GenotypeMixturesDirichletAlphas alphas(num_clones, alpha);
    result.emplace(sample, std::move(alphas));
    return result;
}

auto fit_sublone_model(const std::vector<Haplotype>& haplotypes, const HaplotypeLikelihoodCache& haplotype_likelihoods,
                       const GenotypePriorModel& genotype_prior_model, const SampleName& sample, const unsigned max_clones,
                       const double haploid_model_evidence)
{
    std::vector<Genotype<Haplotype>> polyploid_genotypes {};
    model::SubcloneModel::InferredLatents sublonal_inferences;
    for (unsigned num_clones {2}; num_clones <= max_clones; ++num_clones) {
        auto genotypes = generate_all_full_rank_genotypes(haplotypes, num_clones);
        if (genotypes.empty()) break;
        model::SubcloneModel::Priors subclonal_model_priors {genotype_prior_model, make_sublone_model_mixture_prior_map(sample, num_clones)};
        model::SubcloneModel subclonal_model {{sample}, subclonal_model_priors};
        auto inferences = subclonal_model.evaluate(genotypes, haplotype_likelihoods);
        if (num_clones == 2) {
            polyploid_genotypes = std::move(genotypes);
            sublonal_inferences = std::move(inferences);
            if (sublonal_inferences.approx_log_evidence < haploid_model_evidence) break;
        } else {
            if (inferences.approx_log_evidence <= sublonal_inferences.approx_log_evidence)  break;
            polyploid_genotypes = std::move(genotypes);
            sublonal_inferences = std::move(inferences);
        }
    }
    return std::make_pair(std::move(polyploid_genotypes), std::move(sublonal_inferences));
}

} // namespace

std::unique_ptr<ProkaryoteCaller::Caller::Latents>
ProkaryoteCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                                const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    auto haploid_genotypes = generate_all_genotypes(haplotypes, 1);
    if (debug_log_) stream(*debug_log_) << "There are " << haploid_genotypes.size() << " candidate haploid genotypes";
    const auto genotype_prior_model = make_prior_model(haplotypes);
    const model::IndividualModel haploid_model {*genotype_prior_model, debug_log_};
    haplotype_likelihoods.prime(sample());
    auto haploid_inferences = haploid_model.evaluate(haploid_genotypes, haplotype_likelihoods);
    std::vector<Genotype<Haplotype>> polyploid_genotypes; model::SubcloneModel::InferredLatents sublonal_inferences;
    std::tie(polyploid_genotypes, sublonal_inferences) = fit_sublone_model(haplotypes, haplotype_likelihoods, *genotype_prior_model,
                                                                           sample(), parameters_.max_clones, haploid_inferences.log_evidence);
    if (debug_log_) stream(*debug_log_) << "There are " << polyploid_genotypes.size() << " candidate polyploid genotypes";
    using std::move;
    return std::make_unique<Latents>(move(haploid_genotypes), move(polyploid_genotypes),
                                     move(haploid_inferences), move(sublonal_inferences),
                                     sample());
}

boost::optional<double>
ProkaryoteCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                            const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                            const Caller::Latents& latents) const
{
    return calculate_model_posterior(haplotypes, haplotype_likelihoods, dynamic_cast<const Latents&>(latents));
}

boost::optional<double>
ProkaryoteCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                            const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                            const Latents& latents) const
{
    return boost::none;
}

std::vector<std::unique_ptr<octopus::VariantCall>>
ProkaryoteCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

std::vector<std::unique_ptr<octopus::VariantCall>>
ProkaryoteCaller::call_variants(const std::vector<Variant>& candidates,
                                const Latents& latents) const
{
    return {};
}

std::vector<std::unique_ptr<ReferenceCall>>
ProkaryoteCaller::call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents, const ReadMap& reads) const
{
    return call_reference(alleles, dynamic_cast<const Latents&>(latents), reads);
}

std::vector<std::unique_ptr<ReferenceCall>>
ProkaryoteCaller::call_reference(const std::vector<Allele>& alleles, const Latents& latents, const ReadMap& reads) const
{
    return {};
}

const SampleName& ProkaryoteCaller::sample() const noexcept
{
    return samples_.front();
}

std::unique_ptr<GenotypePriorModel> ProkaryoteCaller::make_prior_model(const std::vector<Haplotype>& haplotypes) const
{
    if (parameters_.prior_model_params) {
        return std::make_unique<CoalescentGenotypePriorModel>(CoalescentModel {
        Haplotype {mapped_region(haplotypes.front()), reference_},
        *parameters_.prior_model_params, haplotypes.size(), CoalescentModel::CachingStrategy::address
        });
    } else {
        return std::make_unique<UniformGenotypePriorModel>();
    }
}

} // namespace octopus
