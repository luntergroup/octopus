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
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "utils/maths.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"
#include "containers/probability_matrix.hpp"
#include "logging/logging.hpp"
#include "core/types/calls/germline_variant_call.hpp"
#include "core/types/calls/reference_call.hpp"

#include "core/models/genotype/uniform_genotype_prior_model.hpp"
#include "core/models/genotype/coalescent_genotype_prior_model.hpp"

#include "timers.hpp"

namespace octopus {

ProkaryoteCaller::ProkaryoteCaller(Caller::Components&& components,
                                   Caller::Parameters general_parameters,
                                   Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {std::move(specific_parameters)}
{
    if (parameters_.max_clones == 0) {
        throw std::logic_error {"ProkaryoteCaller: max_clones must be > 0"};
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

std::shared_ptr<ProkaryoteCaller::Latents::HaplotypeProbabilityMap>
ProkaryoteCaller::Latents::haplotype_posteriors() const noexcept
{
    return haplotype_posteriors_;
}

std::shared_ptr<ProkaryoteCaller::Latents::GenotypeProbabilityMap>
ProkaryoteCaller::Latents::genotype_posteriors() const noexcept
{
    return genotype_posteriors_;
}

// ProkaryoteCaller::Latents private methods

std::unique_ptr<ProkaryoteCaller::Caller::Latents>
ProkaryoteCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                                const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    return nullptr;
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
