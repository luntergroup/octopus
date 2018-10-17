// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cell_caller.hpp"

#include <typeinfo>
#include <unordered_map>
#include <deque>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <utility>
#include <stdexcept>
#include <iostream>

#include "basics/genomic_region.hpp"
#include "containers/probability_matrix.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/calls/germline_variant_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "logging/logging.hpp"

namespace octopus {

CellCaller::CellCaller(Caller::Components&& components,
                       Caller::Parameters general_parameters,
                       Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_{std::move(specific_parameters)}
{}

std::string CellCaller::do_name() const
{
    return "cell";
}

CellCaller::CallTypeSet CellCaller::do_call_types() const
{
    return {std::type_index(typeid(GermlineVariantCall))};
}

unsigned CellCaller::do_min_callable_ploidy() const
{
    return parameters_.ploidy;
}

unsigned CellCaller::do_max_callable_ploidy() const
{
    return parameters_.ploidy;
}

std::size_t CellCaller::do_remove_duplicates(std::vector<Haplotype>& haplotypes) const
{
    if (parameters_.deduplicate_haplotypes_with_prior_model) {
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

// CellCaller::Latents public methods

std::shared_ptr<CellCaller::Latents::HaplotypeProbabilityMap>
CellCaller::Latents::haplotype_posteriors() const noexcept
{
    if (haplotype_posteriors_ == nullptr) {
        // TODO
    }
    return haplotype_posteriors_;
}

std::shared_ptr<CellCaller::Latents::GenotypeProbabilityMap>
CellCaller::Latents::genotype_posteriors() const noexcept
{
    if (genotype_posteriors_ == nullptr) {
        // TODO
    }
    return genotype_posteriors_;
}

// CellCaller::Latents private methods

std::unique_ptr<CellCaller::Caller::Latents>
CellCaller::infer_latents(const std::vector<Haplotype>& haplotypes, const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    return nullptr;
}

boost::optional<double>
CellCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                      const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                      const Caller::Latents& latents) const
{
    return calculate_model_posterior(haplotypes, haplotype_likelihoods, dynamic_cast<const Latents&>(latents));
}

boost::optional<double>
CellCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                      const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                      const Latents& latents) const
{
    return boost::none;
}

std::vector<std::unique_ptr<octopus::VariantCall>>
CellCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

std::vector<std::unique_ptr<octopus::VariantCall>>
CellCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    return {};
}

std::vector<std::unique_ptr<ReferenceCall>>
CellCaller::call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents, const ReadPileupMap& pileup) const
{
    return call_reference(alleles, dynamic_cast<const Latents&>(latents), pileup);
}

std::vector<std::unique_ptr<ReferenceCall>>
CellCaller::call_reference(const std::vector<Allele>& alleles, const Latents& latents, const ReadPileupMap& pileup) const
{
    return {};
}

} // namespace octopus
