// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "trio_caller.hpp"

#include <typeinfo>
#include <iterator>
#include <algorithm>
#include <numeric>

#include "basics/genomic_region.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/genotype/trio_model.hpp"
#include "utils/germline_variant_call.hpp"
#include "utils/reference_call.hpp"

namespace octopus {

TrioCaller::TrioCaller(Caller::Components&& components,
                       Caller::Parameters general_parameters,
                       Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {std::move(specific_parameters)}
{
    if (parameters_.maternal_ploidy == 0) {
        throw std::logic_error {"IndividualCaller: ploidy must be > 0"};
    }
}

Caller::CallTypeSet TrioCaller::do_get_call_types() const
{
    return {std::type_index(typeid(GermlineVariantCall))};
}

namespace {

using JointProbability = model::TrioModel::Latents::JointProbability;

struct GenotypeProbabilityPair
{
    Genotype<Haplotype> genotype;
    double probability;
};

using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;

bool operator==(const GenotypeReference lhs, const GenotypeReference rhs)
{
    return lhs.get() == rhs.get();
}

bool operator<(const GenotypeReference lhs, const GenotypeReference rhs)
{
    return lhs.get() < rhs.get();
}

template <typename Function>
auto marginalise(std::vector<JointProbability>& joint_posteriors, Function who)
{
    std::sort(std::begin(joint_posteriors), std::end(joint_posteriors),
              [who](const auto& lhs, const auto& rhs) {
                  return who(lhs) < who(rhs);
              });
    std::vector<GenotypeProbabilityPair> result {};
    result.reserve(joint_posteriors.size());
    for (auto iter = std::cbegin(joint_posteriors), end = std::cend(joint_posteriors); iter != end;) {
        const auto next = std::find_if_not(std::next(iter), end,
                                           [iter, who](const auto& p) {
                                               return who(p) == who(*iter);
                                           });
        result.push_back({who(*iter),
                          std::accumulate(iter, next, 0.0,
                                          [](const auto curr, const auto& p) {
                                              return curr + p.probability;
                                          })});
        iter = next;
    }
    result.shrink_to_fit();
    return result;
}

auto marginalise_mother(std::vector<JointProbability>& joint_posteriors)
{
    return marginalise(joint_posteriors, [](const JointProbability& p) { return p.maternal; });
}

auto marginalise_father(std::vector<JointProbability>& joint_posteriors)
{
    return marginalise(joint_posteriors, [](const JointProbability& p) { return p.paternal; });
}

auto marginalise_child(std::vector<JointProbability>& joint_posteriors)
{
    return marginalise(joint_posteriors, [](const JointProbability& p) { return p.child; });
}

//void fill_missing_genotypes(std::vector<GenotypeProbabilityPair>& posteriors,
//                            std::vector<Genotype<Haplotype>>& genotypes)
//{
//    std::sort(std::begin(posteriors), std::end(posteriors),
//              [] (const auto& lhs, const auto& rhs) { return lhs.genotype < rhs.genotype; });
//    std::sort(std::begin(genotypes), std::end(genotypes));
//    std::vector<Genotype<Haplotype>> missing {};
//    missing.reserve(genotypes.size());
//    std::set_difference(std::cbegin(genotypes), std::cend(genotypes),
//                        std::cbegin(posteriors), std::cend(posteriors),
//                        std::back_inserter(missing),
//                        [] (const auto))
//}
    
} // namespace

std::unique_ptr<Caller::Latents>
TrioCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    const CoalescentModel germline_prior_model {
        Haplotype {mapped_region(haplotypes.front()), reference_},
        parameters_.germline_prior_model_params
    };
    const DeNovoModel denovo_model {parameters_.denovo_model_params};
    const model::TrioModel model {parameters_.trio, germline_prior_model, denovo_model};
    auto genotypes = generate_all_genotypes(haplotypes, parameters_.maternal_ploidy);
    auto latents = model.evaluate(genotypes, genotypes, genotypes, haplotype_likelihoods);
    return nullptr;
}

boost::optional<double>
TrioCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                      const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                      const Caller::Latents& latents) const
{
    return calculate_model_posterior(haplotypes, haplotype_likelihoods, dynamic_cast<const Latents&>(latents));
}

boost::optional<double>
TrioCaller::calculate_model_posterior(const std::vector<Haplotype>& haplotypes,
                                      const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                      const Latents& latents) const
{
    return boost::none;
}

std::vector<std::unique_ptr<VariantCall>>
TrioCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

std::vector<std::unique_ptr<VariantCall>>
TrioCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    return {};
}

std::vector<std::unique_ptr<ReferenceCall>>
TrioCaller::call_reference(const std::vector<Allele>& alleles, const Caller::Latents& latents,
                           const ReadMap& reads) const
{
    return call_reference(alleles, dynamic_cast<const Latents&>(latents), reads);
}

std::vector<std::unique_ptr<ReferenceCall>>
TrioCaller::call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                            const ReadMap& reads) const
{
    return {};
}

} // namespace octopus
