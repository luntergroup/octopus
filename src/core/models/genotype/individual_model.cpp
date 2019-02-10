// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "individual_model.hpp"

#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

#include "utils/maths.hpp"
#include "constant_mixture_genotype_likelihood_model.hpp"

namespace octopus { namespace model {

IndividualModel::IndividualModel(const GenotypePriorModel& genotype_prior_model,
                                 boost::optional<logging::DebugLogger> debug_log,
                                 boost::optional<logging::TraceLogger> trace_log)
: genotype_prior_model_ {genotype_prior_model}
, haplotypes_ {}
, debug_log_ {debug_log}
, trace_log_ {trace_log}
{}

const GenotypePriorModel& IndividualModel::prior_model() const noexcept
{
    return genotype_prior_model_;
}

void IndividualModel::prime(const std::vector<Haplotype>& haplotypes)
{
    haplotypes_ = std::addressof(haplotypes);
}

void IndividualModel::unprime() noexcept
{
    haplotypes_ = nullptr;
}

bool IndividualModel::is_primed() const noexcept
{
    return haplotypes_;
}

namespace debug {

template <typename S>
void print_genotype_priors(S&& stream, const std::vector<Genotype<Haplotype>>& genotypes,
                           const std::vector<double>& priors, std::size_t n = 5);
void print_genotype_priors(const std::vector<Genotype<Haplotype>>& genotypes,
                           const std::vector<double>& priors, std::size_t n = 5);
void log_genotype_likelihoods(boost::optional<logging::DebugLogger>& debug_log,
                              boost::optional<logging::TraceLogger>& trace_log,
                              const std::vector<Genotype<Haplotype>>& genotypes,
                              const std::vector<double>& likelihoods);
template <typename S>
void print_genotype_likelihoods(S&& stream, const std::vector<Genotype<Haplotype>>& genotypes,
                                const std::vector<double>& likelihoods, std::size_t n = 5);
void print_genotype_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                const std::vector<double>& likelihoods, std::size_t n = 5);

} // namespace debug

IndividualModel::InferredLatents
IndividualModel::evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    const ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods};
    auto posteriors = octopus::model::evaluate(genotypes, likelihood_model);
    debug::log_genotype_likelihoods(debug_log_, trace_log_, genotypes, posteriors);
    octopus::evaluate(genotypes, genotype_prior_model_, posteriors, false, true);
    const auto log_evidence = maths::normalise_exp(posteriors);
    return {{std::move(posteriors)}, log_evidence};
}

IndividualModel::InferredLatents
IndividualModel::evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
                          const std::vector<GenotypeIndex>& genotype_indices,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    assert(genotypes.size() == genotype_indices.size());
    ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods};
    InferredLatents result {};
    if (is_primed()) {
        likelihood_model.prime(*haplotypes_);
        result.posteriors.genotype_probabilities = octopus::model::evaluate(genotype_indices, likelihood_model);
    } else {
        result.posteriors.genotype_probabilities = octopus::model::evaluate(genotypes, likelihood_model);
    }
    debug::log_genotype_likelihoods(debug_log_, trace_log_, genotypes, result.posteriors.genotype_probabilities);
    octopus::evaluate(genotype_indices, genotype_prior_model_, result.posteriors.genotype_probabilities, false, true);
    result.log_evidence = maths::normalise_exp(result.posteriors.genotype_probabilities);
    return result;
}

namespace debug {

using octopus::debug::print_variant_alleles;

template <typename S>
void print_genotype_priors(S&& stream, const std::vector<Genotype<Haplotype>>& genotypes,
                           const std::vector<double>& priors, const std::size_t n)
{
    assert(genotypes.size() == priors.size());
    const auto m = std::min(n, genotypes.size());
    if (m == genotypes.size()) {
        stream << "Printing all genotype priors " << '\n';
    } else {
        stream << "Printing top " << m << " genotype priors " << '\n';
    }
    using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
    std::vector<std::pair<GenotypeReference, double>> v {};
    v.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(priors),
                   std::back_inserter(v), [] (const auto& g, const auto& p) {
                       return std::make_pair(std::cref(g), p);
                   });
    const auto mth = std::next(std::begin(v), m);
    std::partial_sort(std::begin(v), mth, std::end(v),
                      [] (const auto& lhs, const auto& rhs) {
                          return lhs.second > rhs.second;
                      });
    std::for_each(std::begin(v), mth,
                  [&] (const auto& p) {
                      print_variant_alleles(stream, p.first);
                      stream << " " << p.second << '\n';
                  });
}

void print_genotype_priors(const std::vector<Genotype<Haplotype>>& genotypes,
                           const std::vector<double>& priors, const std::size_t n)
{
    print_genotype_priors(std::cout, genotypes, priors, n);
}

void log_genotype_likelihoods(boost::optional<logging::DebugLogger>& debug_log,
                              boost::optional<logging::TraceLogger>& trace_log,
                              const std::vector<Genotype<Haplotype>>& genotypes,
                              const std::vector<double>& likelihoods)
{
    if (debug_log) debug::print_genotype_likelihoods(stream(*debug_log), genotypes, likelihoods);
    if (trace_log) debug::print_genotype_likelihoods(stream(*trace_log), genotypes, likelihoods, genotypes.size());
}

template <typename S>
void print_genotype_likelihoods(S&& stream, const std::vector<Genotype<Haplotype>>& genotypes,
                                const std::vector<double>& likelihoods, std::size_t n)
{
    assert(genotypes.size() == likelihoods.size());
    const auto m = std::min(n, genotypes.size());
    if (m == genotypes.size()) {
        stream << "Printing all genotype likelihoods " << '\n';
    } else {
        stream << "Printing top " << m << " genotype likelihoods " << '\n';
    }
    using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
    std::vector<std::pair<GenotypeReference, double>> v {};
    v.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(likelihoods),
                   std::back_inserter(v), [] (const auto& g, const auto& p) {
                       return std::make_pair(std::cref(g), p);
                   });
    const auto mth = std::next(std::begin(v), m);
    std::partial_sort(std::begin(v), mth, std::end(v),
                      [] (const auto& lhs, const auto& rhs) {
                          return lhs.second > rhs.second;
                      });
    std::for_each(std::begin(v), mth,
                  [&] (const auto& p) {
                      print_variant_alleles(stream, p.first);
                      stream << " " << p.second << '\n';
                  });
}

void print_genotype_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                const std::vector<double>& likelihoods, std::size_t n)
{
    print_genotype_likelihoods(std::cout, genotypes, likelihoods, n);
}

} // namespace debug
} // namesapce model
} // namespace octopus
