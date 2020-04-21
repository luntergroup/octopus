// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "constant_mixture_genotype_likelihood_model.hpp"

#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <array>
#include <limits>
#include <cassert>

#include "utils/maths.hpp"

namespace octopus { namespace model {

ConstantMixtureGenotypeLikelihoodModel::ConstantMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods)
: likelihoods_ {likelihoods}
{}

const HaplotypeLikelihoodArray& ConstantMixtureGenotypeLikelihoodModel::cache() const noexcept
{
    return likelihoods_;
}

// ln p(read | genotype)  = ln sum {haplotype in genotype} p(read | haplotype) - ln ploidy
// ln p(reads | genotype) = sum {read in reads} ln p(read | genotype)
ConstantMixtureGenotypeLikelihoodModel::LogProbability
ConstantMixtureGenotypeLikelihoodModel::evaluate(const Genotype<Haplotype>& genotype) const
{
    assert(likelihoods_.is_primed());
    // These cases are just for optimisation
    switch (genotype.ploidy()) {
        case 0:
            return 0.0;
        case 1:
            return evaluate_haploid(genotype);
        case 2:
            return evaluate_diploid(genotype);
        case 3:
            return evaluate_triploid(genotype);
        default:
            return evaluate_polyploid(genotype);
    }
}

namespace {

template <typename T = double>
static constexpr auto ln(const unsigned n)
{
    constexpr std::array<T, 11> lnLookup {
    std::numeric_limits<T>::infinity(),
    0.0,
    0.693147180559945309417232121458176568075500134360255254120,
    1.098612288668109691395245236922525704647490557822749451734,
    1.386294361119890618834464242916353136151000268720510508241,
    1.609437912434100374600759333226187639525601354268517721912,
    1.791759469228055000812477358380702272722990692183004705855,
    1.945910149055313305105352743443179729637084729581861188459,
    2.079441541679835928251696364374529704226500403080765762362,
    2.197224577336219382790490473845051409294981115645498903469,
    2.302585092994045684017991454684364207601101488628772976033
    };
    return lnLookup[n];
}

} // namespace

ConstantMixtureGenotypeLikelihoodModel::LogProbability
ConstantMixtureGenotypeLikelihoodModel::evaluate(const Genotype<IndexedHaplotype<>>& genotype) const
{
    switch (genotype.ploidy()) {
        case 0:
            return 0.0;
        case 1:
            return evaluate_haploid(genotype);
        case 2:
            return evaluate_diploid(genotype);
        default:
            return evaluate_polyploid(genotype);
    }
}

// private methods

ConstantMixtureGenotypeLikelihoodModel::LogProbability
ConstantMixtureGenotypeLikelihoodModel::evaluate_haploid(const Genotype<Haplotype>& genotype) const
{
    const auto& log_likelihoods = likelihoods_[genotype[0]];
    return std::accumulate(std::cbegin(log_likelihoods), std::cend(log_likelihoods), LogProbability {0});
}

ConstantMixtureGenotypeLikelihoodModel::LogProbability
ConstantMixtureGenotypeLikelihoodModel::evaluate_diploid(const Genotype<Haplotype>& genotype) const
{
    const auto& log_likelihoods1 = likelihoods_[genotype[0]];
    if (is_homozygous(genotype)) {
        return std::accumulate(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1), LogProbability {0});
    }
    const auto& log_likelihoods2 = likelihoods_[genotype[1]];
    return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                              std::cbegin(log_likelihoods2), LogProbability {0}, std::plus<> {},
                              [] (const auto a, const auto b) -> LogProbability {
                                  return maths::log_sum_exp(a, b) - ln<decltype(a)>(2);
                              });
}

ConstantMixtureGenotypeLikelihoodModel::LogProbability
ConstantMixtureGenotypeLikelihoodModel::evaluate_triploid(const Genotype<Haplotype>& genotype) const
{
    using std::cbegin; using std::cend;
    const auto& log_likelihoods1 = likelihoods_[genotype[0]];
    if (is_homozygous(genotype)) {
        return std::accumulate(cbegin(log_likelihoods1), cend(log_likelihoods1), LogProbability {0});
    }
    if (zygosity(genotype) == 3) {
        const auto& log_likelihoods2 = likelihoods_[genotype[1]];
        const auto& log_likelihoods3 = likelihoods_[genotype[2]];
        return maths::inner_product(cbegin(log_likelihoods1), cend(log_likelihoods1),
                                    cbegin(log_likelihoods2), cbegin(log_likelihoods3),
                                    LogProbability {0}, std::plus<> {},
                                    [] (const auto a, const auto b, const auto c) -> LogProbability {
                                        return maths::log_sum_exp(a, b, c) - ln<decltype(a)>(3);
                                    });
    }
    if (genotype[0] != genotype[1]) {
        const auto& log_likelihoods2 = likelihoods_[genotype[1]];
        return std::inner_product(cbegin(log_likelihoods1), cend(log_likelihoods1),
                                  cbegin(log_likelihoods2), LogProbability {0}, std::plus<> {},
                                  [] (const auto a, const auto b) -> LogProbability {
                                      return maths::log_sum_exp(a, ln<decltype(a)>(2) + b) - ln<decltype(a)>(3);
                                  });
    }
    const auto& log_likelihoods3 = likelihoods_[genotype[2]];
    return std::inner_product(cbegin(log_likelihoods1), cend(log_likelihoods1),
                              cbegin(log_likelihoods3), LogProbability {0}, std::plus<> {},
                              [] (const auto a, const auto b) -> LogProbability {
                                  return maths::log_sum_exp(ln<decltype(a)>(2) + a, b) - ln<decltype(a)>(3);
                              });
}

ConstantMixtureGenotypeLikelihoodModel::LogProbability
ConstantMixtureGenotypeLikelihoodModel::evaluate_polyploid(const Genotype<Haplotype>& genotype) const
{
    const auto ploidy = genotype.ploidy();
    const auto ln_ploidy = std::log(ploidy);
    const auto& log_likelihoods1 = likelihoods_[genotype[0]];
    if (is_homozygous(genotype)) {
        return std::accumulate(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1), LogProbability {0});
    }
    likelihood_refs_.reserve(ploidy);
    likelihood_refs_.push_back(log_likelihoods1);
    std::transform(std::next(std::cbegin(genotype)), std::cend(genotype), std::back_inserter(likelihood_refs_),
                   [this] (const auto& haplotype) -> const HaplotypeLikelihoodArray::LikelihoodVector& {
                       return likelihoods_[haplotype]; });
    LogProbability result {0};
    const auto num_likelihoods = likelihood_refs_.front().get().size();
    buffer_.resize(ploidy);
    for (std::size_t read_idx {0}; read_idx < num_likelihoods; ++read_idx) {
        std::transform(std::cbegin(likelihood_refs_), std::cend(likelihood_refs_), std::begin(buffer_),
                       [read_idx] (const auto& likelihoods) noexcept { return likelihoods.get()[read_idx]; });
        result += maths::log_sum_exp(buffer_) - ln_ploidy;
    }
    likelihood_refs_.clear();
    return result;
}

ConstantMixtureGenotypeLikelihoodModel::LogProbability
ConstantMixtureGenotypeLikelihoodModel::evaluate_haploid(const Genotype<IndexedHaplotype<>>& genotype) const
{
    const auto& log_likelihoods = likelihoods_[genotype[0]];
    return std::accumulate(std::cbegin(log_likelihoods), std::cend(log_likelihoods), LogProbability {0});
}

ConstantMixtureGenotypeLikelihoodModel::LogProbability
ConstantMixtureGenotypeLikelihoodModel::evaluate_diploid(const Genotype<IndexedHaplotype<>>& genotype) const
{
    const auto& log_likelihoods1 = likelihoods_[genotype[0]];
    if (is_homozygous(genotype)) {
        return std::accumulate(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1), LogProbability {0});
    } else {
        const auto& log_likelihoods2 = likelihoods_[genotype[1]];
        return std::inner_product(std::cbegin(log_likelihoods1), std::cend(log_likelihoods1),
                                  std::cbegin(log_likelihoods2), LogProbability {0}, std::plus<> {},
                                  [] (const auto a, const auto b) -> LogProbability {
                                      return maths::log_sum_exp(a, b) - ln<decltype(a)>(2);
                                  });
    }
}

ConstantMixtureGenotypeLikelihoodModel::LogProbability 
ConstantMixtureGenotypeLikelihoodModel::evaluate_polyploid(const Genotype<IndexedHaplotype<>>& genotype) const
{
    assert(likelihoods_.is_primed());
    const auto ln_ploidy = std::log(genotype.ploidy());
    buffer_.resize(genotype.ploidy());
    LogProbability result {0};
    const auto num_likelihoods = likelihoods_.num_likelihoods();
    for (std::size_t read_idx {0}; read_idx < num_likelihoods; ++read_idx) {
        std::transform(std::cbegin(genotype), std::cend(genotype), std::begin(buffer_),
                       [&] (const auto& haplotype) noexcept { return likelihoods_[haplotype][read_idx]; });
        result += maths::log_sum_exp(buffer_) - ln_ploidy;
    }
    return result;
}

} // namespace model
} // namespace octopus
