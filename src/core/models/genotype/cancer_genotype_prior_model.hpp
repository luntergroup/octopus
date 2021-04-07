// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cancer_genotype_prior_model_hpp
#define cancer_genotype_prior_model_hpp

#include <functional>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <cassert>

#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/cancer_genotype.hpp"
#include "genotype_prior_model.hpp"
#include "../mutation/somatic_mutation_model.hpp"
#include "utils/maths.hpp"

namespace octopus {

class CancerGenotypePriorModel
{
public:
    using LogProbability = double;
    
    CancerGenotypePriorModel() = delete;
    
    CancerGenotypePriorModel(const GenotypePriorModel& germline_model,
                             SomaticMutationModel mutation_model);
    
    CancerGenotypePriorModel(const CancerGenotypePriorModel&)            = default;
    CancerGenotypePriorModel& operator=(const CancerGenotypePriorModel&) = default;
    CancerGenotypePriorModel(CancerGenotypePriorModel&&)                 = default;
    CancerGenotypePriorModel& operator=(CancerGenotypePriorModel&&)      = default;
    
    ~CancerGenotypePriorModel() = default;
    
    const GenotypePriorModel& germline_model() const noexcept;
    SomaticMutationModel& mutation_model() noexcept;
    const SomaticMutationModel& mutation_model() const noexcept;
    
    LogProbability evaluate(const CancerGenotype<IndexedHaplotype<>>& genotype) const;

private:
    std::reference_wrapper<const GenotypePriorModel> germline_model_;
    SomaticMutationModel mutation_model_;
    
    // p(somatic | germline)
    template <typename G, typename H>
    LogProbability ln_probability_of_somatic_given_genotype(const H& somatic, const G& germline) const;
    LogProbability ln_probability_of_somatic_given_haplotype(const IndexedHaplotype<>& somatic, const IndexedHaplotype<>& germline) const;
};

// p(somatic | germline) = 1 / M sum [k = 1 -> M] p(somatic | germline_k) (M = germline ploidy)
template <typename G, typename H>
CancerGenotypePriorModel::LogProbability
CancerGenotypePriorModel::ln_probability_of_somatic_given_genotype(const H& somatic, const G& germline) const
{
    assert(ploidy(germline) > 0);
    switch (ploidy(germline)) {
        case 1: return ln_probability_of_somatic_given_haplotype(somatic, germline[0]);
        case 2:
        {
            const static LogProbability ln2 {std::log(2)};
            const auto a = ln_probability_of_somatic_given_haplotype(somatic, germline[0]);
            const auto b = ln_probability_of_somatic_given_haplotype(somatic, germline[1]);
            return maths::log_sum_exp(a, b) - ln2;
        }
        case 3:
        {
            const static LogProbability ln3 {std::log(3)};
            const auto a = ln_probability_of_somatic_given_haplotype(somatic, germline[0]);
            const auto b = ln_probability_of_somatic_given_haplotype(somatic, germline[1]);
            const auto c = ln_probability_of_somatic_given_haplotype(somatic, germline[3]);
            return maths::log_sum_exp(a, b, c) - ln3;
        }
        default:
        {
            std::vector<LogProbability> tmp(ploidy(germline));
            std::transform(std::cbegin(germline), std::cend(germline), std::begin(tmp),
                           [&] (const auto& haplotype) {
                               return ln_probability_of_somatic_given_haplotype(somatic, haplotype);
                           });
            return maths::log_sum_exp(tmp) - std::log(ploidy(germline));
        }
    }
}

// non-member methods

template <typename Container1, typename Container2>
Container2&
evaluate(const Container1& genotypes, const CancerGenotypePriorModel& model, Container2& result,
         const bool normalise = false, const bool add = false)
{
    if (add) {
        assert(result.size() == genotypes.size());
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(result), std::begin(result),
                       [&] (const auto& genotype, auto curr) { return curr + model.evaluate(genotype); });
    } else {
        result.resize(genotypes.size());
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [&] (const auto& genotype) { return model.evaluate(genotype); });
    }
    if (normalise) maths::normalise_logs(result);
    return result;
}

template <typename Container>
auto evaluate(const Container& genotypes, const CancerGenotypePriorModel& model, const bool normalise = false)
{
    std::vector<CancerGenotypePriorModel::LogProbability> result(genotypes.size());
    evaluate(genotypes, model, result, normalise, false);
    return result;
}

} // namespace octopus

#endif
