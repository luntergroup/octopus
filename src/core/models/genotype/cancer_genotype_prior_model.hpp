// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cancer_genotype_prior_model_hpp
#define cancer_genotype_prior_model_hpp

#include <type_traits>
#include <functional>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <cassert>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/cancer_genotype.hpp"
#include "genotype_prior_model.hpp"
#include "../mutation/somatic_mutation_model.hpp"
#include "utils/maths.hpp"

namespace octopus {

class CancerGenotypePriorModel
{
public:
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
    
    double evaluate(const CancerGenotype<Haplotype>& genotype) const;
    double evaluate(const std::vector<unsigned>& germline_indices, unsigned somatic_index) const;

private:
    std::reference_wrapper<const GenotypePriorModel> germline_model_;
    SomaticMutationModel mutation_model_;
    
    // p(somatic | germline)
    template <typename G, typename H>
    double ln_probability_of_somatic_given_genotype(const H& somatic, const G& germline) const;
    double ln_probability_of_somatic_given_haplotype(const Haplotype& somatic, const Haplotype& germline) const;
    double ln_probability_of_somatic_given_haplotype(unsigned somatic_index, unsigned germline_index) const;
};

namespace detail {

inline auto get_ploidy(const Genotype<Haplotype>& genotype) noexcept
{
    return genotype.ploidy();
}

inline auto get_ploidy(const std::vector<unsigned>& genotype_indices) noexcept
{
    return static_cast<unsigned>(genotype_indices.size());
}

} // namespace detail

// p(somatic | germline) = 1 / M sum [k = 1 -> M] p(somatic | germline_k) (M = germline ploidy)
template <typename G, typename H>
double CancerGenotypePriorModel::ln_probability_of_somatic_given_genotype(const H& somatic, const G& germline) const
{
    const auto ploidy = detail::get_ploidy(germline);
    assert(ploidy > 0);
    switch (ploidy) {
        case 1: return ln_probability_of_somatic_given_haplotype(somatic, germline[0]);
        case 2:
        {
            const static double ln2 {std::log(2)};
            const auto a = ln_probability_of_somatic_given_haplotype(somatic, germline[0]);
            const auto b = ln_probability_of_somatic_given_haplotype(somatic, germline[1]);
            return maths::log_sum_exp(a, b) - ln2;
        }
        case 3:
        {
            const static double ln3 {std::log(3)};
            const auto a = ln_probability_of_somatic_given_haplotype(somatic, germline[0]);
            const auto b = ln_probability_of_somatic_given_haplotype(somatic, germline[1]);
            const auto c = ln_probability_of_somatic_given_haplotype(somatic, germline[3]);
            return maths::log_sum_exp(a, b, c) - ln3;
        }
        default:
        {
            std::vector<double> tmp(ploidy);
            std::transform(std::cbegin(germline), std::cend(germline), std::begin(tmp),
                           [this, &somatic] (const auto& haplotype) {
                               return ln_probability_of_somatic_given_haplotype(somatic, haplotype);
                           });
            return maths::log_sum_exp(tmp) - std::log(ploidy);
        }
    }
}

// non-member methods

template <typename Container>
auto calculate_log_priors(const Container& genotypes, const CancerGenotypePriorModel& model)
{
    static_assert(std::is_same<typename Container::value_type, CancerGenotype<Haplotype>>::value,
                  "genotypes must contain CancerGenotype<Haplotype>'s");
    std::vector<double> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&model] (const auto& genotype) { return model.evaluate(genotype); });
    maths::normalise_logs(result);
    return result;
}

} // namespace octopus

#endif
