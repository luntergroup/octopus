// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "trio_model.hpp"

#include <functional>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "utils/maths.hpp"
#include "germline_likelihood_model.hpp"

namespace octopus { namespace model {

TrioModel::TrioModel(const Trio& trio,
              const CoalescentModel& genotype_prior_model,
              const DeNovoModel& mutation_model,
              boost::optional<logging::DebugLogger> debug_log)
: trio_ {trio}
, genotype_prior_model_ {genotype_prior_model}
, mutation_model_ {mutation_model}
, debug_log_ {debug_log}
{}

namespace {

template <typename Container>
void clear(Container& c)
{
    c.clear();
    c.shrink_to_fit();
}

using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;

bool operator==(const GenotypeReference lhs, const GenotypeReference rhs)
{
    return lhs.get() == rhs.get();
}

bool operator<(const GenotypeReference lhs, const GenotypeReference rhs)
{
    return GenotypeLess{}(lhs.get(), rhs.get());
}

struct GenotypeRefProbabilityPair
{
    GenotypeReference genotype;
    double probability;
};

bool operator<(const GenotypeRefProbabilityPair& lhs, const GenotypeRefProbabilityPair& rhs) noexcept
{
    return lhs.probability < rhs.probability;
}
bool operator>(const GenotypeRefProbabilityPair& lhs, const GenotypeRefProbabilityPair& rhs) noexcept
{
    return lhs.probability > rhs.probability;
}

struct ParentsProbabilityPair
{
    GenotypeReference maternal, paternal;
    double probability;
};

bool operator<(const ParentsProbabilityPair& lhs, const ParentsProbabilityPair& rhs) noexcept
{
    return lhs.probability < rhs.probability;
}
bool operator>(const ParentsProbabilityPair& lhs, const ParentsProbabilityPair& rhs) noexcept
{
    return lhs.probability > rhs.probability;
}

auto compute_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                         const GermlineLikelihoodModel& model)
{
    std::vector<GenotypeRefProbabilityPair> result {};
    result.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::back_inserter(result),
                   [&model] (const auto& genotype) {
                       return GenotypeRefProbabilityPair {genotype, model.evaluate(genotype)};
                   });
    return result;
}

template <typename T>
void reduce(std::vector<T>& zipped, const std::size_t max_keep)
{
    const auto first_erase = std::next(std::begin(zipped), std::min(zipped.size(), max_keep));
    std::partial_sort(std::begin(zipped), first_erase, std::end(zipped), std::greater<> {});
    zipped.erase(first_erase, std::end(zipped));
}

double probability_of_parents(const Genotype<Haplotype>& mother,
                              const Genotype<Haplotype>& father,
                              const CoalescentModel& model)
{
    thread_local std::vector<std::reference_wrapper<const Haplotype>> parental_haplotypes;
    parental_haplotypes.reserve(mother.ploidy() + father.ploidy());
    parental_haplotypes.assign(std::cbegin(mother), std::cend(mother));
    parental_haplotypes.insert(std::cbegin(parental_haplotypes), std::cbegin(father), std::cend(father));
    return model.evaluate(parental_haplotypes);
}

auto join(const std::vector<GenotypeRefProbabilityPair>& maternal,
          const std::vector<GenotypeRefProbabilityPair>& paternal,
          const CoalescentModel& model)
{
    std::vector<ParentsProbabilityPair> result {};
    result.reserve(maternal.size() * paternal.size());
    for (const auto& m : maternal) {
        for (const auto& p : paternal) {
            result.push_back({m.genotype, p.genotype,
                                m.probability + p.probability
                                + probability_of_parents(m.genotype, p.genotype, model)});
        }
    }
    return result;
}

bool all_diploid(const Genotype<Haplotype>& child,
                 const Genotype<Haplotype>& mother,
                 const Genotype<Haplotype>& father)
{
    return is_diploid(child) && is_diploid(mother) && is_diploid(father);
}

double probability_of_child_given_parent(const Haplotype& child,
                                         const Genotype<Haplotype>& parent,
                                         const DeNovoModel& model)
{
    static const double ln2 {std::log(2)};
    const auto p1 = model.evaluate(child, parent[0]);
    const auto p2 = model.evaluate(child, parent[1]);
    return maths::log_sum_exp(p1, p2) - ln2;
}

double probability_of_child_given_parents(const Haplotype& child_from_mother,
                                          const Haplotype& child_from_father,
                                          const Genotype<Haplotype>& mother,
                                          const Genotype<Haplotype>& father,
                                          const DeNovoModel& model)
{
    return probability_of_child_given_parent(child_from_mother, mother, model)
            + probability_of_child_given_parent(child_from_father, father, model);
}

double probability_of_child_given_parents(const Genotype<Haplotype>& child,
                                          const Genotype<Haplotype>& mother,
                                          const Genotype<Haplotype>& father,
                                          const DeNovoModel& model)
{
    if (all_diploid(child, mother, father)) {
        static const double ln2 {std::log(2)};
        const auto p1 = probability_of_child_given_parents(child[0], child[1], mother, father, model);
        const auto p2 = probability_of_child_given_parents(child[1], child[0], mother, father, model);
//        std::cout << "Child = ";
//        debug::print_variant_alleles(child);
//        std::cout << "\nMother = ";
//        debug::print_variant_alleles(mother);
//        std::cout << "\nFather = ";
//        debug::print_variant_alleles(father);
//        std::cout << "\n" << (maths::log_sum_exp(p1, p2) - ln2) << '\n';
        return maths::log_sum_exp(p1, p2) - ln2;
    }
    return 0; // TODO
}

using JointProbability = TrioModel::Latents::JointProbability;

auto join(const std::vector<ParentsProbabilityPair>& parents,
          const std::vector<GenotypeRefProbabilityPair>& child,
          const DeNovoModel& model)
{
    std::vector<JointProbability> result {};
    result.reserve(parents.size() * child.size());
    for (const auto& p : parents) {
        for (const auto& c : child) {
            result.push_back({p.maternal, p.paternal, c.genotype,
                                p.probability + c.probability
                                + probability_of_child_given_parents(p.maternal, p.paternal,
                                                                     c.genotype, model)});
        }
    }
    return result;
}

auto extract_probabilities(const std::vector<JointProbability>& joint_likelihoods)
{
    std::vector<double> result(joint_likelihoods.size());
    std::transform(std::cbegin(joint_likelihoods), std::cend(joint_likelihoods),
                   std::begin(result), [] (const auto& p) { return p.probability; });
    return result;
}

auto normalise_exp(std::vector<JointProbability>& joint_likelihoods)
{
    auto likelihoods = extract_probabilities(joint_likelihoods);
    const auto norm = maths::normalise_exp(likelihoods);
    auto iter = std::cbegin(likelihoods);
    for (auto& p : joint_likelihoods) p.probability = *iter++;
    return norm;
}
    
} // namespace

void print(const std::vector<GenotypeRefProbabilityPair>& ps)
{
    for (const auto& p : ps) {
        debug::print_variant_alleles(p.genotype);
        std::cout << " " << p.probability << "\n";
    }
}
void print(const std::vector<ParentsProbabilityPair>& ps)
{
    for (const auto& p : ps) {
        debug::print_variant_alleles(p.maternal);
        std::cout << " | ";
        debug::print_variant_alleles(p.paternal);
        std::cout << " " << p.probability << "\n";
    }
}
void print(const std::vector<JointProbability>& ps)
{
    for (const auto& p : ps) {
        debug::print_variant_alleles(p.maternal);
        std::cout << " | ";
        debug::print_variant_alleles(p.paternal);
        std::cout << " | ";
        debug::print_variant_alleles(p.child);
        std::cout << " " << p.probability << "\n";
    }
}
void print_top(std::vector<JointProbability> ps)
{
    std::sort(std::begin(ps), std::end(ps),
              [] (const auto& lhs, const auto& rhs) {
                  return lhs.probability > rhs.probability;
              });
    debug::print_variant_alleles(ps.front().maternal);
    std::cout << " | ";
    debug::print_variant_alleles(ps.front().paternal);
    std::cout << " | ";
    debug::print_variant_alleles(ps.front().child);
    std::cout << " " << ps.front().probability << "\n";
}

TrioModel::InferredLatents
TrioModel::evaluate(const GenotypeVector& maternal_genotypes,
                    const GenotypeVector& paternal_genotypes,
                    const GenotypeVector& child_genotypes,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!maternal_genotypes.empty() && !paternal_genotypes.empty() && !child_genotypes.empty());
    const GermlineLikelihoodModel likelihood_model {haplotype_likelihoods};
    
    haplotype_likelihoods.prime(trio_.mother());
    auto maternal_likelihoods = compute_likelihoods(maternal_genotypes, likelihood_model);
    reduce(maternal_likelihoods, max_search_size_);
    
    haplotype_likelihoods.prime(trio_.father());
    auto paternal_likelihoods = compute_likelihoods(paternal_genotypes, likelihood_model);
    reduce(paternal_likelihoods, max_search_size_);
    
//    std::cout << "maternal" << '\n';
//    print(maternal_likelihoods);
//    std::cout << "paternal" << '\n';
//    print(paternal_likelihoods);
    
    auto parents_joint_likelihoods = join(maternal_likelihoods, paternal_likelihoods, genotype_prior_model_);
    reduce(parents_joint_likelihoods, max_search_size_);
    clear(maternal_likelihoods);
    clear(paternal_likelihoods);
    
//    std::cout << "joint parents" << '\n';
//    print(parents_joint_likelihoods);
    
    haplotype_likelihoods.prime(trio_.child());
    auto child_likelihoods = compute_likelihoods(child_genotypes, likelihood_model);
    reduce(child_likelihoods, max_search_size_);
    
//    std::cout << "child" << '\n';
//    print(child_likelihoods);
    
    auto joint_likelihoods = join(parents_joint_likelihoods, child_likelihoods, mutation_model_);
    clear(parents_joint_likelihoods);
    clear(child_likelihoods);
    const auto evidence = normalise_exp(joint_likelihoods);
    
//    std::cout << "trio joint" << '\n';
//    print(joint_likelihoods);
//    print_top(joint_likelihoods);
//    std::cout << '\n';
    
    return {std::move(joint_likelihoods), evidence};
}

} // namespace model
} // namespace octopus
