// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "inheritance_model.hpp"

#include <utility>

#include "utils/maths.hpp"

namespace octopus {

InheritanceModel::InheritanceModel(DeNovoModel mutation_model)
: mutation_model_ {std::move(mutation_model)}
{}

void InheritanceModel::prime(const MappableBlock<Haplotype>& haplotypes)
{
    mutation_model_.prime(haplotypes);
}

void InheritanceModel::unprime() noexcept
{
    mutation_model_.unprime();
}

bool InheritanceModel::is_primed() const noexcept
{
    return mutation_model_.is_primed();
}

namespace {

template <typename G>
bool all_haploid(const G& child, const G& mother, const G& father) noexcept
{
    return is_haploid(child) && is_haploid(mother) && is_haploid(father);
}

template <typename G>
bool all_diploid(const G& child, const G& mother, const G& father) noexcept
{
    return is_diploid(child) && is_diploid(mother) && is_diploid(father);
}

template <typename G>
bool all_triploid(const G& child, const G& mother, const G& father) noexcept
{
    return is_triploid(child) && is_triploid(mother) && is_triploid(father);
}

} // namespace

// p(offspring | parent)
InheritanceModel::LogProbability 
InheritanceModel::evaluate(const Genotype<Haplotype>& offspring, const Genotype<Haplotype>& parent) const
{
    return 0;
}

InheritanceModel::LogProbability 
InheritanceModel::evaluate(const Genotype<IndexedHaplotype<>>& offspring, const Genotype<IndexedHaplotype<>>& parent) const
{
    if (is_diploid(parent)) {
        if (is_diploid(offspring)) {
            static const double ln4 {std::log(4)};
            const auto p1 = mutation_model_.evaluate(offspring[0], parent[0]);
            const auto p2 = mutation_model_.evaluate(offspring[0], parent[1]);
            const auto p3 = mutation_model_.evaluate(offspring[1], parent[0]);
            const auto p4 = mutation_model_.evaluate(offspring[1], parent[1]);
            return maths::log_sum_exp({p1, p2, p3, p4}) - ln4;
        } else if (is_haploid(offspring)) {
            static const double ln2 {std::log(2)};
            const auto p1 = mutation_model_.evaluate(offspring[0], parent[0]);
            const auto p2 = mutation_model_.evaluate(offspring[0], parent[1]);
            return maths::log_sum_exp(p1, p2) - ln2;
        } else {
            return 0;
        }
    } else if (is_haploid(parent)) {
        if (is_diploid(offspring)) {
            static const double ln2 {std::log(2)};
            const auto p1 = mutation_model_.evaluate(offspring[0], parent[0]);
            const auto p2 = mutation_model_.evaluate(offspring[1], parent[0]);
            return maths::log_sum_exp(p1, p2) - ln2;
        } else if (is_haploid(offspring)) {
            return 0; // TODO: what to do here?
        } else {
            return 0;
        }
    } else {
        return 0;
    }
}

namespace {

template <typename H, typename G>
double probability_of_child_given_haploid_parent(const H& child, const G& parent,
                                                 const DeNovoModel& mutation_model)
{
    return mutation_model.evaluate(child, parent[0]);
}

template <typename H, typename G>
double probability_of_child_given_diploid_parent(const H& child, const G& parent,
                                                 const DeNovoModel& mutation_model)
{
    static const double ln2 {std::log(2)};
    const auto p1 = mutation_model.evaluate(child, parent[0]);
    const auto p2 = mutation_model.evaluate(child, parent[1]);
    return maths::log_sum_exp(p1, p2) - ln2;
}

template <typename H, typename G>
double probability_of_child_given_diploid_parents(const H& child_from_mother, const H& child_from_father,
                                                  const G& mother, const G& father,
                                                  const DeNovoModel& mutation_model)
{
    return probability_of_child_given_diploid_parent(child_from_mother, mother, mutation_model)
            + probability_of_child_given_diploid_parent(child_from_father, father, mutation_model);
}

template <typename H, typename G>
double probability_of_child_given_triploid_parent(const H& child, const G& parent,
                                                  const DeNovoModel& mutation_model)
{
    static const double ln3 {std::log(3)};
    const auto p1 = mutation_model.evaluate(child, parent[0]);
    const auto p2 = mutation_model.evaluate(child, parent[1]);
    const auto p3 = mutation_model.evaluate(child, parent[2]);
    return maths::log_sum_exp(p1, p2, p3) - ln3;
}

template <typename H, typename G>
double probability_of_child_given_triploid_parents(const H& child1_from_mother,
                                                   const H& child2_from_mother,
                                                   const H& child_from_father,
                                                   const G& mother, const G& father,
                                                   const DeNovoModel& mutation_model)
{
    return probability_of_child_given_triploid_parent(child1_from_mother, mother, mutation_model)
           + probability_of_child_given_triploid_parent(child2_from_mother, mother, mutation_model)
           + probability_of_child_given_triploid_parent(child_from_father, father, mutation_model);
}

template <unsigned ChildPloidy, unsigned MotherPloidy, unsigned FatherPloidy>
struct ProbabilityOfChildGivenParents
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        return 0;
    }
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<2, 2, 2>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        static const double ln2 {std::log(2)};
        const auto p1 = probability_of_child_given_diploid_parents(child[0], child[1], mother, father, mutation_model);
        const auto p2 = probability_of_child_given_diploid_parents(child[1], child[0], mother, father, mutation_model);
        return maths::log_sum_exp(p1, p2) - ln2;
    }
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<3, 3, 3>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        static const double ln6 {std::log(6)};
        const auto p1 = probability_of_child_given_triploid_parents(child[0], child[1], child[2], mother, father, mutation_model);
        const auto p2 = probability_of_child_given_triploid_parents(child[0], child[2], child[1], mother, father, mutation_model);
        const auto p3 = probability_of_child_given_triploid_parents(child[1], child[0], child[2], mother, father, mutation_model);
        const auto p4 = probability_of_child_given_triploid_parents(child[1], child[2], child[0], mother, father, mutation_model);
        const auto p5 = probability_of_child_given_triploid_parents(child[2], child[0], child[1], mother, father, mutation_model);
        const auto p6 = probability_of_child_given_triploid_parents(child[2], child[1], child[0], mother, father, mutation_model);
        return maths::log_sum_exp({p1, p2, p3, p4, p5, p6}) - ln6;
    }
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<2, 2, 1>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        static const double ln2 {std::log(2)};
        const auto p1 = probability_of_child_given_diploid_parent(child[0], mother, mutation_model);
        const auto p2 = probability_of_child_given_diploid_parent(child[1], mother, mutation_model);
        const auto p3 = probability_of_child_given_haploid_parent(child[0], father, mutation_model);
        const auto p4 = probability_of_child_given_haploid_parent(child[1], father, mutation_model);
        return maths::log_sum_exp(p1 + p4, p2 + p3) - ln2;
    }
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<1, 2, 1>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        return probability_of_child_given_diploid_parent(child[0], mother, mutation_model);
    }
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<1, 0, 1>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        return probability_of_child_given_haploid_parent(child[0], father, mutation_model);
    }
    const DeNovoModel& mutation_model;
};

template <> struct ProbabilityOfChildGivenParents<1, 1, 1>
{
    ProbabilityOfChildGivenParents(const DeNovoModel& mutation_model) : mutation_model {mutation_model} {}
    template <typename G>
    double operator()(const G& child, const G& mother, const G& father)
    {
        static const double ln2 {std::log(2)};
        const auto p1 = probability_of_child_given_haploid_parent(child[0], mother, mutation_model);
        const auto p2 = probability_of_child_given_haploid_parent(child[0], father, mutation_model);
        return maths::log_sum_exp(p1, p2) - ln2;
    }
    const DeNovoModel& mutation_model;
};

} // namespace

// p(offspring | mother, father)
InheritanceModel::LogProbability 
InheritanceModel::evaluate(const Genotype<Haplotype>& offspring, const Genotype<Haplotype>& mother, const Genotype<Haplotype>& father) const
{
    return 0;
}

InheritanceModel::LogProbability 
InheritanceModel::evaluate(const Genotype<IndexedHaplotype<>>& offspring, const Genotype<IndexedHaplotype<>>& mother, const Genotype<IndexedHaplotype<>>& father) const
{
    const auto mother_ploidy    = mother.ploidy();
    const auto father_ploidy    = father.ploidy();
    const auto offspring_ploidy = offspring.ploidy();
    if (offspring_ploidy == 1) {
        if (father_ploidy == 1) {
            if (mother_ploidy == 0) {
                return ProbabilityOfChildGivenParents<1, 0, 1>{mutation_model_}(offspring, mother, father);
            }
            if (mother_ploidy == 1) {
                return ProbabilityOfChildGivenParents<1, 1, 1>{mutation_model_}(offspring, mother, father);
            }
            if (mother_ploidy == 2) {
                return ProbabilityOfChildGivenParents<1, 2, 1>{mutation_model_}(offspring, mother, father);
            }
        }
    } else if (offspring_ploidy == 2) {
        if (mother_ploidy == 2) {
            if (father_ploidy == 1) {
                return ProbabilityOfChildGivenParents<2, 2, 1>{mutation_model_}(offspring, mother, father);
            }
            if (father_ploidy == 2) {
                return ProbabilityOfChildGivenParents<2, 2, 2>{mutation_model_}(offspring, mother, father);
            }
        }
    }
    throw std::runtime_error {"TrioModel: unimplemented joint probability function"};
}

// p(twin1, twin2 | parent)
InheritanceModel::LogProbability 
InheritanceModel::evaluate_twins(const Genotype<Haplotype>& twin1, const Genotype<Haplotype>& twin2, const Genotype<Haplotype>& parent) const
{
    return 0;
}

InheritanceModel::LogProbability 
InheritanceModel::evaluate_twins(const Genotype<IndexedHaplotype<>>& twin1, const Genotype<IndexedHaplotype<>>& twin2, const Genotype<IndexedHaplotype<>>& parent) const
{
    return 0;
}

// p(twin1, twin2 | mother, father)
InheritanceModel::LogProbability 
InheritanceModel::evaluate_twins(const Genotype<Haplotype>& twin1, const Genotype<Haplotype>& twin2, const Genotype<Haplotype>& mother, const Genotype<Haplotype>& father) const
{
    return 0;
}

InheritanceModel::LogProbability 
InheritanceModel::evaluate_twins(const Genotype<IndexedHaplotype<>>& twin1, const Genotype<IndexedHaplotype<>>& twin2, const Genotype<IndexedHaplotype<>>& mother, const Genotype<IndexedHaplotype<>>& father) const
{
    return 0;
}

} // namespace octopus
