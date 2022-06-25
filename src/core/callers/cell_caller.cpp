// Copyright (c) 2015-2021 Daniel Cooke
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

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "basics/genomic_region.hpp"
#include "containers/probability_matrix.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/phylogeny.hpp"
#include "core/types/calls/cell_variant_call.hpp"
#include "core/types/calls/reference_call.hpp"
#include "core/models/genotype/uniform_genotype_prior_model.hpp"
#include "core/models/genotype/coalescent_genotype_prior_model.hpp"
#include "core/models/genotype/coalescent_population_prior_model.hpp"
#include "core/models/mutation/denovo_model.hpp"
#include "core/models/genotype/single_cell_prior_model.hpp"
#include "utils/maths.hpp"
#include "logging/logging.hpp"

namespace octopus {

using GenotypeBlock = MappableBlock<Genotype<IndexedHaplotype<>>>;

CellCaller::CellCaller(Caller::Components&& components,
                       Caller::Parameters general_parameters,
                       Parameters specific_parameters)
: Caller {std::move(components), std::move(general_parameters)}
, parameters_ {std::move(specific_parameters)}
{
    parameters_.max_copy_loss = std::min(parameters_.max_copy_loss, parameters_.ploidy - 1);
    std::sort(std::begin(parameters_.normal_samples), std::end(parameters_.normal_samples));
}

std::string CellCaller::do_name() const
{
    return "cell";
}

CellCaller::CallTypeSet CellCaller::do_call_types() const
{
    return {std::type_index(typeid(CellVariantCall))};
}

unsigned CellCaller::do_min_callable_ploidy() const
{
    return parameters_.ploidy;
}

unsigned CellCaller::do_max_callable_ploidy() const
{
    return parameters_.ploidy;
}

std::size_t CellCaller::do_remove_duplicates(HaplotypeBlock& haplotypes) const
{
    if (parameters_.deduplicate_haplotypes_with_prior_model) {
        if (haplotypes.size() < 2) return 0;
        CoalescentModel::Parameters model_params {};
        if (parameters_.prior_model_params) model_params = *parameters_.prior_model_params;
        Haplotype reference {mapped_region(haplotypes), reference_.get()};
        CoalescentModel model {std::move(reference), model_params, haplotypes.size(), CoalescentModel::CachingStrategy::none};
        const CoalescentProbabilityGreater cmp {std::move(model)};
        return octopus::remove_duplicates(haplotypes, cmp);
    } else {
        return Caller::do_remove_duplicates(haplotypes);
    }
}

// CellCaller::Latents public methods

CellCaller::Latents::Latents(const CellCaller& caller,
                             IndexedHaplotypeBlock haplotypes,
                             GenotypeBlock genotypes,
                             std::vector<model::SingleCellModel::Inferences> inferences)
: caller_ {caller}
, haplotypes_ {std::move(haplotypes)}
, genotypes_ {std::move(genotypes)}
, phylogeny_inferences_ {std::move(inferences)}
{
    phylogeny_posteriors_.resize(phylogeny_inferences_.size());
    std::transform(std::cbegin(phylogeny_inferences_), std::cend(phylogeny_inferences_), std::begin(phylogeny_posteriors_),
                   [] (const auto& inferences) { return inferences.log_evidence; });
    maths::normalise_exp(phylogeny_posteriors_);
    auto map_phylogeny_posteriors_itr = std::max_element(std::cbegin(phylogeny_posteriors_), std::cend(phylogeny_posteriors_));
    map_phylogeny_idx_ = std::distance(std::cbegin(phylogeny_posteriors_), map_phylogeny_posteriors_itr);
    phylogeny_size_posteriors_.resize(phylogeny_inferences_.back().phylogeny.size() + 1);
    for (std::size_t idx {0}; idx < phylogeny_inferences_.size(); ++idx) {
        phylogeny_size_posteriors_[phylogeny_inferences_[idx].phylogeny.size()] += phylogeny_posteriors_[idx];
    }
}

namespace {

using InverseGenotypeTable = std::vector<std::vector<std::size_t>>;

auto make_inverse_genotype_table(const MappableBlock<Genotype<IndexedHaplotype<>>>& genotypes,
                                 const std::size_t num_haplotypes)
{
    InverseGenotypeTable result(num_haplotypes);
    for (auto& indices : result) indices.reserve(genotypes.size() / num_haplotypes);
    for (std::size_t genotype_idx {0}; genotype_idx < genotypes.size(); ++genotype_idx) {
        for (const auto& haplotype : genotypes[genotype_idx]) {
            result[index_of(haplotype)].push_back(genotype_idx);
        }
    }
    for (auto& indices : result) {
        std::sort(std::begin(indices), std::end(indices));
        indices.erase(std::unique(std::begin(indices), std::end(indices)), std::end(indices));
        indices.shrink_to_fit();
    }
    return result;
}

using GenotypeMarginalPosteriorVector = std::vector<double>;
using GenotypeMarginalPosteriorMatrix = std::vector<GenotypeMarginalPosteriorVector>;

auto calculate_haplotype_posteriors(const MappableBlock<IndexedHaplotype<>>& haplotypes,
                                    const MappableBlock<Genotype<IndexedHaplotype<>>>& genotypes,
                                    const GenotypeMarginalPosteriorMatrix& genotype_posteriors)
{
    const auto inverse_genotypes = make_inverse_genotype_table(genotypes, haplotypes.size());
    std::vector<double> result(haplotypes.size());
    auto itr = std::cbegin(inverse_genotypes);
    std::vector<std::size_t> genotype_indices(genotypes.size());
    std::iota(std::begin(genotype_indices), std::end(genotype_indices), 0);
    // noncontaining genotypes are genotypes that do not contain a particular haplotype.
    const auto num_noncontaining_genotypes = genotypes.size() - itr->size();
    std::vector<std::size_t> noncontaining_genotype_indices(num_noncontaining_genotypes);
    for (const auto& haplotype : haplotypes) {
        std::set_difference(std::cbegin(genotype_indices), std::cend(genotype_indices),
                            std::cbegin(*itr), std::cend(*itr),
                            std::begin(noncontaining_genotype_indices));
        double prob_not_observed {1};
        for (const auto& sample_genotype_posteriors : genotype_posteriors) {
            prob_not_observed *= std::accumulate(std::cbegin(noncontaining_genotype_indices), std::cend(noncontaining_genotype_indices),
                                                 0.0, [&] (const auto total, const auto i) {
                return total + sample_genotype_posteriors[i];
            });
        }
        result[index_of(haplotype)] = 1.0 - prob_not_observed;
        ++itr;
    }
    return result;
}

} // namespace

std::shared_ptr<CellCaller::Latents::HaplotypeProbabilityMap>
CellCaller::Latents::haplotype_posteriors() const noexcept
{
    if (haplotype_posteriors_ == nullptr) {
        if (genotype_posteriors_ == nullptr) {
            this->genotype_posteriors(); // init genotype_posteriors_
        }
        auto haplotype_posteriors = calculate_haplotype_posteriors(haplotypes_, genotypes_, genotype_posteriors_array_);
        haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>();
        haplotype_posteriors_->reserve(haplotypes_.size());
        for (const auto& haplotype : haplotypes_) {
            haplotype_posteriors_->emplace(haplotype, haplotype_posteriors[index_of(haplotype)]);
        }
    }
    return haplotype_posteriors_;
}

namespace {

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
    auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
    auto zip_end   = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
    return boost::make_iterator_range(zip_begin, zip_end);
}

} // namespace

std::shared_ptr<CellCaller::Latents::GenotypeProbabilityMap>
CellCaller::Latents::genotype_posteriors() const noexcept
{
    if (genotype_posteriors_ == nullptr) {
        genotype_posteriors_array_.resize(caller_.samples_.size());
        for (std::size_t sample_idx {0}; sample_idx < caller_.samples_.size(); ++sample_idx) {
            std::vector<double> marginal_genotype_posteriors(genotypes_.size());
            for (const auto& p : zip(phylogeny_inferences_, phylogeny_posteriors_)) {
                const auto& phylogeny = p.get<0>().phylogeny;
                for (unsigned t {0}; t < phylogeny.size(); ++t) {
                    const auto& group = phylogeny.group(t).value;
                    for (std::size_t genotype_idx {0}; genotype_idx < group.genotype_posteriors.size(); ++genotype_idx) {
                        marginal_genotype_posteriors[genotype_idx] += p.get<1>()
                                * group.sample_attachment_posteriors[sample_idx]
                                * group.genotype_posteriors[genotype_idx];
                    }
                }
            }
            genotype_posteriors_array_[sample_idx] = std::move(marginal_genotype_posteriors);
        }
        genotype_posteriors_ = std::make_unique<GenotypeProbabilityMap>(std::cbegin(genotypes_), std::cend(genotypes_));
        for (std::size_t sample_idx {0}; sample_idx < caller_.samples_.size(); ++sample_idx) {
            insert_sample(caller_.samples_[sample_idx], genotype_posteriors_array_[sample_idx], *genotype_posteriors_);
        }
    }
    return genotype_posteriors_;
}

// CellCaller::Latents private methods

Phylogeny<std::size_t, Genotype<IndexedHaplotype<>>>
compute_map_phylogeny(const model::SingleCellModel::Inferences::InferedPhylogeny& phylogeny, const GenotypeBlock& genotypes)
{
    const auto get_map_genotype = [&genotypes] (const model::SingleCellModel::Inferences::GroupInferences& group) -> Genotype<IndexedHaplotype<>> {
        auto map_itr = std::max_element(std::cbegin(group.genotype_posteriors), std::cend(group.genotype_posteriors));
        auto map_idx = static_cast<std::size_t>(std::distance(std::cbegin(group.genotype_posteriors), map_itr));
        return genotypes[map_idx];
    };
    return phylogeny.transform(get_map_genotype);
}

template <typename S>
void log(const model::SingleCellModel::Inferences& inferences,
         const std::vector<SampleName>& samples,
         const GenotypeBlock& genotypes,
         S&& logger)
{
    std::ostringstream ss {};
    inferences.phylogeny.serialise(ss, [] (std::ostream& os, const auto& group) { os << group.id; });
    logger << "Phylogeny: " << ss.str() << '\n';
    std::vector<std::size_t> map_genotypes {};
    map_genotypes.reserve(inferences.phylogeny.size());
    std::vector<std::pair<std::size_t, double>> map_sample_assignments(samples.size());
    for (std::size_t group_id {0}; group_id < inferences.phylogeny.size(); ++group_id) {
        const auto& group = inferences.phylogeny.group(group_id).value;
        auto map_itr = std::max_element(std::cbegin(group.genotype_posteriors), std::cend(group.genotype_posteriors));
        auto map_idx = static_cast<std::size_t>(std::distance(std::cbegin(group.genotype_posteriors), map_itr));
        map_genotypes.push_back(map_idx);
        for (std::size_t sample_idx {0}; sample_idx < samples.size(); ++sample_idx) {
            if (group.sample_attachment_posteriors[sample_idx] > map_sample_assignments[sample_idx].second) {
                map_sample_assignments[sample_idx].first = group_id;
                map_sample_assignments[sample_idx].second = group.sample_attachment_posteriors[sample_idx];
            }
        }
    }
    logger << "MAP genotypes: " << '\n';
    for (std::size_t group_id {0}; group_id < map_genotypes.size(); ++group_id) {
        logger << group_id << ": "; debug::print_variant_alleles(logger, genotypes[map_genotypes[group_id]]); logger << '\n';
    }
    logger << "Sample MAP assignments:" << '\n';
    for (std::size_t sample_idx {0}; sample_idx < samples.size(); ++sample_idx) {
        logger << samples[sample_idx] << ": " << map_sample_assignments[sample_idx].first
               << " (" << map_sample_assignments[sample_idx].second << ")\n";
    }
    logger << "Evidence: " << inferences.log_evidence << '\n';
}

void log(const model::SingleCellModel::Inferences& inferences,
         const std::vector<SampleName>& samples,
         const GenotypeBlock& genotypes,
         boost::optional<logging::DebugLogger>& logger)
{
    if (logger) {
        log(inferences, samples, genotypes, stream(*logger));
    }
}

struct SingleCellModelInferencesEvidenceLess
{
    bool operator()(const model::SingleCellModel::Inferences& lhs, const model::SingleCellModel::Inferences& rhs) const noexcept
    {
        return lhs.log_evidence < rhs.log_evidence;
    }
};

std::vector<model::SingleCellPriorModel::CellPhylogeny>
propose_next_phylogenies(const std::vector<std::vector<model::SingleCellModel::Inferences>>& prev_inferences)
{
    using CellPhylogeny = model::SingleCellPriorModel::CellPhylogeny;
    std::vector<CellPhylogeny> result {};
    if (prev_inferences.empty()) {
        result = {CellPhylogeny {CellPhylogeny::Group {0}}};
    } else if (prev_inferences.size() == 1) {
        CellPhylogeny two_group_phylogeny {{0}};
        two_group_phylogeny.add_descendant({1}, 0);
        result = {std::move(two_group_phylogeny)};
    } else if (prev_inferences.size() == 2) {
        CellPhylogeny linear_three_group_phylogeny{{0}};
        linear_three_group_phylogeny.add_descendant({1}, 0);
        linear_three_group_phylogeny.add_descendant({2}, 1);
        CellPhylogeny forking_three_group_phylogeny{{0}};
        forking_three_group_phylogeny.add_descendant({1}, 0);
        forking_three_group_phylogeny.add_descendant({2}, 0);
        result = {std::move(linear_three_group_phylogeny), std::move(forking_three_group_phylogeny)};
    } else {
        using CellInferredPhylogeny = model::SingleCellModel::Inferences::InferedPhylogeny;
        const auto best_prev_inferences_itr = std::max_element(std::cbegin(prev_inferences.back()), std::cend(prev_inferences.back()), SingleCellModelInferencesEvidenceLess {});
        const CellInferredPhylogeny best_prev_phylogeny {best_prev_inferences_itr->phylogeny};
        const CellPhylogeny template_phylogeny {best_prev_phylogeny.transform([] (const auto&) -> CellPhylogeny::ValueType { return {}; })};
        result.reserve(best_prev_phylogeny.size());
        const auto prev_groups = best_prev_phylogeny.groups();
        for (const auto& group : prev_groups) {
            if (template_phylogeny.num_descendants(group.id) < 2) {
                auto extended_phylogeny = template_phylogeny;
                extended_phylogeny.add_descendant({prev_groups.size()}, group.id);
                result.push_back(std::move(extended_phylogeny));
            }
        }
    }
    return result;
}

namespace {

bool includes(const std::vector<SampleName>& samples, const SampleName& sample)
{
    return std::binary_search(std::cbegin(samples), std::cend(samples), sample);
}

} // namespace

std::unique_ptr<CellCaller::Caller::Latents>
CellCaller::infer_latents(const HaplotypeBlock& haplotypes, 
                          const HaplotypeLikelihoodArray& haplotype_likelihoods,
                          OptionalThreadPool workers) const
{
    const auto indexed_haplotypes = index(haplotypes);
    auto genotypes = generate_all_genotypes(indexed_haplotypes, parameters_.ploidy);
    
    if (debug_log_) stream(*debug_log_) << "There are " << genotypes.size() << " candidate genotypes";
    
    GenotypeBlock copy_change_genotypes {};
    std::size_t default_ploidy_idx {0};
    bool copy_number_change_detection_enabled {false};
    if (parameters_.max_copy_loss > 0 || parameters_.max_copy_gain > 0) {
        copy_number_change_detection_enabled = true;
        for (unsigned loss {1}; loss <= parameters_.max_copy_loss; ++loss) {
            auto copy_loss_genotypes = generate_all_genotypes(indexed_haplotypes, parameters_.ploidy - loss);
            default_ploidy_idx += copy_loss_genotypes.size();
            utils::append(std::move(copy_loss_genotypes), copy_change_genotypes);
        }
        utils::append(genotypes, copy_change_genotypes);
        for (unsigned gain {1}; gain <= parameters_.max_copy_gain; ++gain) {
            auto copy_gain_genotypes = generate_all_genotypes(indexed_haplotypes, parameters_.ploidy + gain);
            utils::append(std::move(copy_gain_genotypes), copy_change_genotypes);
        }
    }
    
    const auto genotype_prior_model = make_prior_model(haplotypes);
    genotype_prior_model->prime(haplotypes);
    DeNovoModel mutation_model {parameters_.mutation_model_parameters};
    mutation_model.prime(haplotypes);
    model::SingleCellPriorModel::Parameters cell_prior_params {};
    cell_prior_params.copy_number_prior = parameters_.somatic_cnv_prior;
    model::SingleCellModel::Parameters model_parameters {};
    model_parameters.dropout_concentration = parameters_.dropout_concentration;
    if (!parameters_.sample_dropout_concentrations.empty()) {
        model_parameters.sample_dropout_concentrations.resize(samples_.size(), parameters_.dropout_concentration);
        for (std::size_t s {0}; s < samples_.size(); ++s) {
            if (parameters_.sample_dropout_concentrations.count(samples_[s]) == 1) {
                model_parameters.sample_dropout_concentrations[s] = parameters_.sample_dropout_concentrations.at(samples_[s]);
            }
        }
    }
    model_parameters.group_concentration = parameters_.clone_concentration;
    model::SingleCellModel::AlgorithmParameters config {};
    config.max_genotype_combinations = *parameters_.max_genotype_combinations;
    config.execution_policy = this->exucution_policy();
    if (parameters_.max_vb_seeds) config.max_seeds = *parameters_.max_vb_seeds;
    CoalescentPopulationPriorModel population_prior_model {{Haplotype {mapped_region(haplotypes), reference_}, {}}};
    population_prior_model.prime(haplotypes);
    using SingleCellModelInferences = model::SingleCellModel::Inferences;
    std::vector<std::vector<SingleCellModelInferences>> inferences {};
    double max_log_evidence {};
    bool copy_change_predicted {false};
    const auto max_clones = std::min({parameters_.max_clones, static_cast<unsigned>(genotypes.size()), static_cast<unsigned>(samples_.size())});
    
    for (unsigned clones {1}; clones <= max_clones; ++clones) {
        auto phylogenies = propose_next_phylogenies(inferences);
        if (!phylogenies.empty()) {
            std::vector<SingleCellModelInferences> clone_inferences {};
            clone_inferences.reserve(phylogenies.size());
            for (auto& phylogeny : phylogenies) {
                assert(phylogeny.size() == clones);
                if (parameters_.normal_samples.empty()) {
                    model_parameters.group_priors = boost::none;
                } else {
                    std::vector<double> normal_group_priors(phylogeny.size(), parameters_.normal_not_founder_prior / (phylogeny.size() - 1));
                    normal_group_priors[0] = 1 - parameters_.normal_not_founder_prior;
                    model_parameters.group_priors = model::SingleCellModel::Parameters::GroupOptionalPriorArray {};
                    model_parameters.group_priors->reserve(samples_.size());
                    for (const auto& sample : samples_) {
                        if (includes(parameters_.normal_samples, sample)) {
                            model_parameters.group_priors->push_back(normal_group_priors);
                        } else {
                            model_parameters.group_priors->push_back(boost::none);
                        }
                    }
                }
                model::SingleCellPriorModel phylogeny_prior_model {std::move(phylogeny), *genotype_prior_model, mutation_model, cell_prior_params};
                const model::SingleCellModel phylogeny_model {samples_, std::move(phylogeny_prior_model), model_parameters, config, population_prior_model};
                auto phylogeny_inferences = phylogeny_model.evaluate(genotypes, haplotype_likelihoods);
                log(phylogeny_inferences, samples_, genotypes, debug_log_);
                
                if (clones > 1 && copy_number_change_detection_enabled) {
                    std::vector<unsigned> phylogeny_ploidy_assignments((1 + parameters_.max_copy_loss + parameters_.max_copy_gain) * (clones - 1));
                    auto assignment_itr = std::begin(phylogeny_ploidy_assignments);
                    for (auto ploidy = parameters_.ploidy - parameters_.max_copy_loss; ploidy <= parameters_.ploidy + parameters_.max_copy_gain; ++ploidy) {
                        assignment_itr = std::fill_n(assignment_itr, clones - 1, ploidy);
                    }
                    std::unordered_map<std::size_t, unsigned> phylogeny_ploidies {};
                    phylogeny_ploidies.reserve(clones);
                    bool can_ignore_future_copy_changes {false};
                    using PloidyLabeledPhylogeny = Phylogeny<std::size_t, unsigned>;
                    std::deque<PloidyLabeledPhylogeny> evaluated_phylogenies {};
                    do {
                        if (phylogeny_ploidy_assignments[0] == parameters_.ploidy) {
                            for (std::size_t id {0}; id < clones; ++id) {
                                phylogeny_ploidies[id] = phylogeny_ploidy_assignments[id];
                            }
                            auto labeled_phylogeny = phylogeny_model.prior_model().phylogeny().transform([&] (const auto&) { return 0u; });
                            for (const auto& p : phylogeny_ploidies) {
                                labeled_phylogeny.group(p.first).value = p.second;
                            }
                            const auto is_isomorphic = [&] (const auto& other) { return labeled_phylogeny.is_isomorphism(other); };
                            if (std::none_of(std::cbegin(evaluated_phylogenies), std::cend(evaluated_phylogenies), is_isomorphic)) {
                                try {
                                    auto phylogeny_copy_inferences = phylogeny_model.evaluate(phylogeny_ploidies, copy_change_genotypes, haplotype_likelihoods);
                                    log(phylogeny_copy_inferences, samples_, copy_change_genotypes, debug_log_);
                                    if (phylogeny_copy_inferences.log_evidence > phylogeny_inferences.log_evidence) {
                                        phylogeny_inferences = std::move(phylogeny_copy_inferences);
                                        copy_change_predicted = true;
                                        can_ignore_future_copy_changes = false;
                                    } else if (phylogeny_copy_inferences.log_evidence > max_log_evidence) {
                                        can_ignore_future_copy_changes = false;
                                    }
                                    phylogeny_ploidies.clear();
                                    evaluated_phylogenies.push_back(std::move(labeled_phylogeny));
                                } catch (const model::SingleCellModel::NoViableGenotypeCombinationsError&) {
                                    can_ignore_future_copy_changes = true;
                                    break;
                                }
                            }
                        }
                    } while (std::next_permutation(std::begin(phylogeny_ploidy_assignments), std::end(phylogeny_ploidy_assignments)));
                    if (can_ignore_future_copy_changes) copy_number_change_detection_enabled = false;
                }
                clone_inferences.push_back(std::move(phylogeny_inferences));
            }
            if (clones == 1) {
                max_log_evidence = clone_inferences.front().log_evidence;
            } else {
                const static auto evidence_less = [] (const auto& lhs, const auto& rhs) { return lhs.log_evidence < rhs.log_evidence; };
                const auto max_clone_log_evidence = std::max_element(std::cbegin(clone_inferences), std::cend(clone_inferences), evidence_less)->log_evidence;
                if (max_clone_log_evidence < max_log_evidence) {
                    inferences.push_back(std::move(clone_inferences));
                    break;
                } else {
                    max_log_evidence = max_clone_log_evidence;
                }
            }
            inferences.push_back(std::move(clone_inferences));
        }
    }
    std::vector<SingleCellModelInferences> flat_inferences {};
    for (auto& clone_inferences : inferences) {
        utils::append(std::move(clone_inferences), flat_inferences);
    }
    if (copy_change_predicted) {
        for (SingleCellModelInferences& phylogeny_inferences : flat_inferences) {
            for (std::size_t id {0}; id < phylogeny_inferences.phylogeny.size(); ++id) {
                auto& genotype_posteriors = phylogeny_inferences.phylogeny.group(id).value.genotype_posteriors;
                if (genotype_posteriors.size() < copy_change_genotypes.size()) {
                    genotype_posteriors.resize(copy_change_genotypes.size());
                    std::rotate(std::rbegin(genotype_posteriors), std::next(std::rbegin(genotype_posteriors), default_ploidy_idx), std::rend(genotype_posteriors));
                }
            }
        }
        genotypes = std::move(copy_change_genotypes);
    }
    return std::make_unique<Latents>(*this, std::move(indexed_haplotypes), std::move(genotypes), std::move(flat_inferences));
}

std::vector<std::unique_ptr<octopus::VariantCall>>
CellCaller::call_variants(const std::vector<Variant>& candidates, const Caller::Latents& latents, OptionalThreadPool workers) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}

namespace {

using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<IndexedHaplotype<>>>::InnerMap;
using PopulationGenotypeProbabilityMap = ProbabilityMatrix<Genotype<IndexedHaplotype<>>>;

using VariantReference = std::reference_wrapper<const Variant>;
using VariantPosteriorVector = std::vector<std::pair<VariantReference, std::vector<Phred<double>>>>;

struct VariantCall : Mappable<VariantCall>
{
    VariantCall() = delete;
    VariantCall(const std::pair<VariantReference, std::vector<Phred<double>>>& p)
    : variant {p.first}
    , posteriors {p.second}
    {}
    VariantCall(const Variant& variant, std::vector<Phred<double>> posterior)
    : variant {variant}
    , posteriors {posterior}
    {}
    
    const GenomicRegion& mapped_region() const noexcept
    {
        return octopus::mapped_region(variant.get());
    }
    
    VariantReference variant;
    std::vector<Phred<double>> posteriors;
};

using VariantCalls = std::vector<VariantCall>;

struct GenotypeCall
{
    Genotype<Allele> genotype;
    Phred<double> posterior;
};

using GenotypeCalls = std::vector<std::vector<GenotypeCall>>;

// allele posterior calculations

using AlleleBools           = std::deque<bool>; // using std::deque because std::vector<bool> is evil
using GenotypePropertyBools = std::vector<AlleleBools>;

auto marginalise(const GenotypeProbabilityMap& genotype_posteriors,
                 const AlleleBools& contained_alleles)
{
    auto p = std::inner_product(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                std::cbegin(contained_alleles), 0.0, std::plus<> {},
                                [] (const auto& p, const bool is_contained) {
                                    return is_contained ? 0.0 : p.second;
                                });
    return probability_false_to_phred(p);
}

auto compute_sample_allele_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
                                      const GenotypePropertyBools& contained_alleles)
{
    std::vector<Phred<double>> result {};
    result.reserve(contained_alleles.size());
    for (const auto& allele : contained_alleles) {
        result.emplace_back(marginalise(genotype_posteriors, allele));
    }
    return result;
}

auto get_contained_alleles(const PopulationGenotypeProbabilityMap& genotype_posteriors,
                           const std::vector<Allele>& alleles)
{
    const auto num_genotypes = genotype_posteriors.size2();
    GenotypePropertyBools result {};
    if (num_genotypes == 0 || genotype_posteriors.empty1() || alleles.empty()) {
        return result;
    }
    result.reserve(alleles.size());
    const auto& test_sample   = genotype_posteriors.begin()->first;
    const auto genotype_begin = genotype_posteriors.begin(test_sample);
    const auto genotype_end   = genotype_posteriors.end(test_sample);
    for (const auto& allele : alleles) {
        result.emplace_back(num_genotypes);
        transform_contains(genotype_begin, genotype_end, allele, std::begin(result.back()),
                           [] (const auto& p, bool contain) { return contain; },
                           [] (const auto& p) { return p.first; });
    }
    return result;
}

using AllelePosteriorMatrix = std::vector<std::vector<Phred<double>>>;

auto compute_posteriors(const std::vector<SampleName>& samples,
                        const std::vector<Allele>& alleles,
                        const PopulationGenotypeProbabilityMap& genotype_posteriors)
{
    const auto contained_alleles = get_contained_alleles(genotype_posteriors, alleles);
    AllelePosteriorMatrix result {};
    result.reserve(genotype_posteriors.size1());
    for (const auto& sample : samples) {
        result.emplace_back(compute_sample_allele_posteriors(genotype_posteriors[sample], contained_alleles));
    }
    return result;
}

auto extract_ref_alleles(const std::vector<Variant>& variants)
{
    std::vector<Allele> result {};
    result.reserve(variants.size());
    std::transform(std::cbegin(variants), std::cend(variants), std::back_inserter(result),
                   [] (const auto& variant) { return variant.ref_allele(); });
    return result;
}

auto extract_alt_alleles(const std::vector<Variant>& variants)
{
    std::vector<Allele> result {};
    result.reserve(variants.size());
    std::transform(std::cbegin(variants), std::cend(variants), std::back_inserter(result),
                   [] (const auto& variant) { return variant.alt_allele(); });
    return result;
}

auto compute_posteriors(const std::vector<SampleName>& samples,
                        const std::vector<Variant>& variants,
                        const PopulationGenotypeProbabilityMap& genotype_posteriors)
{
    const auto allele_posteriors = compute_posteriors(samples, extract_alt_alleles(variants), genotype_posteriors);
    VariantPosteriorVector result {};
    result.reserve(variants.size());
    for (std::size_t i {0}; i < variants.size(); ++i) {
        std::vector<Phred<double>> sample_posteriors(samples.size());
        std::transform(std::cbegin(allele_posteriors), std::cend(allele_posteriors), std::begin(sample_posteriors),
                       [i] (const auto& ps) { return ps[i]; });
        result.emplace_back(variants[i], std::move(sample_posteriors));
    }
    return result;
}

// haplotype genotype calling

auto call_genotype(const PopulationGenotypeProbabilityMap::InnerMap& genotype_posteriors)
{
    return std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                            [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; })->first;
}

auto call_genotypes(const std::vector<SampleName>& samples, const PopulationGenotypeProbabilityMap& genotype_posteriors)
{
    GenotypeBlock result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.push_back(call_genotype(genotype_posteriors[sample]));
    }
    return result;
}

// variant calling

bool has_above(const std::vector<Phred<double>>& posteriors, const Phred<double> min_posterior)
{
    return std::any_of(std::cbegin(posteriors), std::cend(posteriors), [=] (auto p) { return p >= min_posterior; });
}

bool contains_alt(const Genotype<IndexedHaplotype<>>& genotype_call, const VariantReference& candidate)
{
    return includes(genotype_call, candidate.get().alt_allele());
}

bool contains_alt(const GenotypeBlock& genotype_calls, const VariantReference& candidate)
{
    return std::any_of(std::cbegin(genotype_calls), std::cend(genotype_calls),
                       [&] (const auto& genotype) { return contains_alt(genotype, candidate); });
}

VariantCalls call_candidates(const VariantPosteriorVector& candidate_posteriors,
                             const GenotypeBlock& genotype_calls,
                             const Phred<double> min_posterior)
{
    VariantCalls result {};
    result.reserve(candidate_posteriors.size());
    std::copy_if(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors),
                 std::back_inserter(result),
                 [&genotype_calls, min_posterior] (const auto& p) {
                     return has_above(p.second, min_posterior) && contains_alt(genotype_calls, p.first);
                 });
    return result;
}

// allele genotype calling

auto marginalise(const Genotype<Allele>& genotype, const GenotypeProbabilityMap& genotype_posteriors)
{
    auto p = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), 0.0,
                             [&genotype] (const double curr, const auto& p) {
                                 return curr + (contains(p.first, genotype) ? 0.0 : p.second);
                             });
    return probability_false_to_phred(p);
}

auto call_genotypes(const std::vector<SampleName>& samples,
                    const GenotypeBlock& genotype_calls,
                    const PopulationGenotypeProbabilityMap& genotype_posteriors,
                    const std::vector<GenomicRegion>& variant_regions)
{
    GenotypeCalls result {};
    result.reserve(variant_regions.size());
    for (const auto& region : variant_regions) {
        std::vector<GenotypeCall> region_calls {};
        region_calls.reserve(samples.size());
        for (std::size_t s {0}; s < samples.size(); ++s) {
            auto genotype_chunk = copy<Allele>(genotype_calls[s], region);
            const auto posterior = marginalise(genotype_chunk, genotype_posteriors[samples[s]]);
            region_calls.push_back({std::move(genotype_chunk), posterior});
        }
        result.push_back(std::move(region_calls));
    }
    return result;
}

// output

using PhylogenyInferenceSummary = CellVariantCall::PhylogenyInferenceSummary;

octopus::VariantCall::GenotypeCall convert(GenotypeCall&& call)
{
    return octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
}

std::unique_ptr<octopus::VariantCall>
transform_call(const std::vector<SampleName>& samples,
               VariantCall&& variant_call,
               std::vector<GenotypeCall>&& sample_genotype_calls,
               PhylogenyInferenceSummary phylogeny_summary)
{
    std::vector<std::pair<SampleName, Call::GenotypeCall>> tmp {};
    tmp.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples),
                   std::make_move_iterator(std::begin(sample_genotype_calls)),
                   std::back_inserter(tmp),
                   [] (const auto& sample, auto&& genotype) {
                       return std::make_pair(sample, convert(std::move(genotype)));
                   });
    auto quality = *std::max_element(std::cbegin(variant_call.posteriors), std::cend(variant_call.posteriors));
    return std::make_unique<CellVariantCall>(variant_call.variant.get(), std::move(tmp), quality, std::move(phylogeny_summary));
}

auto transform_calls(const std::vector<SampleName>& samples,
                     VariantCalls&& variant_calls,
                     GenotypeCalls&& genotype_calls,
                     PhylogenyInferenceSummary phylogeny_summary)
{
    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
    result.reserve(variant_calls.size());
    std::transform(std::make_move_iterator(std::begin(variant_calls)), std::make_move_iterator(std::end(variant_calls)),
                   std::make_move_iterator(std::begin(genotype_calls)), std::back_inserter(result),
                   [&] (auto&& variant_call, auto&& genotype_call) {
                       return transform_call(samples, std::move(variant_call), std::move(genotype_call), phylogeny_summary);
                   });
    return result;
}

} // namespace

std::vector<std::unique_ptr<octopus::VariantCall>>
CellCaller::call_variants(const std::vector<Variant>& candidates, const Latents& latents) const
{
    const auto& genotype_posteriors = *(latents.genotype_posteriors());
    const auto sample_candidate_posteriors = compute_posteriors(samples_, candidates, genotype_posteriors);
    const auto genotype_calls = call_genotypes(samples_, genotype_posteriors);
    auto variant_calls = call_candidates(sample_candidate_posteriors, genotype_calls, parameters_.min_variant_posterior);
    const auto called_regions = extract_regions(variant_calls);
    auto allele_genotype_calls = call_genotypes(samples_, genotype_calls, genotype_posteriors, called_regions);
    PhylogenyInferenceSummary summary {};
    using SummaryPhylogeny = decltype(summary.map);
    const auto& map_phylogeny = latents.phylogeny_inferences_[latents.map_phylogeny_idx_].phylogeny;
    summary.map = map_phylogeny.transform([] (const auto&) -> SummaryPhylogeny::ValueType { return {}; });
    const static auto prob_true_to_phred = [] (double p) { return probability_false_to_phred(std::max(1 - p, 0.0)); };
    summary.map_posterior = prob_true_to_phred(latents.phylogeny_posteriors_[latents.map_phylogeny_idx_]);
    summary.size_posteriors.resize(latents.phylogeny_size_posteriors_.size());
    summary.sample_node_posteriors.reserve(samples_.size());
    for (const auto& sample : samples_) {
        summary.sample_node_posteriors[sample].resize(map_phylogeny.size());
    }
    for (std::size_t id {0}; id < map_phylogeny.size(); ++id) {
        const auto& group = map_phylogeny.group(id);
        for (std::size_t sample_idx {0}; sample_idx < samples_.size(); ++sample_idx) {
            summary.sample_node_posteriors[samples_[sample_idx]][id] = prob_true_to_phred(group.value.sample_attachment_posteriors[sample_idx]);
        }
    }
    std::transform(std::cbegin(latents.phylogeny_size_posteriors_), std::cend(latents.phylogeny_size_posteriors_),
                   std::begin(summary.size_posteriors), prob_true_to_phred);
    return transform_calls(samples_, std::move(variant_calls), std::move(allele_genotype_calls), std::move(summary));
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

std::unique_ptr<GenotypePriorModel> CellCaller::make_prior_model(const HaplotypeBlock& haplotypes) const
{
    if (parameters_.prior_model_params) {
        return std::make_unique<CoalescentGenotypePriorModel>(CoalescentModel {
        Haplotype {mapped_region(haplotypes), reference_},
        *parameters_.prior_model_params, haplotypes.size(), CoalescentModel::CachingStrategy::address
        });
    } else {
        return std::make_unique<UniformGenotypePriorModel>();
    }
}

} // namespace octopus
