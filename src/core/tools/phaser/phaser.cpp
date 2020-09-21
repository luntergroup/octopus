// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "phaser.hpp"

#include <deque>
#include <map>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <cmath>
#include <utility>
#include <iostream>

#include <boost/functional/hash.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graphviz.hpp>

#if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-variable"
        #include <boost/graph/bron_kerbosch_all_cliques.hpp>
    #pragma clang diagnostic pop
#elif defined(__GNUG__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
        #include <boost/graph/bron_kerbosch_all_cliques.hpp>
    #pragma GCC diagnostic pop
#endif // defined (__clang__)

#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"
#include "utils/map_utils.hpp"

namespace octopus {

Phaser::Phaser(Config config) : config_ {config} {}

namespace {

using CompressedGenotype = Genotype<IndexedHaplotype<>>;
using SharedGenotype = Genotype<SharedHaplotype>;

template <typename Range>
auto minmax_ploidy(const Range& genotypes) noexcept
{
    assert(!genotypes.empty());
    const static auto ploidy_less = [] (const auto& lhs, const auto& rhs) { return lhs.ploidy() < rhs.ploidy(); };
    auto p = std::minmax_element(std::cbegin(genotypes), std::cend(genotypes), ploidy_less);
    return std::make_pair(p.first->ploidy(), p.second->ploidy());
}

double maximum_entropy(const std::size_t num_elements)
{
    return std::log2(num_elements);
}

auto min_phase_quality(double p)
{
    if (maths::almost_one(p)) {
        static const Phred<double> max_possible_score {Phred<double>::Probability {0.0}};
        return max_possible_score;
    } else {
        p = std::max(p, 0.5);
        const auto e = -(p * std::log2(p) + (1.0 - p) * std::log2(1.0 - p)) / maximum_entropy(2);
        return Phred<double> {Phred<double>::Probability {std::max(e, 0.0)}};
    }
}

auto min_phase_quality(const CompressedGenotype& called_genotype, const Phaser::SampleGenotypePosteriorMap& genotype_posteriors)
{
    return min_phase_quality(genotype_posteriors[called_genotype]);
}

Phaser::GenotypePosteriorMap marginalise_collapsed_genotypes(const Phaser::GenotypePosteriorMap& genotype_posteriors)
{
    auto genotypes = extract_key_refs(genotype_posteriors);
    const auto num_genotypes = genotypes.size();
    std::vector<CompressedGenotype> collapsed_genotypes {};
    collapsed_genotypes.reserve(num_genotypes);
    std::unordered_map<CompressedGenotype, std::size_t> collapsed_genotype_indices {};
    collapsed_genotype_indices.reserve(num_genotypes);
    std::vector<std::size_t> collapsed_genotype_index_table {};
    collapsed_genotype_index_table.reserve(num_genotypes);
    for (const auto& genotype : genotypes) {
        auto collapsed = collapse(genotype.get());
        const auto itr = collapsed_genotype_indices.find(collapsed);
        if (itr == std::cend(collapsed_genotype_indices)) {
            const auto index = collapsed_genotypes.size();
            collapsed_genotype_indices.emplace(collapsed, index);
            collapsed_genotypes.push_back(std::move(collapsed));
            collapsed_genotype_index_table.push_back(index);
        } else {
            collapsed_genotype_index_table.push_back(itr->second);
        }
    }
    genotypes.clear();
    genotypes.shrink_to_fit();
    collapsed_genotype_indices.clear();
    const auto num_collapsed_genotypes = collapsed_genotypes.size();
    Phaser::GenotypePosteriorMap result {std::make_move_iterator(std::begin(collapsed_genotypes)),
                                         std::make_move_iterator(std::end(collapsed_genotypes))};
    collapsed_genotypes.clear();
    collapsed_genotypes.shrink_to_fit();
    for (const auto& sample_posteriors : genotype_posteriors) {
        std::vector<double> posteriors(num_collapsed_genotypes);
        std::size_t genotype_idx {0};
        std::for_each(std::cbegin(sample_posteriors.second), std::cend(sample_posteriors.second), [&] (const auto& p) {
            posteriors[collapsed_genotype_index_table[genotype_idx]] += p.second;
            ++genotype_idx;
        });
        insert_sample(sample_posteriors.first, std::move(posteriors), result);
    }
    return result;
}

void collapse_each(Phaser::GenotypeCallMap& genotypes)
{
    for (auto& p : genotypes) collapse(p.second);
}

auto copy_unique(std::vector<GenomicRegion> regions)
{
    assert(std::is_sorted(std::cbegin(regions), std::cend(regions)));
    regions.erase(std::unique(std::begin(regions), std::end(regions)), std::end(regions));
    return regions;
}

} // namespace

Phaser::PhaseSetMap
Phaser::phase(const MappableBlock<Haplotype>& haplotypes,
              const GenotypePosteriorMap& genotype_posteriors,
              const std::vector<GenomicRegion>& variation_sites,
              boost::optional<GenotypeCallMap> genotype_calls) const
{
    assert(!haplotypes.empty());
    assert(!genotype_posteriors.empty1() && !genotype_posteriors.empty2());
    assert(std::is_sorted(std::cbegin(variation_sites), std::cend(variation_sites)));
    const auto unique_variation_sites = copy_unique(variation_sites);
    assert(!unique_variation_sites.empty() && unique_variation_sites.size() <= variation_sites.size());
    auto genotypes = extract_keys(genotype_posteriors);
    unsigned min_genotype_ploidy, max_genotype_ploidy; std::tie(min_genotype_ploidy, max_genotype_ploidy) = minmax_ploidy(genotypes);
    PhaseSetMap result {};
    const auto num_samples = genotype_posteriors.size1();
    result.reserve(num_samples);
    std::vector<std::size_t> max_phase_set(variation_sites.size());
    std::iota(std::begin(max_phase_set), std::end(max_phase_set), std::size_t {0});
    static const Phred<double> max_possible_quality {Phred<double>::Probability {0.0}};
    if (max_genotype_ploidy == 1 || unique_variation_sites.size() == 1) {
        for (const auto& p : genotype_posteriors) {
            const SampleName& sample {p.first};
            result[sample].push_back({max_phase_set, max_possible_quality});
        }
    } else {
        boost::optional<GenotypePosteriorMap> collapsed_genotype_posteriors {};
        for (const auto& p : genotype_posteriors) {
            const SampleName& sample {p.first};
            if (!collapsed_genotype_posteriors && genotype_calls && config_.max_phase_quality
             && min_phase_quality(genotype_calls->at(sample), p.second) >= *config_.max_phase_quality) {
                result[sample].push_back({max_phase_set, *config_.max_phase_quality});
            } else {
                if (!collapsed_genotype_posteriors && (max_genotype_ploidy > 2 || min_genotype_ploidy != max_genotype_ploidy)) {
                    collapsed_genotype_posteriors = marginalise_collapsed_genotypes(genotype_posteriors);
                    genotypes = extract_keys(*collapsed_genotype_posteriors);
                    std::tie(min_genotype_ploidy, max_genotype_ploidy) = minmax_ploidy(genotypes);
                    if (genotype_calls) collapse_each(*genotype_calls);
                }
                PhaseSetVector phase_sets {};
                if (collapsed_genotype_posteriors) {
                    const auto& collapsed_sample_genotype_posteriors = (*collapsed_genotype_posteriors)[sample];
                    if (genotype_calls && config_.max_phase_quality
                     && min_phase_quality(genotype_calls->at(sample), collapsed_sample_genotype_posteriors) >= *config_.max_phase_quality) {
                        result[sample].push_back({max_phase_set, *config_.max_phase_quality});
                    } else {
                        phase_sets = phase_sample(unique_variation_sites, genotypes, collapsed_sample_genotype_posteriors);
                    }
                } else {
                    phase_sets = phase_sample(unique_variation_sites, genotypes, p.second);
                }
                result.emplace(sample, std::move(phase_sets));
            }
        }
    }
    if (config_.max_phase_quality) {
        for (auto& p : result) {
            for (auto& phase_set : p.second) {
                phase_set.quality = std::min(phase_set.quality, *config_.max_phase_quality);
            }
        }
    }
    if (unique_variation_sites.size() < variation_sites.size()) {
        std::vector<std::vector<std::size_t>> unique_index_to_input_index_table {};
        for (auto& p : result) {
            assert(!p.second.empty());
            if (p.second.size() > 1 || p.second.front().site_indices.size() < variation_sites.size()) {
                // Need to convert unique variation site indices back to input indices
                unique_index_to_input_index_table.resize(unique_variation_sites.size());
                for (auto& phase_set : p.second) {
                    auto& unique_site_indices = phase_set.site_indices;
                    for (auto site_index_itr = std::begin(unique_site_indices);
                         site_index_itr != std::end(unique_site_indices);
                         ++site_index_itr) {
                        auto& site_indices = unique_index_to_input_index_table[*site_index_itr];
                        if (site_indices.empty()) {
                            // lazy evaluate this
                            const auto& target_site = unique_variation_sites[*site_index_itr];
                            const auto overlapped_sites = std::equal_range(std::cbegin(variation_sites), std::cend(variation_sites), target_site);
                            auto first_site_index = static_cast<std::size_t>(std::distance(std::cbegin(variation_sites), overlapped_sites.first));
                            site_indices.resize(std::distance(overlapped_sites.first, overlapped_sites.second));
                            std::iota(std::begin(site_indices), std::end(site_indices), first_site_index);
                        }
                        assert(!site_indices.empty());
                        assert(site_indices.front() >= *site_index_itr);
                        if (site_indices.size() == 1) {
                            *site_index_itr = site_indices.front();
                        } else {
                            *site_index_itr = site_indices.front();
                            ++site_index_itr;
                            site_index_itr = unique_site_indices.insert(site_index_itr,
                                                                      std::next(std::cbegin(site_indices)), std::cend(site_indices));
                            std::advance(site_index_itr, site_indices.size() - 2);
                        }
                    }
                }
            }
        }
    }
    return result;
}

namespace {

using GenotypeChunk = SharedGenotype;
using GenotypeChunkPosteriorMap = std::unordered_map<GenotypeChunk, double>;

struct GenotypeInfo
{
    using AlleleIndex = std::uint8_t;
    using AlleleIndexSet = std::vector<AlleleIndex>;
    AlleleIndexSet alleles;
};

using GenotypeInfoVector = std::vector<GenotypeInfo>;
using GenotypeInfoMatrix = std::vector<GenotypeInfoVector>;

using AlleleVector = std::vector<Allele>;

GenotypeInfo compute_genotype_info(const Genotype<Allele>& genotype, AlleleVector& alleles)
{
    GenotypeInfo::AlleleIndexSet indices(genotype.ploidy());
    const auto get_allele_index = [&] (const Allele& allele) -> GenotypeInfo::AlleleIndex {
        const auto allele_itr = std::find(std::cbegin(alleles), std::cend(alleles), allele);
        if (allele_itr != std::cend(alleles)) {
            return std::distance(std::cbegin(alleles), allele_itr);
        }
        alleles.push_back(allele);
        return alleles.size() - 1;
    };
    std::transform(std::cbegin(genotype), std::cend(genotype), std::begin(indices), get_allele_index);
    std::sort(std::begin(indices), std::end(indices));
    indices.erase(std::unique(std::begin(indices), std::end(indices)), std::end(indices));
    return {std::move(indices)};
}

GenotypeInfoVector compute_genotype_info(const std::vector<CompressedGenotype>& genotypes, const GenomicRegion& region)
{
    GenotypeInfoVector result(genotypes.size());
    AlleleVector alleles {};
    alleles.reserve(5);
    transform_each(std::cbegin(genotypes), std::cend(genotypes),
                  [&] (const auto& element) { return copy<Allele>(element, region); }, 
                  [&] (const auto& genotype) { return compute_genotype_info(genotype, alleles); },
                  std::begin(result));
    return result;
}

auto compute_genotype_info(const std::vector<CompressedGenotype>& genotypes, const std::vector<GenomicRegion>& regions)
{
    GenotypeInfoMatrix result {};
    result.reserve(regions.size());
    for (const auto& region : regions) {
        result.push_back(compute_genotype_info(genotypes, region));
    }
    return result;
}

bool is_heterozygous(std::size_t genotype_index, const GenotypeInfoVector& info) noexcept
{
    return info[genotype_index].alleles.size() > 1;
}

bool is_very_likely_homozygous(const GenomicRegion& region, 
                               const Phaser::SampleGenotypePosteriorMap& genotype_posteriors, 
                               const GenotypeInfoVector& info)
{
    assert(!genotype_posteriors.empty());
    const static auto posterior_less = [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; };
    const auto map_posterior_itr = std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), posterior_less);
    const auto map_posterior = map_posterior_itr->second;
    return map_posterior > 0.9999 && !is_heterozygous(std::distance(std::cbegin(genotype_posteriors), map_posterior_itr), info);
}

auto copy_each_helper(const std::vector<CompressedGenotype>& genotypes,
                      const std::vector<GenomicRegion>& sites,
                      const std::size_t lhs, const std::size_t rhs)
{
    if (!are_adjacent(sites[lhs], sites[rhs])) {
        const std::vector<GenomicRegion> sites_tmp {sites[lhs], sites[rhs]};
        return copy_each<GenotypeChunk::ElementType>(genotypes, sites_tmp);
    } else {
        assert(is_before(sites[lhs], sites[rhs]));
        const auto site = closed_region(sites[lhs], sites[rhs]);
        return copy_each<GenotypeChunk::ElementType>(genotypes, site);
    }
}

auto
compute_chunk_set_posteriors(const std::vector<CompressedGenotype>& genotypes,
                             const std::vector<GenomicRegion>& sites,
                             const std::size_t lhs, const std::size_t rhs,
                             const Phaser::SampleGenotypePosteriorMap& genotype_posteriors,
                             const GenotypeInfoMatrix& info)
{
    auto chunks = copy_each_helper(genotypes, sites, lhs, rhs);
    using AlleleIndexSetRef = std::reference_wrapper<const GenotypeInfo::AlleleIndexSet>;
    using AlleleIndexSetRefPair = std::pair<AlleleIndexSetRef, AlleleIndexSetRef>;
    const static auto allele_index_pair_hasher = [] (const auto& p) {
        std::size_t result {};
        boost::hash_combine(result, p.first.get());
        boost::hash_combine(result, p.second.get());
        return result;
    };
    const static auto allele_index_pair_equal = [] (const auto& lhs, const auto& rhs) {
        return lhs.first.get() == rhs.first.get() && lhs.second.get() == rhs.second.get();
    };
    using AlleleIndexSetPairPosteriorMap = std::unordered_map<AlleleIndexSetRefPair, GenotypeChunkPosteriorMap,
                                                              decltype(allele_index_pair_hasher), decltype(allele_index_pair_equal)>;
    AlleleIndexSetPairPosteriorMap chunk_set_posteriors {10, allele_index_pair_hasher, allele_index_pair_equal};
    for (std::size_t g {0}; g < genotypes.size(); ++g) {
        if (is_heterozygous(g, info[lhs]) && is_heterozygous(g, info[rhs])) {
            collapse(chunks[g]);
            const auto index_pair = std::make_pair(std::cref(info[lhs][g].alleles), std::cref(info[rhs][g].alleles));
            chunk_set_posteriors[index_pair][chunks[g]] += genotype_posteriors[genotypes[g]];
        }
    }
    std::vector<std::vector<double>> result {};
    result.reserve(chunk_set_posteriors.size());
    for (const auto& p : chunk_set_posteriors) {
        result.push_back(extract_values(p.second));
    }
    return result;
}

auto compute_phase_quality(const std::vector<CompressedGenotype>& genotypes,
                           const std::vector<GenomicRegion>& sites,
                           const std::size_t lhs, const std::size_t rhs,
                           const Phaser::SampleGenotypePosteriorMap& genotype_posteriors,
                           const GenotypeInfoMatrix& info)
{
    if (overlaps(sites[lhs], sites[rhs])
     || is_very_likely_homozygous(sites[lhs], genotype_posteriors, info[lhs])
     || is_very_likely_homozygous(sites[rhs], genotype_posteriors, info[rhs])) {
        return probability_false_to_phred(0.0); // maximum quality
    }
    auto chunk_set_posteriors = compute_chunk_set_posteriors(genotypes, sites, lhs, rhs, genotype_posteriors, info);
    std::vector<double> set_weights(chunk_set_posteriors.size());
    const static auto sum_probabilities = [] (const auto& probs) {
         return std::accumulate(std::cbegin(probs), std::cend(probs), 0.0); };
    std::transform(std::cbegin(chunk_set_posteriors), std::cend(chunk_set_posteriors), 
                   std::begin(set_weights), sum_probabilities);
    double total_not_map_posterior {0};
    // subnormal numbers can cause divide by zero problems here when ffast-math is used.
    const auto heterozygous_mass = maths::normalise(set_weights);
    if (!maths::is_subnormal(heterozygous_mass) && heterozygous_mass > 0) {
        for (std::size_t set_idx {0}; set_idx < chunk_set_posteriors.size(); ++set_idx) {
            auto& posteriors = chunk_set_posteriors[set_idx];
            if (posteriors.size() > 1 && !maths::is_subnormal(maths::normalise(posteriors))) {
                for (auto& p : posteriors) p *= set_weights[set_idx];
                assert(!posteriors.empty());
                const auto map_posterior_itr = std::max_element(std::cbegin(posteriors), std::cend(posteriors));
                const auto not_map_posterior = std::accumulate(std::cbegin(posteriors), map_posterior_itr,    
                                               std::accumulate(std::next(map_posterior_itr), std::cend(posteriors), 0.0));
                total_not_map_posterior += not_map_posterior;
            }
        }
        total_not_map_posterior *= heterozygous_mass;
    }
    return probability_false_to_phred(total_not_map_posterior);
}

} // namespace

template <typename Graph>
struct MaxCliqueVisitor
{
    template <typename Clique>
    void clique(const Clique& clique, const Graph& g)
    {
        std::cout << "clique: "; for (auto v : clique) std::cout << v << ' '; std::cout << std::endl;
        if (clique.size() > max_clique.size()) {
            max_clique.assign(std::cbegin(clique), std::cend(clique));
        }
    }
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    std::vector<Vertex>& max_clique;
};

template <typename Graph>
struct CliqueRecorder
{
    template <typename Clique_>
    void clique(const Clique_& clique, const Graph& g)
    {
        cliques.emplace_back(std::cbegin(clique), std::cend(clique));
    }
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using Clique = std::deque<Vertex>;
    using CliqueList = std::deque<Clique>;
    CliqueList& cliques;
};

using PhaseQualityTable = std::vector<std::vector<Phred<double>>>;

Phred<double>
calculate_phase_quality(const std::vector<std::size_t>& phase_set,
                        const PhaseQualityTable& pairwise_phase_qualities)
{
    auto result = probability_false_to_phred(0.0);
    if (phase_set.size() > 1) {
        for (auto lhs_itr = std::cbegin(phase_set), penultimate_itr = std::prev(std::cend(phase_set)); lhs_itr != penultimate_itr; ++lhs_itr) {
            for (auto rhs_itr = std::next(lhs_itr); rhs_itr != std::cend(phase_set); ++rhs_itr) {
                result = std::min(result, pairwise_phase_qualities[*lhs_itr][*rhs_itr]);
            }
        }
    }
    return result;
}

namespace debug {

void write_phase_graph(const std::vector<GenomicRegion>& sites,
                       const PhaseQualityTable& pairwise_phase_qualities,
                       const Phred<double> min_phase_quality,
                       std::ofstream& graph_out)
{
    using PhaseGraph = boost::adjacency_matrix<boost::undirectedS, std::size_t, Phred<double>>;
    PhaseGraph phase_graph(pairwise_phase_qualities.size());
    for (std::size_t lhs {0}; lhs < pairwise_phase_qualities.size() - 1; ++lhs) {
        for (auto rhs = lhs + 1; rhs < pairwise_phase_qualities.size(); ++rhs) {
            boost::add_edge(lhs, rhs, pairwise_phase_qualities[lhs][rhs], phase_graph);
        }
    }
    const auto vertex_writer = [&] (std::ostream& out, auto v) {
        out << " [label=\"" << sites[v] << "\"]" << std::endl;
    };
    const auto edge_writer = [&] (std::ostream& out, auto e) {
        if (phase_graph[e] >= min_phase_quality) {
            out << " [color=green]" << std::endl;
        } else {
            out << " [style=dotted,color=red]" << std::endl;
        }
        out << " [label=\"" << phase_graph[e] << "\"]" << std::endl;
    };
    boost::write_graphviz(graph_out, phase_graph, vertex_writer, edge_writer);
}

} // namespace

Phaser::PhaseSetVector
Phaser::phase_sample(const std::vector<GenomicRegion>& sites,
                     const std::vector<CompressedGenotype>& genotypes,
                     const SampleGenotypePosteriorMap& genotype_posteriors) const
{
    using std::cbegin; using std::cend;
    using CompletePhaseGraph = boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, std::size_t>;
    using CompletePhaseGraphVertex = boost::graph_traits<CompletePhaseGraph>::vertex_descriptor;
    const auto genotype_info = compute_genotype_info(genotypes, sites);
    CompletePhaseGraph phase_graph {};
    std::vector<CompletePhaseGraphVertex> vertices(sites.size());
    for (std::size_t idx {0}; idx < sites.size(); ++idx) {
        vertices[idx] = boost::add_vertex(idx, phase_graph);
    }
    PhaseQualityTable pairwise_phase_qualities(sites.size(), PhaseQualityTable::value_type(sites.size()));
    for (std::size_t lhs_region_idx {0}; lhs_region_idx < sites.size() - 1; ++lhs_region_idx) {
        for (auto rhs_region_idx = lhs_region_idx + 1; rhs_region_idx < sites.size(); ++rhs_region_idx) {
            const auto phase_quality = compute_phase_quality(genotypes, sites, lhs_region_idx, rhs_region_idx,
                                                             genotype_posteriors, genotype_info);
            if (phase_quality >= config_.min_phase_quality) {
                boost::add_edge(vertices[lhs_region_idx], vertices[rhs_region_idx], phase_graph);
            }
            pairwise_phase_qualities[lhs_region_idx][rhs_region_idx] = phase_quality;
            pairwise_phase_qualities[rhs_region_idx][lhs_region_idx] = phase_quality;
        }
    }
//    std::string phase_graph_dot_filename {"/Users/dcooke/Genomics/octopus/scratch/phase_graph"};
//    phase_graph_dot_filename += "_" + to_string(contig_region(encompassing_region(sites)));
//    phase_graph_dot_filename += ".dot";
//    std::ofstream phase_graph_dot {phase_graph_dot_filename};
//    debug::write_phase_graph(sites, pairwise_phase_qualities, config_.min_phase_quality, phase_graph_dot);
    std::vector<std::size_t> fully_connected_vertices {}, not_fully_connected_vertices {};
    fully_connected_vertices.reserve(sites.size());
    not_fully_connected_vertices.reserve(sites.size());
    for (std::size_t vertex_idx {0}; vertex_idx < sites.size(); ++vertex_idx) {
        if (boost::degree(vertices[vertex_idx], phase_graph) == boost::num_vertices(phase_graph) - 1) {
            fully_connected_vertices.push_back(vertex_idx);
            boost::clear_vertex(vertices[vertex_idx], phase_graph);
            boost::remove_vertex(vertices[vertex_idx], phase_graph);
        } else {
            not_fully_connected_vertices.push_back(vertex_idx);
        }
    }
    PhaseSetVector result {};
    if (not_fully_connected_vertices.empty()) {
        // everything phased
        auto phase_quality = calculate_phase_quality(fully_connected_vertices, pairwise_phase_qualities);
        result.push_back({std::move(fully_connected_vertices), phase_quality});
    } else {
        std::vector<std::size_t> singleton_vertices {}, partially_connected_vertices {};
        singleton_vertices.reserve(not_fully_connected_vertices.size());
        partially_connected_vertices.reserve(not_fully_connected_vertices.size());
        for (const auto vertex_idx : not_fully_connected_vertices) {
            if (boost::degree(vertices[vertex_idx], phase_graph) == 0) {
                singleton_vertices.push_back(vertex_idx);
                boost::remove_vertex(vertices[vertex_idx], phase_graph);
            } else {
                partially_connected_vertices.push_back(vertex_idx);
            }
        }
        not_fully_connected_vertices.clear();
        not_fully_connected_vertices.shrink_to_fit();
        fully_connected_vertices.shrink_to_fit();
        partially_connected_vertices.shrink_to_fit();
        singleton_vertices.shrink_to_fit();
        using PartialPhaseGraph = boost::adjacency_matrix<boost::undirectedS>;
        CliqueRecorder<PartialPhaseGraph>::CliqueList cliques {};
        if (!partially_connected_vertices.empty()) {
            PartialPhaseGraph partial_phase_graph(partially_connected_vertices.size());
            for (std::size_t lhs_idx {0}; lhs_idx < partially_connected_vertices.size() - 1; ++lhs_idx) {
                for (std::size_t rhs_idx {0}; rhs_idx < partially_connected_vertices.size(); ++rhs_idx) {
                    const auto lhs_vertex = vertices[partially_connected_vertices[lhs_idx]];
                    const auto rhs_vertex = vertices[partially_connected_vertices[rhs_idx]];
                    if (boost::edge(lhs_vertex, rhs_vertex, phase_graph).second) {
                        boost::add_edge(lhs_idx, rhs_idx, partial_phase_graph);
                    }
                }
            }
            CliqueRecorder<PartialPhaseGraph> clique_recorder {cliques};
            boost::bron_kerbosch_all_cliques(partial_phase_graph, clique_recorder);
            for (auto& clique : cliques) {
                for (auto& idx : clique) {
                    idx = partially_connected_vertices[idx];
                }
                std::sort(std::begin(clique), std::end(clique));
            }
        }
        std::vector<std::vector<std::size_t>> region_possible_cliques(sites.size());
        for (const auto region_idx : partially_connected_vertices) {
            for (std::size_t clique_idx {0}; clique_idx < cliques.size(); ++clique_idx) {
                const auto& clique = cliques[clique_idx];
                if (std::binary_search(std::cbegin(clique), std::cend(clique), region_idx)) {
                    region_possible_cliques[region_idx].push_back(clique_idx);
                }
            }
        }
        for (std::size_t i {0}; i < singleton_vertices.size(); ++i) {
            region_possible_cliques[singleton_vertices[i]].push_back(cliques.size() + i);
        }
        for (const auto idx : singleton_vertices) {
            cliques.push_back({idx});
        }
        for (const auto region_idx : fully_connected_vertices) {
            region_possible_cliques[region_idx].resize(cliques.size());
            std::iota(std::begin(region_possible_cliques[region_idx]), std::end(region_possible_cliques[region_idx]), 0);
        }
        std::vector<std::vector<std::size_t>> phase_sets {};
        phase_sets.resize(cliques.size());
        for (std::size_t region_idx {0}; region_idx < sites.size(); ++region_idx) {
            const auto& possible_cliques = region_possible_cliques[region_idx];
            if (possible_cliques.size() == 1) {
                phase_sets[possible_cliques[0]].push_back(region_idx);
            } else {
                auto selected_phase_set_idx = possible_cliques.front();
                GenomicRegion::Distance min_region_distance {-1};
                for (const auto clique_idx : possible_cliques) {
                    for (const auto other_region_idx : cliques[clique_idx]) {
                        if (region_possible_cliques[other_region_idx].size() == 1) {
                            const auto region_distance = std::abs(inner_distance(sites[region_idx], sites[other_region_idx]));
                            if (min_region_distance < 0 || region_distance < min_region_distance) {
                                selected_phase_set_idx = clique_idx;
                                min_region_distance = region_distance;
                            }
                        }
                    }
                }
                phase_sets[selected_phase_set_idx].push_back(region_idx);
            }
        }
        phase_sets.erase(std::remove_if(std::begin(phase_sets), std::end(phase_sets),
                                        [] (const auto& c) { return c.empty(); }), std::end(phase_sets));
        const static auto front_less = [] (const auto& lhs, const auto& rhs) { return lhs.front() < rhs.front(); };
        std::sort(std::begin(phase_sets), std::end(phase_sets), front_less);
        result.reserve(phase_sets.size());
        for (auto& phase_set : phase_sets) {
            auto phase_quality = calculate_phase_quality(phase_set, pairwise_phase_qualities);
            result.push_back({std::move(phase_set), phase_quality});
        }
    }
    return result;
}

// non-member methods
    
namespace debug {

void print_phase_sets(const Phaser::PhaseSetMap& phasings, const std::vector<GenomicRegion>& variation_sites)
{
    print_phase_sets(std::cout, phasings, variation_sites);
}

} // namespace debug

} // namespace Ocotpus
