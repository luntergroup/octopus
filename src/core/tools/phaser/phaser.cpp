// Copyright (c) 2015-2019 Daniel Cooke
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

} // namespace

Phaser::PhaseSetMap
Phaser::phase(const MappableBlock<Haplotype>& haplotypes,
              const GenotypePosteriorMap& genotype_posteriors,
              const std::vector<GenomicRegion>& variation_regions,
              boost::optional<GenotypeCallMap> genotype_calls) const
{
    assert(!haplotypes.empty());
    assert(!genotype_posteriors.empty1() && !genotype_posteriors.empty2());
    assert(std::is_sorted(std::cbegin(variation_regions), std::cend(variation_regions)));
    const auto partitions = extract_covered_regions(variation_regions);
    assert(!partitions.empty() && partitions.size() <= variation_regions.size());
    auto genotypes = extract_keys(genotype_posteriors);
    unsigned min_genotype_ploidy, max_genotype_ploidy; std::tie(min_genotype_ploidy, max_genotype_ploidy) = minmax_ploidy(genotypes);
    PhaseSetMap result {};
    const auto num_samples = genotype_posteriors.size1();
    result.reserve(num_samples);
    std::vector<std::size_t> max_phase_set(variation_regions.size());
    std::iota(std::begin(max_phase_set), std::end(max_phase_set), std::size_t {0});
    static const Phred<double> max_possible_quality {Phred<double>::Probability {0.0}};
    if (max_genotype_ploidy == 1 || partitions.size() == 1) {
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
                        phase_sets = phase_sample(partitions, genotypes, collapsed_sample_genotype_posteriors);
                    }
                } else {
                    phase_sets = phase_sample(partitions, genotypes, p.second);
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
    if (partitions.size() < variation_regions.size()) {
        std::vector<std::vector<std::size_t>> parition_index_to_variation_site_index_table {};
        for (auto& p : result) {
            assert(!p.second.empty());
            if (p.second.size() > 1 || p.second.front().site_indices.size() < variation_regions.size()) {
                // Need to convert partition indices back to site indices
                parition_index_to_variation_site_index_table.resize(partitions.size());
                for (auto& phase_set : p.second) {
                    auto& partition_indices = phase_set.site_indices;
                    for (auto partition_index_itr = std::begin(partition_indices);
                         partition_index_itr != std::end(partition_indices);
                         ++partition_index_itr) {
                        auto& site_indices = parition_index_to_variation_site_index_table[*partition_index_itr];
                        if (site_indices.empty()) {
                            // lazy evaluate this
                            const auto overlapped_sites = bases(overlap_range(variation_regions, partitions[*partition_index_itr]));
                            auto first_site_index = static_cast<std::size_t>(std::distance(std::cbegin(variation_regions), std::cbegin(overlapped_sites)));
                            site_indices.resize(size(overlapped_sites));
                            std::iota(std::begin(site_indices), std::end(site_indices), first_site_index);
                        }
                        assert(!site_indices.empty());
                        assert(site_indices.front() >= *partition_index_itr);
                        if (site_indices.size() == 1) {
                            *partition_index_itr = site_indices.front();
                        } else {
                            *partition_index_itr = site_indices.front();
                            ++partition_index_itr;
                            partition_index_itr = partition_indices.insert(partition_index_itr,
                                                                      std::next(std::cbegin(site_indices)), std::cend(site_indices));
                            std::advance(partition_index_itr, site_indices.size() - 2);
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

bool is_very_likely_homozygous(const GenomicRegion& region, const Phaser::SampleGenotypePosteriorMap& genotype_posteriors)
{
    assert(!genotype_posteriors.empty());
    const static auto posterior_less = [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; };
    const auto map_posterior_itr = std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), posterior_less);
    const auto map_posterior = map_posterior_itr->second;
    return map_posterior > 0.9999 && is_homozygous(map_posterior_itr->first, region);
}

auto marginalise(const std::vector<CompressedGenotype>& genotypes,
                 const std::vector<GenomicRegion>& regions,
                 const Phaser::SampleGenotypePosteriorMap& genotype_posteriors)
{
    auto chunks = copy_each<GenotypeChunk::ElementType>(genotypes, regions);
    for (auto& genotype : chunks) {
        collapse(genotype);
    }
    GenotypeChunkPosteriorMap chunk_posteriors {genotypes.size()};
    for (std::size_t g {0}; g < genotypes.size(); ++g) {
        chunk_posteriors[chunks[g]] += genotype_posteriors[genotypes[g]];
    }
    auto result = extract_values(chunk_posteriors);
    maths::normalise(result);
    return result;
}

auto compute_phase_quality(const std::vector<CompressedGenotype>& genotypes,
                           const std::vector<GenomicRegion>& regions,
                           const Phaser::SampleGenotypePosteriorMap& genotype_posteriors)
{
    assert(regions.size() == 2);
    if (is_very_likely_homozygous(regions.front(), genotype_posteriors)
     || is_very_likely_homozygous(regions.back(), genotype_posteriors)) {
        return probability_false_to_phred(0.0);
    }
    const auto marginal_posteriors = marginalise(genotypes, regions, genotype_posteriors);
    assert(!marginal_posteriors.empty());
    const auto map_marginal_posterior_itr = std::max_element(std::cbegin(marginal_posteriors), std::cend(marginal_posteriors));
    const auto not_map_posterior = std::accumulate(std::cbegin(marginal_posteriors), map_marginal_posterior_itr,
                                   std::accumulate(std::next(map_marginal_posterior_itr), std::cend(marginal_posteriors), 0.0));
    return probability_false_to_phred(not_map_posterior);
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
    if (phase_set.size() < 2) {
        return probability_false_to_phred(0.0); // max quality
    } else {
        double pairwise_quality_sum {0};
        unsigned num_pairs {0};
        for (std::size_t lhs {0}; lhs < phase_set.size() - 1; ++lhs) {
            for (auto rhs = lhs + 1; rhs < phase_set.size(); ++rhs) {
                pairwise_quality_sum += pairwise_phase_qualities[lhs][rhs].score();
                ++num_pairs;
            }
        }
        return Phred<double> {pairwise_quality_sum / num_pairs};
    }
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
Phaser::phase_sample(const std::vector<GenomicRegion>& partitions,
                     const std::vector<CompressedGenotype>& genotypes,
                     const SampleGenotypePosteriorMap& genotype_posteriors) const
{
    using std::cbegin; using std::cend;
    std::vector<GenomicRegion> partition_pair {};
    using CompletePhaseGraph = boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, std::size_t>;
    using CompletePhaseGraphVertex = boost::graph_traits<CompletePhaseGraph>::vertex_descriptor;
    CompletePhaseGraph phase_graph {};
    std::vector<CompletePhaseGraphVertex> vertices(partitions.size());
    for (std::size_t idx {0}; idx < partitions.size(); ++idx) {
        vertices[idx] = boost::add_vertex(idx, phase_graph);
    }
    PhaseQualityTable pairwise_phase_qualities(partitions.size(), PhaseQualityTable::value_type(partitions.size()));
    for (std::size_t lhs_partition_idx {0}; lhs_partition_idx < partitions.size() - 1; ++lhs_partition_idx) {
        for (auto rhs_partition_idx = lhs_partition_idx + 1; rhs_partition_idx < partitions.size(); ++rhs_partition_idx) {
            partition_pair.assign({partitions[lhs_partition_idx], partitions[rhs_partition_idx]});
            const auto phase_quality = compute_phase_quality(genotypes, partition_pair, genotype_posteriors);
            if (phase_quality >= config_.min_phase_quality) {
                boost::add_edge(vertices[lhs_partition_idx], vertices[rhs_partition_idx], phase_graph);
            }
            pairwise_phase_qualities[lhs_partition_idx][rhs_partition_idx] = phase_quality;
            pairwise_phase_qualities[rhs_partition_idx][lhs_partition_idx] = phase_quality;
        }
    }
//    std::string phase_graph_dot_filename {"/Users/dcooke/Genomics/octopus/scratch/phase_graph"};
//    phase_graph_dot_filename += "_" + to_string(contig_region(encompassing_region(partitions)));
//    phase_graph_dot_filename += ".dot";
//    std::ofstream phase_graph_dot {phase_graph_dot_filename};
//    debug::write_phase_graph(partitions, pairwise_phase_qualities, config_.min_phase_quality, phase_graph_dot);
    std::vector<std::size_t> fully_connected_vertices {}, not_fully_connected_vertices {};
    fully_connected_vertices.reserve(partitions.size());
    not_fully_connected_vertices.reserve(partitions.size());
    for (std::size_t vertex_idx {0}; vertex_idx < partitions.size(); ++vertex_idx) {
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
        std::vector<std::vector<std::size_t>> partition_possible_cliques(partitions.size());
        for (const auto partition_idx : partially_connected_vertices) {
            for (std::size_t clique_idx {0}; clique_idx < cliques.size(); ++clique_idx) {
                const auto& clique = cliques[clique_idx];
                if (std::binary_search(std::cbegin(clique), std::cend(clique), partition_idx)) {
                    partition_possible_cliques[partition_idx].push_back(clique_idx);
                }
            }
        }
        for (std::size_t i {0}; i < singleton_vertices.size(); ++i) {
            partition_possible_cliques[singleton_vertices[i]].push_back(cliques.size() + i);
        }
        for (const auto idx : singleton_vertices) {
            cliques.push_back({idx});
        }
        for (const auto partition_idx : fully_connected_vertices) {
            partition_possible_cliques[partition_idx].resize(cliques.size());
            std::iota(std::begin(partition_possible_cliques[partition_idx]), std::end(partition_possible_cliques[partition_idx]), 0);
        }
        std::vector<std::vector<std::size_t>> phase_sets {};
        phase_sets.resize(cliques.size());
        for (std::size_t partition_idx {0}; partition_idx < partitions.size(); ++partition_idx) {
            const auto& possible_cliques = partition_possible_cliques[partition_idx];
            if (possible_cliques.size() == 1) {
                phase_sets[possible_cliques[0]].push_back(partition_idx);
            } else {
                auto selected_phase_set_idx = possible_cliques.front();
                GenomicRegion::Distance min_partition_distance {-1};
                for (const auto clique_idx : possible_cliques) {
                    for (const auto other_partition_idx : cliques[clique_idx]) {
                        if (partition_possible_cliques[other_partition_idx].size() == 1) {
                            const auto partition_distance = std::abs(inner_distance(partitions[partition_idx], partitions[other_partition_idx]));
                            if (min_partition_distance < 0 || partition_distance < min_partition_distance) {
                                selected_phase_set_idx = clique_idx;
                                min_partition_distance = partition_distance;
                            }
                        }
                    }
                }
                phase_sets[selected_phase_set_idx].push_back(partition_idx);
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

void print_phase_sets(const Phaser::PhaseSetMap& phasings, const std::vector<GenomicRegion>& variation_regions)
{
    print_phase_sets(std::cout, phasings, variation_regions);
}

} // namespace debug

} // namespace Ocotpus
