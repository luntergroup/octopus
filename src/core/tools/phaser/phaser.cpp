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
    #pragma gcc diagnostic push
    #pragma gcc diagnostic ignored "-Wunused-variable"
        #include <boost/graph/bron_kerbosch_all_cliques.hpp>
    #pragma gcc diagnostic pop
#endif // defined (__clang__)

#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"

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

auto min_phase_score(double p)
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

auto min_phase_score(const CompressedGenotype& called_genotype, const Phaser::SampleGenotypePosteriorMap& genotype_posteriors)
{
    return min_phase_score(genotype_posteriors[called_genotype]);
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
    for (auto& p : genotypes) p.second = collapse(p.second);
}

} // namespace

Phaser::PhaseSet
Phaser::phase(const MappableBlock<Haplotype>& haplotypes,
              const GenotypePosteriorMap& genotype_posteriors,
              const std::vector<GenomicRegion>& variation_regions,
              boost::optional<GenotypeCallMap> genotype_calls) const
{
    assert(!haplotypes.empty());
    assert(!genotype_posteriors.empty1() && !genotype_posteriors.empty2());
    assert(std::is_sorted(std::cbegin(variation_regions), std::cend(variation_regions)));
    const auto& haplotype_region = mapped_region(haplotypes);
    const auto partitions = extract_covered_regions(variation_regions);
    auto genotypes = extract_keys(genotype_posteriors);
    unsigned min_genotype_ploidy, max_genotype_ploidy; std::tie(min_genotype_ploidy, max_genotype_ploidy) = minmax_ploidy(genotypes);
    PhaseSet result {haplotype_region};
    result.phase_regions.reserve(genotype_posteriors.size1());
    if (max_genotype_ploidy == 1 || partitions.size() == 1) {
        for (const auto& p : genotype_posteriors) {
            if (config_.max_phase_score) {
                result.phase_regions[p.first].emplace_back(haplotype_region, *config_.max_phase_score);
            } else {
                static const Phred<double> max_possible_score {Phred<double>::Probability {0.0}};
                result.phase_regions[p.first].emplace_back(haplotype_region, max_possible_score);
            }
        }
    } else {
        boost::optional<GenotypePosteriorMap> collapsed_genotype_posteriors {};
        for (const auto& p : genotype_posteriors) {
            const SampleName& sample {p.first};
            if (!collapsed_genotype_posteriors && genotype_calls && config_.max_phase_score
             && min_phase_score(genotype_calls->at(sample), p.second) >= *config_.max_phase_score) {
                result.phase_regions[sample].emplace_back(haplotype_region, *config_.max_phase_score);
            } else {
                if (!collapsed_genotype_posteriors && (max_genotype_ploidy > 2 || min_genotype_ploidy != max_genotype_ploidy)) {
                    collapsed_genotype_posteriors = marginalise_collapsed_genotypes(genotype_posteriors);
                    genotypes = extract_keys(*collapsed_genotype_posteriors);
                    std::tie(min_genotype_ploidy, max_genotype_ploidy) = minmax_ploidy(genotypes);
                    if (genotype_calls) collapse_each(*genotype_calls);
                }
                PhaseSet::SamplePhaseRegions phases;
                if (collapsed_genotype_posteriors) {
                    const auto& collapsed_sample_genotype_posteriors = (*collapsed_genotype_posteriors)[sample];
                    if (genotype_calls && config_.max_phase_score
                     && min_phase_score(genotype_calls->at(sample), collapsed_sample_genotype_posteriors) >= *config_.max_phase_score) {
                        result.phase_regions[sample].emplace_back(haplotype_region, *config_.max_phase_score);
                    } else {
                        phases = phase_sample(haplotype_region, partitions, genotypes, collapsed_sample_genotype_posteriors);
                    }
                } else {
                    phases = phase_sample(haplotype_region, partitions, genotypes, p.second);
                }
                if (config_.max_phase_score) {
                    for (auto& phase : phases) phase.score = std::min(phase.score, *config_.max_phase_score);
                }
                result.phase_regions.emplace(sample, std::move(phases));
            }
        }
    }
    return result;
}

namespace {

using GenotypeChunk        = SharedGenotype;
using PhaseComplementSet   = std::deque<GenotypeChunk>;
using PhaseComplementSets  = std::vector<PhaseComplementSet>;
using PartitionIterator    = std::vector<GenomicRegion>::const_iterator;
using GenotypeChunkVector  = std::vector<GenotypeChunk>;

struct GenotypeChunkPartitions
{
    GenotypeChunk genotype;
    GenotypeChunkVector chunks;
};

using GenotypeChunkTable = std::vector<GenotypeChunkPartitions>;

template <typename Container>
GenotypeChunkTable
make_chunk_table(const Container& genotypes, PartitionIterator first_parition, PartitionIterator last_partition)
{
    const auto num_partitions = std::distance(first_parition, last_partition);
    GenotypeChunkTable result {};
    result.reserve(genotypes.size());
    for (const auto& genotype : genotypes) {
        result.push_back({genotype, GenotypeChunkVector {}});
        result.back().chunks.reserve(num_partitions);
    }
    std::for_each(first_parition, last_partition, [&] (const auto& region) {
        auto chunks = copy_each<GenotypeChunk::ElementType>(genotypes, region);
        for (std::size_t i {0}; i < genotypes.size(); ++i) {
            result[i].chunks.push_back(std::move(chunks[i]));
        }
    });
    return result;
}

auto max_ploidy(const GenotypeChunkTable& table) noexcept
{
    assert(!table.empty());
    const static auto ploidy_less = [] (const auto& lhs, const auto& rhs) { return lhs.genotype.ploidy() < rhs.genotype.ploidy(); };
    return std::max_element(std::cbegin(table), std::cend(table), ploidy_less)->genotype.ploidy();
}

struct GenotypeChunkVectorExactHash
{
    auto operator()(const GenotypeChunkVector& chunks) const noexcept
    {
        std::size_t result {};
        for (const auto& genotype : chunks) {
            using boost::hash_combine;
            hash_combine(result, std::hash<GenotypeChunk>{}(genotype));
        }
        return result;
    }
};
struct GenotypeUniqueHash
{
    std::size_t operator()(const GenotypeChunk& genotype) const
    {
        std::size_t result {0};
        for (const auto& haplotype : collapse(genotype)) {
            using boost::hash_combine;
            hash_combine(result, std::hash<GenotypeChunk::ElementType>{}(haplotype));
        }
        return result;
    }
};
struct GenotypeChunkVectorUniqueHash
{
    auto operator()(const GenotypeChunkVector& chunks) const noexcept
    {
        std::size_t result {};
        for (const auto& genotype : chunks) {
            using boost::hash_combine;
            hash_combine(result, GenotypeUniqueHash{}(genotype));
        }
        return result;
    }
};
struct GenotypeChunkVectorExactEqual
{
    bool operator()(const GenotypeChunkVector& lhs, const GenotypeChunkVector& rhs) const
    {
        return lhs == rhs;
    }
};

template <typename MappableSet>
bool have_same_element_set(const Genotype<MappableSet>& lhs, const Genotype<MappableSet>& rhs) noexcept
{
    if (lhs.ploidy() < 3 && rhs.ploidy() < 3) {
        return lhs == rhs;
    } else {
        return collapse(lhs) == collapse(rhs);
    }
}

struct GenotypeChunkUniqueEqual
{
    bool operator()(const GenotypeChunk& lhs, const GenotypeChunk& rhs) const
    {
        return have_same_element_set(lhs, rhs);
    }
};
struct GenotypeChunkVectorUniqueEqual
{
    bool operator()(const GenotypeChunkVector& lhs, const GenotypeChunkVector& rhs) const
    {
        return lhs.size() == rhs.size() && std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), GenotypeChunkUniqueEqual {});
    }
};

template <typename Hash, typename Equal>
using PhaseSetMap = std::unordered_map<GenotypeChunkVector, PhaseComplementSet, Hash, Equal>;

template <typename Map>
PhaseComplementSets
make_phase_sets_helper(const GenotypeChunkTable& table, Map& phase_sets)
{
    phase_sets.reserve(table.size());
    for (const auto& row : table) {
        phase_sets[row.chunks].push_back(row.genotype);
    }
    PhaseComplementSets result {};
    result.reserve(phase_sets.size());
    for (auto&& p : phase_sets) {
        result.push_back(std::move(p.second));
    }
    return result;
}

PhaseComplementSets
make_phase_sets(const GenotypeChunkTable& table, const Phaser::GenotypeMatchType match_type)
{
    if (match_type == Phaser::GenotypeMatchType::exact || max_ploidy(table) < 3) {
        // if haploid or diploid genotypes then always use exact match as exact and unique match are
        // identical and exact match is faster
        PhaseSetMap<GenotypeChunkVectorExactHash, GenotypeChunkVectorExactEqual> phase_sets {};
        return make_phase_sets_helper(table, phase_sets);
    } else {
        PhaseSetMap<GenotypeChunkVectorUniqueHash, GenotypeChunkVectorUniqueEqual> phase_sets {};
        return make_phase_sets_helper(table, phase_sets);
    }
}

template <typename Container>
PhaseComplementSets
generate_phase_complement_sets(const Container& genotypes,
                               PartitionIterator first_parition, PartitionIterator last_partition,
                               const Phaser::GenotypeMatchType match_type)
{
    const auto table = make_chunk_table(genotypes, first_parition, last_partition);
    return make_phase_sets(table, match_type);
}

template <typename Map>
double marginalise(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    return std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                           [&] (const auto curr, const auto& genotype) {
                               return curr + genotype_posteriors.at(genotype);
                           });
}

template <typename Map>
double calculate_entropy(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    const auto norm = marginalise(phase_set, genotype_posteriors);
    if (norm <= 0.0) {
        // if norm ~= 0 then every element in the phase must must have probability ~= 0, so it
        // just looks like a uniform distirbution
        return maximum_entropy(phase_set.size());
    }
    return std::max(0.0, -std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                                          [&genotype_posteriors, norm] (const auto curr, const auto& genotype) {
                                              const auto p = genotype_posteriors.at(genotype) / norm;
                                              return curr + p * std::log2(std::max(p, std::numeric_limits<double>::min()));
                                          }));
}

template <typename Map>
double calculate_relative_entropy(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    if (phase_set.size() < 2) return 1.0;
    return 1.0 - std::min(calculate_entropy(phase_set, genotype_posteriors) / maximum_entropy(2), 1.0);
}

template <typename Map>
auto calculate_phase_score(const PhaseComplementSet& phase_set, const Map& genotype_posteriors)
{
    return marginalise(phase_set, genotype_posteriors) * calculate_relative_entropy(phase_set, genotype_posteriors);
}

template <typename Map>
Phred<double> calculate_phase_score(const PhaseComplementSets& phase_sets, const Map& genotype_posteriors)
{
    return Phred<double> { Phred<double>::Probability {
    std::max(0.0, 1.0 - std::accumulate(std::cbegin(phase_sets), std::cend(phase_sets), 0.0,
                                        [&] (const auto curr, const auto& phase_set) {
                                            return curr + calculate_phase_score(phase_set, genotype_posteriors);
                                        }))
    }};
}

using GenotypeChunkPosteriorMap = std::unordered_map<GenotypeChunk, double>;

auto copy_and_marginalise(const std::vector<CompressedGenotype>& genotypes,
                          const Phaser::SampleGenotypePosteriorMap& genotype_posteriors,
                          const GenomicRegion& region)
{
    auto chunks = copy_each<GenotypeChunk::ElementType>(genotypes, region);
    GenotypeChunkPosteriorMap chunk_posteriors {genotypes.size()};
    for (std::size_t g {0}; g < genotypes.size(); ++g) {
        chunk_posteriors[chunks[g]] += genotype_posteriors[genotypes[g]];
    }
    std::sort(std::begin(chunks), std::end(chunks), GenotypeLess {});
    chunks.erase(std::unique(std::begin(chunks), std::end(chunks)), std::end(chunks));
    return std::make_pair(std::move(chunks), std::move(chunk_posteriors));
}

template <typename Range>
auto marginalise(const Range& genotypes,
                 const std::vector<GenomicRegion>& regions,
                 const Phaser::SampleGenotypePosteriorMap& genotype_posteriors)
{
    auto chunks = copy_each(genotypes, regions);
    for (auto& genotype : chunks) {
        genotype = genotype.collapse();
    }
    GenotypeChunkPosteriorMap chunk_posteriors {genotypes.size()};
    for (std::size_t g {0}; g < genotypes.size(); ++g) {
        chunk_posteriors[chunks[g]] += genotype_posteriors[genotypes[g]];
    }
    std::sort(std::begin(chunks), std::end(chunks), GenotypeLess {});
    chunks.erase(std::unique(std::begin(chunks), std::end(chunks)), std::end(chunks));
    std::vector<double> posteriors(chunks.size());
    std::transform(std::cbegin(chunks), std::cend(chunks), std::begin(posteriors),
                   [&] (const auto& chunk) { return chunk_posteriors[chunk]; });
    maths::normalise(posteriors);
    return posteriors;
}

template <typename Range>
auto compute_phase_quality(const Range& genotypes,
                           const std::vector<GenomicRegion>& regions,
                           const Phaser::SampleGenotypePosteriorMap& genotype_posteriors)
{
    const auto marginal_posteriors = marginalise(genotypes, regions, genotype_posteriors);
    assert(!marginal_posteriors.empty());
    const auto map_marginal_posterior_itr = std::max_element(std::cbegin(marginal_posteriors), std::cend(marginal_posteriors));
    const auto not_map_posterior = std::accumulate(std::cbegin(marginal_posteriors), map_marginal_posterior_itr,
                                   std::accumulate(std::next(map_marginal_posterior_itr), std::cend(marginal_posteriors), 0.0));
    assert(not_map_posterior > 0);
    return probability_false_to_phred(not_map_posterior);
}

template <typename Map>
void print(const PhaseComplementSets& phase_sets, const Map& genotype_posteriors)
{
    for (std::size_t i {0}; i < phase_sets.size(); ++i) {
        const auto& phase_set = phase_sets[i];
        std::cout << "Phase set " << i << ":" << std::endl;;
        for (const auto& genotype : phase_set) {
            octopus::debug::print_variant_alleles(genotype); std::cout << ' ' << genotype_posteriors.at(genotype) << std::endl;
        }
        std::cout << "entropy = " << calculate_entropy(phase_set, genotype_posteriors) << std::endl;
        std::cout << "probability mass = " << marginalise(phase_set, genotype_posteriors) << std::endl;
    }
}

template <typename T, typename Container, typename BinaryPredicate = std::less<>>
auto
insert_sorted(T value, Container& values, BinaryPredicate compare = std::less<> {})
{
    auto position = std::upper_bound(std::begin(values), std::end(values), value, compare);
    return values.insert(position, std::move(value));
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

template <typename UnaryPredicate, typename Graph>
void remove_vertex_if(UnaryPredicate&& pred, Graph& g)
{
    typename boost::graph_traits<Graph>::vertex_iterator vi, vi_end, vi_next;
    std::tie(vi, vi_end) = boost::vertices(g);
    for (vi_next = vi; vi != vi_end; vi = vi_next) {
        ++vi_next;
        if (pred(*vi)) {
            boost::remove_vertex(*vi, g);
        }
    }
}

Phaser::PhaseSet::SamplePhaseRegions
Phaser::phase_sample(const GenomicRegion& region,
                     const std::vector<GenomicRegion>& partitions,
                     const std::vector<CompressedGenotype>& genotypes,
                     const SampleGenotypePosteriorMap& genotype_posteriors) const
{
    using std::cbegin; using std::cend;
    std::vector<GenomicRegion> partition_pair {};
    using CompletePhaseGraph = boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, std::size_t, Phred<double>>;
    using CompletePhaseGraphVertex = boost::graph_traits<CompletePhaseGraph>::vertex_descriptor;
    CompletePhaseGraph phase_graph {};
    std::vector<CompletePhaseGraphVertex> vertices(partitions.size());
    for (std::size_t idx {0}; idx < partitions.size(); ++idx) {
        vertices[idx] = boost::add_vertex(idx, phase_graph);
    }
    for (std::size_t lhs_partition_idx {0}; lhs_partition_idx < partitions.size() - 1; ++lhs_partition_idx) {
        for (auto rhs_partition_idx = lhs_partition_idx + 1; rhs_partition_idx < partitions.size(); ++rhs_partition_idx) {
            partition_pair.assign({partitions[lhs_partition_idx], partitions[rhs_partition_idx]});
            const auto phase_quality = compute_phase_quality(genotypes, partition_pair, genotype_posteriors);
            if (phase_quality >= config_.min_phase_score) {
                boost::add_edge(vertices[lhs_partition_idx], vertices[rhs_partition_idx], phase_quality, phase_graph);
            }
        }
    }
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
    std::vector<std::vector<std::size_t>> phase_sets {};
    if (not_fully_connected_vertices.empty()) {
        phase_sets.push_back(std::move(fully_connected_vertices)); // everything phased
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
    }
    
    std::cout << "phase sets" << std::endl;
    for (const auto& ps : phase_sets) {
        for (auto partition_idx : ps) {
            std::cout << partition_idx << " (" << partitions[partition_idx] << ") ";
        }
        std::cout << std::endl;
    }
    
//    std::ofstream phase_graph_dot {"/Users/dcooke/Genomics/octopus/scratch/phase_graph.dot"};
//    const auto edge_writer = [&phase_graph] (std::ostream& out, auto e) {
//        if (phase_graph[e].score() >= 5) {
//            out << " [color=green]" << std::endl;
//        } else {
//            out << " [style=dotted,color=red]" << std::endl;
//        }
//        out << " [label=\"" << phase_graph[e] << "\"]" << std::endl;
//    };
//    boost::write_graphviz(phase_graph_dot, phase_graph, boost::default_writer(), edge_writer);
////    boost::write_graphviz(phase_graph_dot, phase_graph, boost::default_writer(), boost::make_label_writer(boost::get(boost::edge_bundle, phase_graph)));
    
    auto first_partition = std::cbegin(partitions);
    auto last_partition  = std::cend(partitions);
    Phaser::PhaseSet::SamplePhaseRegions result {};
    std::vector<GenotypeChunk> chunks;
    GenotypeChunkPosteriorMap chunk_posteriors;
    while (first_partition != std::cend(partitions)) {
        auto curr_region = encompassing_region(first_partition, last_partition);
        std::tie(chunks, chunk_posteriors) = copy_and_marginalise(genotypes, genotype_posteriors, curr_region);
        auto phase_set = generate_phase_complement_sets(chunks, first_partition, last_partition, config_.genotype_match);
        auto phase_score = calculate_phase_score(phase_set, chunk_posteriors);
        if (phase_score >= config_.min_phase_score || std::distance(first_partition, last_partition) == 1) {
            result.emplace_back(encompassing_region(first_partition, last_partition), phase_score);
            first_partition = last_partition;
            last_partition  = std::cend(partitions);
        } else {
            --last_partition;
        }
    }
    return result;
}

// non-member methods

bool is_split_phasing(const Phaser::PhaseSet& phase)
{
    return std::any_of(std::cbegin(phase.phase_regions), std::cend(phase.phase_regions),
                       [] (const auto& p) { return p.second.size() > 1; });
}
    
namespace debug {

void print_phase_sets(const Phaser::PhaseSet& phasings)
{
    print_phase_sets(std::cout, phasings);
}

} // namespace debug

} // namespace Ocotpus
