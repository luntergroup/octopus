// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "assembler.hpp"

#include <iterator>
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <limits>
#include <cassert>
#include <iostream>

#include <boost/functional/hash.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/dominator_tree.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/exception.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/copy.hpp>

#include "ksp/yen_ksp.hpp"

#include "utils/sequence_utils.hpp"
#include "utils/append.hpp"
#include "utils/maths.hpp"

#define _unused(x) ((void)(x))

namespace octopus { namespace coretools {

namespace {

template <typename I>
auto sequence_length(const I num_kmers, const unsigned kmer_size) noexcept
{
    return num_kmers > 0 ? num_kmers + kmer_size - 1 : 0;
}

template <typename T>
std::size_t count_kmers(const T& sequence, const unsigned kmer_size) noexcept
{
    return sequence.size() >= kmer_size ? sequence.size() - kmer_size + 1 : 0;
}

} // namespace

// public methods

Assembler::NonCanonicalReferenceSequence::NonCanonicalReferenceSequence(NucleotideSequence reference_sequence)
: std::invalid_argument {"bad reference sequence"}
, reference_sequence_ {std::move(reference_sequence)}
{}

const char* Assembler::NonCanonicalReferenceSequence::what() const noexcept
{
    return std::invalid_argument::what();
}

Assembler::Assembler(const Parameters params)
: params_ {params}
, reference_kmers_ {}
, reference_head_position_ {0}
, reference_vertices_ {}
{}

Assembler::Assembler(const Parameters params, const NucleotideSequence& reference)
: params_ {params}
, reference_kmers_ {}
, reference_head_position_ {0}
, reference_vertices_ {}
{
    insert_reference_into_empty_graph(reference);
}

Assembler::Assembler(const Assembler& other)
: params_ {other.params_}
, reference_kmers_ {other.reference_kmers_}
, reference_head_position_ {other.reference_head_position_}
{
    std::unordered_map<Vertex, std::size_t> index_map {};
    index_map.reserve(boost::num_vertices(other.graph_));
    const auto p = boost::vertices(other.graph_);
    std::size_t i {0};
    std::for_each(p.first, p.second, [&i, &index_map] (const Vertex& v) { index_map.emplace(v, i++); });
    std::unordered_map<Vertex, Vertex> vertex_copy_map {};
    vertex_copy_map.reserve(boost::num_vertices(other.graph_));
    boost::copy_graph(other.graph_, this->graph_,
                      boost::vertex_index_map(boost::make_assoc_property_map(index_map))
                      .orig_to_copy(boost::make_assoc_property_map(vertex_copy_map)));
    vertex_cache_ = other.vertex_cache_;
    for (auto& p : vertex_cache_) {
        p.second = vertex_copy_map.at(p.second);
    }
    reference_vertices_ = other.reference_vertices_;
    for (auto& v : reference_vertices_) {
        v = vertex_copy_map.at(v);
    }
    for (const auto& other_edge : other.reference_edges_) {
        Edge e; bool e_in_graph;
        const auto u = vertex_copy_map.at(boost::source(other_edge, other.graph_));
        const auto v = vertex_copy_map.at(boost::target(other_edge, other.graph_));
        std::tie(e, e_in_graph) = boost::edge(u, v, graph_);
        assert(e_in_graph);
        reference_edges_.push_back(e);
    }
}

unsigned Assembler::kmer_size() const noexcept
{
    return params_.kmer_size;
}

Assembler::Parameters Assembler::params() const
{
    return params_;
}

void Assembler::insert_reference(const NucleotideSequence& sequence)
{
    if (sequence.size() >= kmer_size()) {
        if (is_empty()) {
            insert_reference_into_empty_graph(sequence);
        } else if (reference_kmers_.empty()) {
            insert_reference_into_populated_graph(sequence);
        } else {
            throw std::runtime_error {"Assembler: only one reference sequence can be inserted into the graph"};
        }
    } else {
        throw std::runtime_error {"Assembler:: reference length must >= kmer_size"};
    }
}

void Assembler::insert_read(const NucleotideSequence& sequence,
                            const BaseQualityVector& base_qualities,
                            const Direction strand,
                            const SampleID sample)
{
    if (sequence.size() >= kmer_size()) {
        const bool is_forward_strand {strand == Direction::forward};
        auto kmer_begin = std::cbegin(sequence);
        auto kmer_end   = std::next(kmer_begin, kmer_size());
        auto base_quality_itr = std::next(std::cbegin(base_qualities), kmer_size());
        Kmer prev_kmer {kmer_begin, kmer_end};
        bool prev_kmer_good {true};
        auto vertex_itr = vertex_cache_.find(prev_kmer);
        auto ref_kmer_itr = std::cbegin(reference_kmers_);
        if (vertex_itr == std::cend(vertex_cache_)) {
            const auto u = add_vertex(prev_kmer);
            if (!u) prev_kmer_good = false;
        } else if (is_reference(vertex_itr->second)) {
            ref_kmer_itr = std::find(std::cbegin(reference_kmers_), std::cend(reference_kmers_), prev_kmer);
            assert(ref_kmer_itr != std::cend(reference_kmers_));
            auto next_kmer_begin = std::next(kmer_begin);
            auto next_kmer_end   = std::next(kmer_end);
            const auto ref_offset = std::distance(std::cbegin(reference_kmers_), ref_kmer_itr);
            auto ref_vertex_itr = std::next(std::cbegin(reference_vertices_), ref_offset);
            auto ref_edge_itr = std::next(std::cbegin(reference_edges_), ref_offset);
            ++ref_kmer_itr;
            for (; next_kmer_end <= std::cend(sequence) && ref_kmer_itr < std::cend(reference_kmers_);
                   ++next_kmer_begin, ++next_kmer_end, ++ref_kmer_itr, ++ref_vertex_itr, ++ref_edge_itr, ++base_quality_itr) {
                if (std::equal(next_kmer_begin, next_kmer_end, std::cbegin(*ref_kmer_itr))) {
                    assert(ref_edge_itr != std::cend(reference_edges_));
                    increment_weight(*ref_edge_itr, is_forward_strand, *base_quality_itr, sample);
                } else {
                    break;
                }
            }
            if (next_kmer_end > std::cend(sequence)) {
                return;
            }
            kmer_begin = std::prev(next_kmer_begin);
            kmer_end   = std::prev(next_kmer_end);
            assert(kmer_end <= std::cend(sequence));
            prev_kmer = Kmer {kmer_begin, kmer_end};
        }
        ++kmer_begin;
        ++kmer_end;
        for (; kmer_end <= std::cend(sequence); ++kmer_begin, ++kmer_end, ++base_quality_itr) {
            Kmer kmer {kmer_begin, kmer_end};
            const auto kmer_itr = vertex_cache_.find(kmer);
            if (kmer_itr == std::cend(vertex_cache_)) {
                const auto v = add_vertex(kmer);
                if (v) {
                    if (prev_kmer_good) {
                        assert(vertex_cache_.count(prev_kmer) == 1);
                        const auto u = vertex_cache_.at(prev_kmer);
                        add_edge(u, *v, 1, is_forward_strand, *base_quality_itr, false, sample);
                    }
                    prev_kmer_good = true;
                } else {
                    prev_kmer_good = false;
                }
            } else {
                if (prev_kmer_good) {
                    const auto u = vertex_cache_.at(prev_kmer);
                    const auto v = kmer_itr->second;
                    Edge e; bool e_in_graph;
                    std::tie(e, e_in_graph) = boost::edge(u, v, graph_);
                    if (e_in_graph) {
                        increment_weight(e, is_forward_strand, *base_quality_itr, sample);
                    } else {
                        add_edge(u, v, 1, is_forward_strand, *base_quality_itr, false, sample);
                    }
                }
                if (is_reference(kmer_itr->second)) {
                    ref_kmer_itr = std::find(ref_kmer_itr, std::cend(reference_kmers_), kmer);
                    if (ref_kmer_itr != std::cend(reference_kmers_)) {
                        auto next_kmer_begin = std::next(kmer_begin);
                        auto next_kmer_end   = std::next(kmer_end);
                        const auto ref_offset = std::distance(std::cbegin(reference_kmers_), ref_kmer_itr);
                        auto ref_vertex_itr = std::next(std::cbegin(reference_vertices_), ref_offset);
                        auto ref_edge_itr = std::next(std::cbegin(reference_edges_), ref_offset);
                        ++ref_kmer_itr;
                        for (; next_kmer_end <= std::cend(sequence) && ref_kmer_itr < std::cend(reference_kmers_);
                               ++next_kmer_begin, ++next_kmer_end, ++ref_kmer_itr, ++ref_vertex_itr, ++ref_edge_itr) {
                            if (std::equal(next_kmer_begin, next_kmer_end, std::cbegin(*ref_kmer_itr))) {
                                assert(ref_edge_itr != std::cend(reference_edges_));
                                increment_weight(*ref_edge_itr, is_forward_strand, *base_quality_itr, sample);
                            } else {
                                break;
                            }
                        }
                        if (next_kmer_end > std::cend(sequence)) {
                            return;
                        }
                        kmer_begin = std::prev(next_kmer_begin);
                        kmer_end   = std::prev(next_kmer_end);
                        assert(kmer_end <= std::cend(sequence));
                        kmer = Kmer {kmer_begin, kmer_end};
                    }
                }
                prev_kmer_good = true;
            }
            prev_kmer = kmer;
        }
    }
}

std::size_t Assembler::num_kmers() const noexcept
{
    return vertex_cache_.size();
}

bool Assembler::is_empty() const noexcept
{
    return vertex_cache_.empty();
}

bool Assembler::is_acyclic() const
{
    return !(graph_has_trivial_cycle() || graph_has_nontrivial_cycle());
}

void Assembler::remove_nonreference_cycles(bool break_chains)
{
    remove_all_nonreference_cycles(break_chains);
}

namespace {

struct CycleDetector : public boost::default_dfs_visitor
{
    struct CycleDetectedException {};
    explicit CycleDetector(bool allow_self_edges = true) : allow_self_edges_ {allow_self_edges} {}
    
    template <typename Graph>
    void back_edge(typename boost::graph_traits<Graph>::edge_descriptor e, const Graph& g)
    {
        if (boost::source(e, g) != boost::target(e, g) || !allow_self_edges_) {
            throw CycleDetectedException {};
        }
    }
protected:
    const bool allow_self_edges_;
};

template <typename Container>
struct CyclicEdgeDetector : public boost::default_dfs_visitor
{
    CyclicEdgeDetector(Container& result, bool include_self_edges = true)
    : include_self_edges_ {include_self_edges}, result_ {result} {}
    
    template <typename Graph>
    void back_edge(typename boost::graph_traits<Graph>::edge_descriptor e, const Graph& g)
    {
        if (boost::source(e, g) != boost::target(e, g) || include_self_edges_) {
            result_.push_back(e);
        }
    }
protected:
    const bool include_self_edges_;
    Container& result_;
};

} // namespace

std::vector<Assembler::SampleID> Assembler::find_cyclic_samples(const unsigned min_weight) const
{
    const auto index_map = boost::get(&GraphNode::index, graph_);
    std::deque<Edge> cyclic_edges {};
    CyclicEdgeDetector<decltype(cyclic_edges)> vis {cyclic_edges};
    boost::depth_first_search(graph_, boost::visitor(vis).root_vertex(reference_head()).vertex_index_map(index_map));
    std::set<SampleID> cyclic_samples {};
    for (const Edge& e : cyclic_edges) {
        if (!is_reference(e)) {
            for (const auto& p : graph_[e].samples) {
                if (p.second.weight >= min_weight) {
                    cyclic_samples.insert(p.first);
                }
            } 
        }
    }
    return {std::begin(cyclic_samples), std::end(cyclic_samples)};
}

bool Assembler::is_all_reference() const
{
    const auto p = boost::edges(graph_);
    return std::all_of(p.first, p.second, [this] (const Edge& e) { return is_reference(e); });
}

bool Assembler::is_unique_reference() const
{
    return is_reference_unique_path();
}

void Assembler::try_recover_dangling_branches()
{
    const auto p = boost::vertices(graph_);
    std::for_each(p.first, p.second, [this] (const Vertex& v) {
        if (is_dangling_branch(v)) {
            const auto joining_kmer = find_joining_kmer(v);
            if (joining_kmer) {
                add_edge(v, *joining_kmer, 1, 0, false);
            }
        }
    });
}

void Assembler::prune(const unsigned min_weight)
{
    remove_low_weight_edges(min_weight);
}

void Assembler::cleanup()
{
    if (!is_reference_unique_path()) {
        throw NonUniqueReferenceSequence {};
    }
    auto old_size = boost::num_vertices(graph_);
    if (old_size < 2) return;
    assert(is_reference_unique_path());
    remove_disconnected_vertices();
    auto new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) return;
        old_size = new_size;
    }
    assert(is_reference_unique_path());
    remove_vertices_that_cant_be_reached_from(reference_head());
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) return;
        old_size = new_size;
    }
    assert(is_reference_unique_path());
    remove_vertices_past(reference_tail());
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) return;
        old_size = new_size;
    }
    assert(is_reference_unique_path());
    remove_vertices_that_cant_reach(reference_tail());
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) return;
        old_size = new_size;
    }
    assert(is_reference_unique_path());
    prune_reference_flanks();
    assert(is_reference_unique_path());
    if (is_reference_empty()) {
        clear();
        return;
    }
    new_size = boost::num_vertices(graph_);
    assert(new_size != 0);
    assert(!(boost::num_edges(graph_) == 0 && new_size > 1));
    assert(is_reference_unique_path());
    if (new_size != old_size) {
        regenerate_vertex_indices();
    }
}

void Assembler::clear(const std::vector<SampleID>& samples)
{
    const auto can_remove_edge = [&] (const Edge& e) {
        if (is_reference(e)) return false;
        if (is_artificial(e)) return false;
        if (graph_[e].samples.size() > samples.size()) return false;
        const auto includes_sample = [&] (const auto& p) {
            return std::find(std::cbegin(samples), std::cend(samples), p.first) != std::cend(samples);
        };
        return std::all_of(std::cbegin(graph_[e].samples), std::cend(graph_[e].samples), includes_sample);
    };
    boost::remove_edge_if(can_remove_edge, graph_);
    const auto p = boost::edges(graph_);
    std::for_each(p.first, p.second, [&] (const Edge& e) {
        auto& edge = graph_[e];
        for (const auto& sample : samples) {
            const auto sample_itr = edge.samples.find(sample);
            if (sample_itr != std::cend(edge.samples)) {
                edge.weight -= sample_itr->second.weight;
                edge.forward_strand_weight -= sample_itr->second.forward_strand_weight;
                edge.base_quality_sum -= sample_itr->second.base_quality_sum;
                edge.samples.erase(sample_itr);
            }
        }
    });
}

void Assembler::clear()
{
    graph_.clear();
    vertex_cache_.clear();
    reference_kmers_.clear();
    reference_kmers_.shrink_to_fit();
    reference_vertices_.clear();
    reference_vertices_.shrink_to_fit();
    reference_edges_.clear();
    reference_edges_.shrink_to_fit();
}

bool operator<(const Assembler::Variant& lhs, const Assembler::Variant& rhs) noexcept
{
    if (lhs.begin_pos == rhs.begin_pos) {
        if (lhs.ref.size() == rhs.ref.size()) {
            return lhs.alt < rhs.alt;
        }
        return lhs.ref.size() < rhs.ref.size();
    }
    return lhs.begin_pos < rhs.begin_pos;
}

std::deque<Assembler::Variant>
Assembler::extract_variants(const unsigned max_bubbles, const BubbleScoreSetter min_bubble_scorer)
{
    if (is_empty() || is_all_reference()) return {};
    set_all_edge_transition_scores_from(reference_head());
    auto result = extract_bubble_paths(max_bubbles, min_bubble_scorer);
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

void Assembler::write_dot(std::ostream& out) const
{
    const auto vertex_writer = [this] (std::ostream& out, Vertex v) {
        if (is_reference(v)) {
            out << " [shape=box,color=blue]" << std::endl;
        } else {
            out << " [shape=box,color=red]" << std::endl;
        }
        out << " [label=\"" << kmer_of(v) << "\"]" << std::endl;
    };
    const auto edge_writer = [this] (std::ostream& out, Edge e) {
        if (is_reference(e)) {
            out << " [color=blue]" << std::endl;
        } else if (is_artificial(e)) {
            out << " [style=dotted,color=grey]" << std::endl;
        } else {
            out << " [color=red]" << std::endl;
        }
        out << " [label=\"" << graph_[e].weight << "\"]" << std::endl;
    };
    const auto graph_writer = [] (std::ostream& out) {
        out << "rankdir=LR" << std::endl;
    };
    boost::write_graphviz(out, graph_, vertex_writer, edge_writer, graph_writer);
}

// Kmer
Assembler::Kmer::Kmer(SequenceIterator first, SequenceIterator last) noexcept
: first_ {first}
, last_ {last}
, hash_ {boost::hash_range(first_, last_)}
{}

char Assembler::Kmer::front() const noexcept
{
    return *first_;
}

char Assembler::Kmer::back() const noexcept
{
    return *std::prev(last_);
}

Assembler::Kmer::SequenceIterator Assembler::Kmer::begin() const noexcept
{
    return first_;
}

Assembler::Kmer::SequenceIterator Assembler::Kmer::end() const noexcept
{
    return last_;
}

Assembler::Kmer::operator NucleotideSequence() const
{
    return NucleotideSequence {first_, last_};
}

std::size_t Assembler::Kmer::hash() const noexcept
{
    return hash_;
}

bool operator==(const Assembler::Kmer& lhs, const Assembler::Kmer& rhs) noexcept
{
    return std::equal(lhs.first_, lhs.last_, rhs.first_);
}

bool operator<(const Assembler::Kmer& lhs, const Assembler::Kmer& rhs) noexcept
{
    return std::lexicographical_compare(lhs.first_, lhs.last_, rhs.first_, rhs.last_);
}
//
// Assembler private methods
//

std::size_t Assembler::EdgeHash::operator()(const Edge& e) const
{
    std::size_t result {};
    using boost::hash_combine;
    hash_combine(result, assembler_.source_kmer_of(e).hash());
    hash_combine(result, assembler_.target_kmer_of(e).hash());
    return result;
}

void Assembler::insert_reference_into_empty_graph(const NucleotideSequence& sequence)
{
    assert(sequence.size() >= kmer_size());
    vertex_cache_.reserve(sequence.size() + std::pow(4, 5));
    auto kmer_begin = std::cbegin(sequence);
    auto kmer_end   = std::next(kmer_begin, kmer_size());
    reference_kmers_.emplace_back(kmer_begin, kmer_end);
    if (!contains_kmer(reference_kmers_.back())) {
        const auto u = add_vertex(reference_kmers_.back(), true);
        if (!u) {
            throw NonCanonicalReferenceSequence {sequence};
        }
        reference_vertices_.push_back(*u);
    } else {
        reference_vertices_.push_back(vertex_cache_.at(reference_kmers_.back()));
    }
    ++kmer_begin;
    ++kmer_end;
    for (; kmer_end <= std::cend(sequence); ++kmer_begin, ++kmer_end) {
        reference_kmers_.emplace_back(kmer_begin, kmer_end);
        const auto& kmer = reference_kmers_.back();
        if (!contains_kmer(kmer)) {
            const auto v = add_vertex(kmer, true);
            if (v) {
                reference_vertices_.push_back(*v);
                const auto u = vertex_cache_.at(std::crbegin(reference_kmers_)[1]);
                const auto e = add_reference_edge(u, *v);
                reference_edges_.push_back(e);
            } else {
                throw NonCanonicalReferenceSequence {sequence};
            }
        } else {
            const auto u = vertex_cache_.at(std::crbegin(reference_kmers_)[1]);
            const auto v = vertex_cache_.at(kmer);
            reference_vertices_.push_back(v);
            const auto e = add_reference_edge(u, v);
            reference_edges_.push_back(e);
        }
    }
    assert(reference_vertices_.size() == reference_kmers_.size());
    assert(reference_edges_.size() == reference_vertices_.size() - 1);
    reference_kmers_.shrink_to_fit();
    reference_vertices_.shrink_to_fit();
    reference_edges_.shrink_to_fit();
}

void Assembler::insert_reference_into_populated_graph(const NucleotideSequence& sequence)
{
    assert(sequence.size() >= kmer_size());
    assert(reference_kmers_.empty());
    vertex_cache_.reserve(vertex_cache_.size() + sequence.size() + std::pow(4, 5));
    auto kmer_begin = std::cbegin(sequence);
    auto kmer_end   = std::next(kmer_begin, kmer_size());
    reference_kmers_.emplace_back(kmer_begin, kmer_end);
    if (!contains_kmer(reference_kmers_.back())) {
        const auto u = add_vertex(reference_kmers_.back(), true);
        if (!u) {
            throw NonCanonicalReferenceSequence {sequence};
        }
        reference_vertices_.push_back(*u);
    } else {
        set_vertex_reference(reference_kmers_.back());
        reference_vertices_.push_back(vertex_cache_.at(reference_kmers_.back()));
    }
    ++kmer_begin;
    ++kmer_end;
    for (; kmer_end <= std::cend(sequence); ++kmer_begin, ++kmer_end) {
        reference_kmers_.emplace_back(kmer_begin, kmer_end);
        if (!contains_kmer(reference_kmers_.back())) {
            const auto v = add_vertex(reference_kmers_.back(), true);
            if (v) {
                reference_vertices_.push_back(*v);
                const auto u = vertex_cache_.at(std::crbegin(reference_kmers_)[1]);
                const auto e = add_reference_edge(u, *v);
                reference_edges_.push_back(e);
            } else {
                throw NonCanonicalReferenceSequence {sequence};
            }
        } else {
            const auto u = vertex_cache_.at(std::crbegin(reference_kmers_)[1]);
            const auto v = vertex_cache_.at(reference_kmers_.back());
            reference_vertices_.push_back(v);
            set_vertex_reference(v);
            Edge e; bool e_in_graph;
            std::tie(e, e_in_graph) = boost::edge(u, v, graph_);
            if (e_in_graph) {
                set_edge_reference(e);
            } else {
                e = add_reference_edge(u, v);
            }
            reference_edges_.push_back(e);
        }
    }
    vertex_cache_.rehash(vertex_cache_.size());
    reference_kmers_.shrink_to_fit();
    reference_vertices_.shrink_to_fit();
    reference_edges_.shrink_to_fit();
    regenerate_vertex_indices();
    reference_head_position_ = 0;
}

bool Assembler::contains_kmer(const Kmer& kmer) const noexcept
{
    return vertex_cache_.count(kmer) == 1;
}

std::size_t Assembler::count_kmer(const Kmer& kmer) const noexcept
{
    return vertex_cache_.count(kmer);
}

std::size_t Assembler::reference_size() const noexcept
{
    return sequence_length(reference_kmers_.size(), kmer_size());
}

void Assembler::regenerate_vertex_indices()
{
    const auto p = boost::vertices(graph_);
    unsigned idx {0};
    std::for_each(p.first, p.second, [this, &idx] (Vertex v) { graph_[v].index = idx++; });
}

bool Assembler::is_reference_unique_path() const
{
    if (is_reference_empty()) {
        return true;
    } else {
        auto u = reference_head();
        const auto tail = reference_tail();
        const auto is_reference_edge = [this] (const Edge e) { return is_reference(e); };
        while (u != tail) {
            const auto p = boost::out_edges(u, graph_);
            const auto itr = std::find_if(p.first, p.second, is_reference_edge);
            assert(itr != p.second);
            if (std::any_of(boost::next(itr), p.second, is_reference_edge)) {
                return false;
            }
            u = boost::target(*itr, graph_);
        }
        const auto p = boost::out_edges(tail, graph_);
        return std::none_of(p.first, p.second, is_reference_edge);
    }
}

Assembler::Vertex Assembler::null_vertex() const
{
    return boost::graph_traits<KmerGraph>::null_vertex();
}

boost::optional<Assembler::Vertex> Assembler::add_vertex(const Kmer& kmer, const bool is_reference)
{
    if (!utils::is_canonical_dna(kmer)) return boost::none;
    const auto u = boost::add_vertex({boost::num_vertices(graph_), kmer, is_reference}, graph_);
    vertex_cache_.emplace(kmer, u);
    return u;
}

void Assembler::remove_vertex(const Vertex v)
{
    const auto c = vertex_cache_.erase(kmer_of(v));
    assert(c == 1);
    _unused(c); // make production build happy
    boost::remove_vertex(v, graph_);
}

void Assembler::clear_and_remove_vertex(const Vertex v)
{
    const auto c = vertex_cache_.erase(kmer_of(v));
    assert(c == 1);
    _unused(c); // make production build happy
    boost::clear_vertex(v, graph_);
    boost::remove_vertex(v, graph_);
}

void Assembler::clear_and_remove_all(const std::unordered_set<Vertex>& vertices)
{
    for (const Vertex v : vertices) {
        clear_and_remove_vertex(v);
    }
}

Assembler::Edge 
Assembler::add_edge(const Vertex u, const Vertex v,
                    const GraphEdge::WeightType weight, 
                    GraphEdge::WeightType forward_weight,
                    const int base_quality_sum,
                    const bool is_reference,
                    const boost::optional<Assembler::SampleID> sample)
{
    decltype(GraphEdge::samples) samples {};
    if (sample) {
        samples[*sample] = {weight, forward_weight, base_quality_sum};
    }
    return boost::add_edge(u, v, {std::move(samples), weight, forward_weight, base_quality_sum, is_reference}, graph_).first;
}

Assembler::Edge Assembler::add_reference_edge(const Vertex u, const Vertex v)
{
    return add_edge(u, v, 0, 0, 0, true);
}

void Assembler::remove_edge(const Vertex u, const Vertex v)
{
    boost::remove_edge(u, v, graph_);
}

void Assembler::remove_edge(const Edge e)
{
    boost::remove_edge(e, graph_);
}

void Assembler::increment_weight(const Edge e, const bool is_forward, const int base_quality, const SampleID sample)
{
    auto& edge = graph_[e];
    auto sample_itr = edge.samples.find(sample);
    if (sample_itr != std::cend(edge.samples)) {
        ++sample_itr->second.weight;
        sample_itr->second.base_quality_sum += base_quality;
        if (is_forward) ++sample_itr->second.base_quality_sum;
    } else {
        edge.samples[sample] = {1, is_forward, base_quality};
    }
    ++edge.weight;
    edge.base_quality_sum += base_quality;
    if (is_forward) ++edge.forward_strand_weight;
}

void Assembler::set_vertex_reference(const Vertex v)
{
    graph_[v].is_reference = true;
}

void Assembler::set_vertex_reference(const Kmer& kmer)
{
    set_vertex_reference(vertex_cache_.at(kmer));
}

void Assembler::set_edge_reference(const Edge e)
{
    graph_[e].is_reference = true;
}

const Assembler::Kmer& Assembler::kmer_of(const Vertex v) const
{
    return graph_[v].kmer;
}

char Assembler::front_base_of(const Vertex v) const
{
    return kmer_of(v).front();
}

char Assembler::back_base_of(const Vertex v) const
{
    return kmer_of(v).back();
}

const Assembler::Kmer& Assembler::source_kmer_of(const Edge e) const
{
    return kmer_of(boost::source(e, graph_));
}

const Assembler::Kmer& Assembler::target_kmer_of(const Edge e) const
{
    return kmer_of(boost::target(e, graph_));
}

bool Assembler::is_reference(const Vertex v) const
{
    return graph_[v].is_reference;
}

bool Assembler::is_source_reference(const Edge e) const
{
    return is_reference(boost::source(e, graph_));
}

bool Assembler::is_target_reference(const Edge e) const
{
    return is_reference(boost::target(e, graph_));
}

bool Assembler::is_reference(const Edge e) const
{
    return graph_[e].is_reference;
}

bool Assembler::is_artificial(const Edge e) const
{
    return graph_[e].samples.empty();
}

bool Assembler::is_reference_empty() const noexcept
{
    return reference_vertices_.empty();
}

Assembler::Vertex Assembler::reference_head() const
{
    return reference_vertices_.front();
}

Assembler::Vertex Assembler::reference_tail() const
{
    return reference_vertices_.back();
}

Assembler::Vertex Assembler::next_reference(const Vertex u) const
{
    const auto p = boost::out_edges(u, graph_);
    const auto itr = std::find_if(p.first, p.second, [this] (const Edge e) { return is_reference(e); });
    assert(itr != p.second);
    return boost::target(*itr, graph_);
}

Assembler::Vertex Assembler::prev_reference(const Vertex v) const
{
    const auto p = boost::in_edges(v, graph_);
    const auto itr = std::find_if(p.first, p.second, [this] (const Edge e) { return is_reference(e); });
    assert(itr != p.second);
    return boost::source(*itr, graph_);
}

std::size_t Assembler::num_reference_kmers() const
{
    const auto p = boost::vertices(graph_);
    return std::count_if(p.first, p.second, [this] (const Vertex& v) { return is_reference(v); });
}

bool Assembler::is_dangling_branch(const Vertex v) const
{
    return !is_reference(v) && boost::in_degree(v, graph_) > 0 && boost::out_degree(v, graph_) == 0;
}

boost::optional<Assembler::Vertex> Assembler::find_joining_kmer(const Vertex v) const
{
    const auto& kmer = kmer_of(v);
    NucleotideSequence adjacent_kmer {std::next(std::cbegin(kmer)), std::cend(kmer)};
    adjacent_kmer.resize(kmer_size());
    constexpr std::array<NucleotideSequence::value_type, 4> bases {'A', 'C', 'G', 'T'};
    for (const auto base : bases) {
        adjacent_kmer.back() = base;
        const Kmer k {std::cbegin(adjacent_kmer), std::cend(adjacent_kmer)};
        const auto itr = vertex_cache_.find(k);
        if (itr != std::cend(vertex_cache_)) {
            return itr->second;
        }
    }
    return boost::none;
}

Assembler::NucleotideSequence Assembler::make_sequence(const Path& path) const
{
    assert(!path.empty());
    NucleotideSequence result(kmer_size() + path.size() - 1, 'N');
    const auto& first_kmer = kmer_of(path.front());
    auto itr = std::copy(std::cbegin(first_kmer), std::cend(first_kmer), std::begin(result));
    std::transform(std::next(std::cbegin(path)), std::cend(path), itr,
                  [this] (const Vertex v) { return back_base_of(v); });
    return result;
}

Assembler::NucleotideSequence Assembler::make_reference(Vertex from, const Vertex to) const
{
    const auto null = null_vertex();
    
    NucleotideSequence result {};
    if (from == to || from == null) {
        return result;
    }
    auto last = to;
    if (last == null) {
        if (from == reference_tail()) {
            return static_cast<NucleotideSequence>(kmer_of(from));
        }
        last = reference_tail();
    }
    result.reserve(2 * kmer_size());
    const auto& first_kmer = kmer_of(from);
    result.assign(std::cbegin(first_kmer), std::cend(first_kmer));
    from = next_reference(from);
    while (from != last) {
        result.push_back(back_base_of(from));
        from = next_reference(from);
    }
    if (to == null) {
        result.push_back(back_base_of(last));
    }
    
    result.shrink_to_fit();
    
    return result;
}

void Assembler::remove_path(const std::deque<Vertex>& path)
{
    assert(!path.empty());
    if (path.size() == 1) {
        clear_and_remove_vertex(path.front());
    } else {
        remove_edge(*boost::in_edges(path.front(), graph_).first);
        auto prev = path.front();
        std::for_each(std::next(std::cbegin(path)), std::cend(path),
                      [this, &prev] (const Vertex v) {
                          remove_edge(prev, v);
                          remove_vertex(prev);
                          prev = v;
                      });
        remove_edge(*boost::out_edges(path.back(), graph_).first);
        remove_vertex(path.back());
    }
}

bool Assembler::is_bridge(const Vertex v) const
{
    return boost::in_degree(v, graph_) == 1 && boost::out_degree(v, graph_) == 1;
}

bool Assembler::is_reference_bridge(const Vertex v) const
{
    return is_reference(v) && is_bridge(v);
}

Assembler::Path::const_iterator
Assembler::is_bridge_until(const Path::const_iterator first, const Path::const_iterator last) const
{
    return std::find_if_not(first, last, [this] (const Vertex v) { return is_bridge(v); });
}

Assembler::Path::const_iterator Assembler::is_bridge_until(const Path& path) const
{
    return is_bridge_until(std::cbegin(path), std::cend(path));
}

bool Assembler::is_bridge(Path::const_iterator first, Path::const_iterator last) const
{
    return std::all_of(first, last, [this] (const Vertex v) { return is_bridge(v); });
}

bool Assembler::is_bridge(const Path& path) const
{
    return is_bridge(std::cbegin(path), std::cend(path));
}

std::pair<bool, Assembler::Vertex> Assembler::is_bridge_to_reference(Vertex from) const
{
    while (is_bridge(from)) {
        from = *boost::adjacent_vertices(from, graph_).first;
        if (is_reference(from)) {
            return std::make_pair(true, from);
        }
    }
    return std::make_pair(false, null_vertex());
}

bool Assembler::joins_reference_only(const Vertex v) const
{
    return boost::out_degree(v, graph_) == 1 && is_reference(*boost::out_edges(v, graph_).first);
}

bool Assembler::joins_reference_only(Path::const_iterator first, Path::const_iterator last) const
{
    const auto itr = std::find_if(first, last, [this] (Vertex v) { return is_reference(v) || boost::out_degree(v, graph_) != 1; });
    return itr == last || is_reference(*itr);
}

bool Assembler::is_trivial_cycle(const Edge e) const
{
    return boost::source(e, graph_) == boost::target(e, graph_);
}

bool Assembler::graph_has_trivial_cycle() const
{
    const auto p = boost::edges(graph_);
    return std::any_of(p.first, p.second, [this] (const Edge& e) { return is_trivial_cycle(e); });
}

bool Assembler::graph_has_nontrivial_cycle() const
{
    const auto index_map = boost::get(&GraphNode::index, graph_);
    try {
        boost::depth_first_search(graph_, boost::visitor(CycleDetector {}).root_vertex(reference_head()).vertex_index_map(index_map));
        return false;
    } catch (const CycleDetector::CycleDetectedException&) {
        return true;
    }
}

void Assembler::remove_trivial_nonreference_cycles()
{
    boost::remove_edge_if([this] (const Edge e) { return !is_reference(e) && is_trivial_cycle(e); }, graph_);
}

void Assembler::remove_nontrivial_nonreference_cycles()
{
    const auto index_map = boost::get(&GraphNode::index, graph_);
    std::deque<Edge> cyclic_edges {};
    CyclicEdgeDetector<decltype(cyclic_edges)> vis {cyclic_edges, false};
    boost::depth_first_search(graph_, boost::visitor(vis).root_vertex(reference_head()).vertex_index_map(index_map));
    for (const Edge& e : cyclic_edges) {
        if (!(is_reference(e) || is_simple_deletion(e))) {
            remove_edge(e);
        }
    }
}

void Assembler::remove_all_nonreference_cycles(const bool break_chains)
{
    const auto index_map = boost::get(&GraphNode::index, graph_);
    std::deque<Edge> cyclic_edges {};
    CyclicEdgeDetector<decltype(cyclic_edges)> vis {cyclic_edges};
    boost::depth_first_search(graph_, boost::visitor(vis).root_vertex(reference_head()).vertex_index_map(index_map));
    if (cyclic_edges.empty()) return;
    std::unordered_set<Vertex> bad_kmers {}, reference_origins {}, reference_sinks {};
    std::deque<std::pair<Vertex, Vertex>> cyclic_reference_segments {};
    if (break_chains) {
        bad_kmers.reserve(std::min(2 * cyclic_edges.size(), num_kmers()));
        reference_origins.reserve(num_reference_kmers());
    }
    for (const Edge& back_edge : cyclic_edges) {
        if (!is_reference(back_edge)) {
            if (break_chains) {
                Vertex cycle_origin {boost::source(back_edge, graph_)};
                while (!is_reference(cycle_origin) && is_bridge(cycle_origin) && bad_kmers.count(cycle_origin) == 0) {
                    bad_kmers.insert(cycle_origin);
                    cycle_origin = *boost::inv_adjacent_vertices(cycle_origin, graph_).first;
                }
                bool is_reference_origin {false};
                if (is_reference(cycle_origin)) {
                    reference_origins.insert(cycle_origin);
                    is_reference_origin = true;
                } else {
                    bad_kmers.insert(cycle_origin);
                }
                Vertex cycle_sink {boost::target(back_edge, graph_)};
                while (!is_reference(cycle_sink) && is_bridge(cycle_sink) && bad_kmers.count(cycle_sink) == 0) {
                    bad_kmers.insert(cycle_origin);
                    cycle_sink = *boost::adjacent_vertices(cycle_sink, graph_).first;
                }
                if (is_reference(cycle_sink)) {
                    reference_sinks.insert(cycle_sink);
                    if (is_reference_origin) {
                        cyclic_reference_segments.emplace_back(cycle_sink, cycle_origin);
                    } else if (boost::out_degree(cycle_origin, graph_) > 1) {
                        const auto p = boost::out_edges(cycle_origin, graph_);
                        std::vector<Vertex> reference_tails {};
                        reference_tails.reserve(std::distance(p.first, p.second));
                        std::for_each(p.first, p.second, [&] (Edge tail_edge) {
                            if (tail_edge != back_edge) {
                                auto tail = boost::target(tail_edge, graph_);
                                while (!is_reference(tail) && boost::out_degree(tail, graph_) == 1
                                       && bad_kmers.count(tail) == 0) {
                                    bad_kmers.insert(tail);
                                    tail = *boost::adjacent_vertices(tail, graph_).first;
                                }
                                if (is_reference(tail)) {
                                    reference_tails.push_back(tail);
                                }
                            }
                        });
                        if (!reference_tails.empty()) {
                            Vertex cycle_tail;
                            if (reference_tails.size() == 1) {
                                cycle_tail = reference_tails.front();
                            } else {
                                // Just add the rightmost reference vertex
                                auto itr = std::find_first_of(std::crbegin(reference_vertices_), std::crend(reference_vertices_),
                                                              std::cbegin(reference_tails), std::cend(reference_tails));
                                assert(itr != std::crend(reference_vertices_));
                                cycle_tail = *itr;
                            }
                            const std::array<Vertex, 2> cycle_vertices {cycle_sink, cycle_tail};
                            auto itr = std::find_first_of(std::cbegin(reference_vertices_), std::cend(reference_vertices_),
                                                          std::cbegin(cycle_vertices), std::cend(cycle_vertices));
                            assert(itr != std::cend(reference_vertices_));
                            if (*itr == cycle_vertices.front()) {
                                cyclic_reference_segments.emplace_back(cycle_sink, cycle_tail);
                            } else {
                                cyclic_reference_segments.emplace_back(cycle_tail, cycle_sink);
                            }
                        }
                    }
                } else {
                    bad_kmers.insert(cycle_sink);
                }
            }
            remove_edge(back_edge);
        }
    }
    bool regenerate_indices {false};
    for (Vertex v : reference_origins) {
        boost::remove_in_edge_if(v, [this] (Edge e) { return !is_reference(e); }, graph_);
        regenerate_indices = true;
    }
    for (Vertex v : reference_sinks) {
        boost::remove_out_edge_if(v, [this] (Edge e) { return !is_reference(e); }, graph_);
        regenerate_indices = true;
    }
    if (!bad_kmers.empty()) {
        clear_and_remove_all(bad_kmers);
        regenerate_indices = true;
    }
    if (!cyclic_reference_segments.empty()) {
        for (const auto& p : cyclic_reference_segments) {
            const auto first_vertex_itr = std::find(std::cbegin(reference_vertices_), std::cend(reference_vertices_), p.first);
            assert(first_vertex_itr != std::cend(reference_vertices_));
            const auto last_vertex_itr  = std::find(first_vertex_itr, std::cend(reference_vertices_), p.second);
            assert(last_vertex_itr != std::cend(reference_vertices_));
            std::for_each(first_vertex_itr, last_vertex_itr, [this] (Vertex v) {
                boost::remove_in_edge_if(v, [this] (Edge e) { return !is_reference(e); }, graph_);
                boost::remove_out_edge_if(v, [this] (Edge e) { return !is_reference(e); }, graph_);
            });
        }
        regenerate_indices = true;
    }
    if (regenerate_indices) {
        regenerate_vertex_indices();
    }
}

bool Assembler::is_simple_deletion(Edge e) const
{
    return !is_reference(e) && is_source_reference(e) && is_target_reference(e);
}

bool Assembler::is_on_path(const Edge e, const Path& path) const
{
    if (path.size() < 2) return false;
    auto first_vertex = std::cbegin(path);
    auto next_vertex = std::next(first_vertex);
    const auto last_vertex = std::cend(path);
    Edge path_edge; bool good;
    for (; next_vertex != last_vertex; ++first_vertex, ++next_vertex) {
        std::tie(path_edge, good) = boost::edge(*first_vertex, *next_vertex, graph_);
        assert(good);
        if (path_edge == e) return true;
    }
    return false;
}

bool Assembler::connects_to_path(Edge e, const Path& path) const
{
    return e == *boost::in_edges(path.front(), graph_).first || e == *boost::out_edges(path.back(), graph_).first;
}

bool Assembler::is_dependent_on_path(Edge e, const Path& path) const
{
    return connects_to_path(e, path) || is_on_path(e, path);
}

Assembler::GraphEdge::WeightType Assembler::weight(const Path& path) const
{
    if (path.size() < 2) return 0;
    return std::inner_product(std::cbegin(path), std::prev(std::cend(path)),
                              std::next(std::cbegin(path)), GraphEdge::WeightType {0},
                              std::plus<> {},
                              [this] (const auto& u, const auto& v) {
                                  Edge e; bool good;
                                  std::tie(e, good) = boost::edge(u, v, graph_);
                                  assert(good);
                                  return graph_[e].weight;
                              });
}

Assembler::GraphEdge::WeightType Assembler::max_weight(const Path& path, const EdgeSet& scored_edges) const
{
    GraphEdge::WeightType result {0};
    for_each_edge(path, [&] (const Edge e) {
        if (scored_edges.count(e) == 0) {
            result = std::max(result, graph_[e].weight);
        }
    });
    return result;
}

unsigned Assembler::count_low_weights(const Path& path, const unsigned low_weight) const
{
    if (path.size() < 2) return 0;
    return std::inner_product(std::cbegin(path), std::prev(std::cend(path)), std::next(std::cbegin(path)), 0u, std::plus<> {},
                              [this, low_weight] (const auto& u, const auto& v) {
                                  Edge e; bool good;
                                  std::tie(e, good) = boost::edge(u, v, graph_);
                                  assert(good);
                                  return graph_[e].weight <= low_weight ? 1 : 0;
                              });
}

bool Assembler::has_low_weight_flanks(const Path& path, const unsigned low_weight) const
{
    if (path.size() < 2) return false;
    Edge e; bool good;
    std::tie(e, good) = boost::edge(path[0], path[1], graph_);
    assert(good);
    if (graph_[e].weight <= low_weight) return true;
    std::tie(e, good) = boost::edge(std::crbegin(path)[1], std::crbegin(path)[0], graph_);
    assert(good);
    return graph_[e].weight <= low_weight;
}

unsigned Assembler::count_low_weight_flanks(const Path& path, unsigned low_weight) const
{
    if (path.size() < 2) return 0;
    const auto is_low_weight = [this, low_weight] (const auto& u, const auto& v) {
        Edge e; bool good;
        std::tie(e, good) = boost::edge(u, v, graph_);
        assert(good);
        return graph_[e].weight > low_weight ? 1 : 0;
    };
    const auto first_head_high_weight = std::adjacent_find(std::cbegin(path), std::cend(path), is_low_weight);
    const auto first_tail_high_weight = std::adjacent_find(std::crbegin(path), std::make_reverse_iterator(first_head_high_weight),
                                                           [&] (const auto& u, const auto& v) {
                                                               return is_low_weight(v, u);
                                                           });
    const auto num_head_low_weight = std::distance(std::cbegin(path), first_head_high_weight);
    const auto num_tail_low_weight = std::distance(std::crbegin(path), first_tail_high_weight);
    return static_cast<unsigned>(num_head_low_weight + num_tail_low_weight);
}

Assembler::PathWeightStats 
Assembler::compute_weight_stats(const Path& path, const EdgeSet& scored_edges) const
{
    PathWeightStats result {};
    if (path.size() > 1) {
        result.max = max_weight(path, scored_edges);
        result.min = result.max;
        result.distribution.resize(result.max + 1);
        std::vector<GraphEdge::WeightType> weights {};
        weights.reserve(path.size() - 1);
        for_each_edge(path, [&] (const Edge e) {
            if (scored_edges.count(e) == 0) {
                const auto weight = graph_[e].weight;
                ++result.distribution[weight];
                result.total += weight;
                result.min = std::min(result.min, weight);
                if (!is_artificial(e)) {
                    const auto forward_weight = graph_[e].forward_strand_weight;
                    result.total_forward += forward_weight;
                    result.total_reverse += weight - forward_weight;
                }
                weights.push_back(weight);
            }
        });
        if (!weights.empty()) {
            for (auto& w : result.distribution) w /= weights.size();
            result.mean = static_cast<double>(result.total) / weights.size();
            result.median = maths::median(weights);
            result.stdev = maths::stdev(weights);
        }
    }
    return result;
}

Assembler::GraphEdge::WeightType Assembler::sum_source_in_edge_weight(const Edge e) const
{
    const auto p = boost::in_edges(boost::source(e, graph_), graph_);
    using Weight = GraphEdge::WeightType;
    return std::accumulate(p.first, p.second, Weight {0},
                           [this] (const Weight curr, const Edge& e) {
                               return curr + graph_[e].weight;
                           });
}

Assembler::GraphEdge::WeightType Assembler::sum_target_out_edge_weight(const Edge e) const
{
    const auto p = boost::out_edges(boost::target(e, graph_), graph_);
    using Weight = GraphEdge::WeightType;
    return std::accumulate(p.first, p.second, Weight {0},
                           [this] (const Weight curr, const Edge& e) {
                               return curr + graph_[e].weight;
                           });
}

bool Assembler::all_in_edges_low_weight(Vertex v, unsigned min_weight) const
{
    const auto p = boost::in_edges(v, graph_);
    return std::all_of(p.first, p.second, [this, min_weight] (Edge e) { return graph_[e].weight < min_weight; });
}

bool Assembler::all_out_edges_low_weight(Vertex v, unsigned min_weight) const
{
    const auto p = boost::out_edges(v, graph_);
    return std::all_of(p.first, p.second, [this, min_weight] (Edge e) { return graph_[e].weight < min_weight; });
}

bool Assembler::is_low_weight(const Vertex v, const unsigned min_weight) const
{
    return !is_reference(v) && all_in_edges_low_weight(v, min_weight) && all_out_edges_low_weight(v, min_weight);
}

std::size_t Assembler::low_weight_out_degree(Vertex v, unsigned min_weight) const
{
    const auto p = boost::out_edges(v, graph_);
    const auto d = std::count_if(p.first, p.second, [this, min_weight] (Edge e) { return graph_[e].weight < min_weight; });
    return static_cast<std::size_t>(d);
}

std::size_t Assembler::low_weight_in_degree(Vertex v, unsigned min_weight) const
{
    const auto p = boost::in_edges(v, graph_);
    const auto d = std::count_if(p.first, p.second, [this, min_weight] (Edge e) { return graph_[e].weight < min_weight; });
    return static_cast<std::size_t>(d);
}

bool Assembler::is_low_weight_source(Vertex v, unsigned min_weight) const
{
    const auto num_low_weight = low_weight_out_degree(v, min_weight);
    return num_low_weight > 0 && num_low_weight < boost::out_degree(v, graph_);
}

bool Assembler::is_low_weight_sink(Vertex v, unsigned min_weight) const
{
    const auto num_low_weight = low_weight_in_degree(v, min_weight);
    return num_low_weight > 0 && num_low_weight < boost::in_degree(v, graph_);
}

namespace {

template <typename Iterator, typename Set>
bool all_in(Iterator first, Iterator last, const Set& values)
{
    return std::all_of(first, last, [&] (const auto& value) { return values.count(value) == 1; });
}

} // namespace

void Assembler::remove_low_weight_edges(const unsigned min_weight)
{
    boost::remove_edge_if([this, min_weight] (const Edge& e) {
        return !is_reference(e) && graph_[e].weight < min_weight
               && sum_source_in_edge_weight(e) < min_weight
               && sum_target_out_edge_weight(e) < min_weight;
    }, graph_);
}

void Assembler::remove_disconnected_vertices()
{
    VertexIterator vi, vi_end, vi_next;
    std::tie(vi, vi_end) = boost::vertices(graph_);
    for (vi_next = vi; vi != vi_end; vi = vi_next) {
        ++vi_next;
        if (boost::degree(*vi, graph_) == 0) {
            remove_vertex(*vi);
        }
    }
}

std::unordered_set<Assembler::Vertex> Assembler::find_reachable_kmers(const Vertex from) const
{
    std::unordered_set<Vertex> result {};
    result.reserve(boost::num_vertices(graph_));
    auto vis = boost::make_bfs_visitor(boost::write_property(boost::typed_identity_property_map<Vertex>(),
                                                             std::inserter(result, std::begin(result)),
                                                             boost::on_discover_vertex()));
    boost::breadth_first_search(graph_, from,
                                boost::visitor(vis).vertex_index_map(boost::get(&GraphNode::index, graph_)));
    return result;
}

std::deque<Assembler::Vertex> Assembler::remove_vertices_that_cant_be_reached_from(const Vertex v)
{
    const auto reachables = find_reachable_kmers(v);
    VertexIterator vi, vi_end, vi_next;
    std::tie(vi, vi_end) = boost::vertices(graph_);
    std::deque<Vertex> result {};
    for (vi_next = vi; vi != vi_end; vi = vi_next) {
        ++vi_next;
        if (reachables.count(*vi) == 0) {
            result.push_back(*vi);
            clear_and_remove_vertex(*vi);
        }
    }
    return result;
}

void Assembler::remove_vertices_that_cant_reach(const Vertex v)
{
    if (!is_reference_empty()) {
        const auto transpose = boost::make_reverse_graph(graph_);
        const auto index_map = boost::get(&GraphNode::index, transpose);
        std::unordered_set<Vertex> reachables {};
        auto vis = boost::make_bfs_visitor(boost::write_property(boost::typed_identity_property_map<Vertex>(),
                                                                 std::inserter(reachables, std::begin(reachables)),
                                                                 boost::on_discover_vertex()));
        boost::breadth_first_search(transpose, v, boost::visitor(vis).vertex_index_map(index_map));
        VertexIterator vi, vi_end, vi_next;
        std::tie(vi, vi_end) = boost::vertices(graph_);
        for (vi_next = vi; vi != vi_end; vi = vi_next) {
            ++vi_next;
            if (reachables.count(*vi) == 0) {
                clear_and_remove_vertex(*vi);
            }
        }
    }
}

void Assembler::remove_vertices_past(const Vertex v)
{
    auto reachables = find_reachable_kmers(v);
    reachables.erase(v);
    boost::clear_out_edges(v, graph_);
    std::deque<Vertex> cycle_tails {};
    // Must check for cycles that lead back to v
    for (auto u : reachables) {
        Edge e; bool present;
        std::tie(e, present) = boost::edge(u, v, graph_);
        if (present) cycle_tails.push_back(u);
    }
    if (!cycle_tails.empty()) {
        // We can check reachable back edges as the links from v were cut previously
        const auto transpose = boost::make_reverse_graph(graph_);
        const auto index_map = boost::get(&GraphNode::index, transpose);
        std::unordered_set<Vertex> back_reachables {};
        auto vis = boost::make_bfs_visitor(boost::write_property(boost::typed_identity_property_map<Vertex>(),
                                                                 std::inserter(back_reachables, std::begin(back_reachables)),
                                                                 boost::on_discover_vertex()));
        for (auto u : cycle_tails) {
            boost::breadth_first_search(transpose, u, boost::visitor(vis).vertex_index_map(index_map));
            reachables.erase(u);
        }
        // The intersection of reachables & back_reachables are vertices part
        // of a cycle past v. The remaining vertices in reachables are safe to
        // remove.
        bool has_intersects {false};
        for (auto u : back_reachables) {
            const auto iter = reachables.find(u);
            if (iter != std::cend(reachables)) {
                reachables.erase(iter);
                has_intersects = true;
            }
        }
        if (has_intersects) {
            const auto removed = remove_vertices_that_cant_be_reached_from(reference_head());
            for (auto u : removed) reachables.erase(u);
        }
    }
    clear_and_remove_all(reachables);
}

bool Assembler::can_prune_reference_flanks() const
{
    return boost::out_degree(reference_head(), graph_) == 1 || boost::in_degree(reference_tail(), graph_) == 1;
}

void Assembler::pop_reference_head()
{
    reference_kmers_.pop_front();
    reference_vertices_.pop_front();
    if (!reference_edges_.empty()) {
        reference_edges_.pop_front();
    }
    ++reference_head_position_;
}

void Assembler::pop_reference_tail()
{
    reference_kmers_.pop_back();
    reference_vertices_.pop_back();
    if (!reference_edges_.empty()) {
        reference_edges_.pop_back();
    }
}

void Assembler::prune_reference_flanks()
{
    if (!is_reference_empty()) {
        auto new_head_itr = std::cbegin(reference_vertices_);
        const auto is_bridge_vertex = [this] (const Vertex v) { return is_bridge(v); };
        if (boost::in_degree(reference_head(), graph_) == 0 && boost::out_degree(reference_head(), graph_) == 1) {
            new_head_itr = std::find_if_not(std::next(new_head_itr), std::cend(reference_vertices_), is_bridge_vertex);
            std::for_each(std::cbegin(reference_vertices_), new_head_itr, [this] (const Vertex u) {
                remove_edge(u, *boost::adjacent_vertices(u, graph_).first);
                remove_vertex(u);
                pop_reference_head();
            });
        }
        if (new_head_itr != std::cend(reference_vertices_) && boost::in_degree(reference_tail(), graph_) == 1
            && boost::out_degree(reference_tail(), graph_) == 0) {
            const auto new_tail_itr = std::find_if_not(std::next(std::crbegin(reference_vertices_)),
                                                       std::make_reverse_iterator(new_head_itr),
                                                       is_bridge_vertex);
            std::for_each(std::crbegin(reference_vertices_), new_tail_itr, [this] (const Vertex u) {
                remove_edge(*boost::inv_adjacent_vertices(u, graph_).first, u);
                remove_vertex(u);
                pop_reference_tail();
            });
        }
    }
}

Assembler::DominatorMap
Assembler::build_dominator_tree(const Vertex from) const
{
    DominatorMap result;
    result.reserve(boost::num_vertices(graph_));
    boost::lengauer_tarjan_dominator_tree(graph_, from,  boost::make_assoc_property_map(result));
    auto it = std::cbegin(result);
    for (; it != std::cend(result);) {
        if (it->second == null_vertex()) {
            it = result.erase(it);
        } else {
            ++it;
        }
    }
    result.rehash(result.size());
    return result;
}

std::unordered_set<Assembler::Vertex> Assembler::extract_nondominants(const Vertex from) const
{
    const auto dom_tree = build_dominator_tree(from);
    std::unordered_set<Vertex> dominators {};
    dominators.reserve(dom_tree.size());
    for (const auto& p : dom_tree) {
        dominators.emplace(p.second);
    }
    std::unordered_set<Vertex> result {};
    result.reserve(dom_tree.size());
    for (const auto& p : dom_tree) {
        if (dominators.count(p.first) == 0) {
            result.emplace(p.first);
        }
    }
    return result;
}

std::deque<Assembler::Vertex> Assembler::extract_nondominant_reference(const DominatorMap& dominator_tree) const
{
    std::unordered_set<Vertex> dominators {};
    dominators.reserve(dominator_tree.size());
    for (const auto& p : dominator_tree) {
        dominators.emplace(p.second);
    }
    std::deque<Vertex> result {};
    for (const auto& p : dominator_tree) {
        if (is_reference(p.first) && p.first != reference_tail() && dominators.count(p.first) == 0) {
            result.push_back(p.first);
        }
    }
    return result;
}

std::pair<Assembler::Vertex, unsigned> Assembler::find_bifurcation(Vertex from, const Vertex to) const
{
    unsigned count {0};
    while (from != to) {
        const auto d = boost::out_degree(from, graph_);
        if (d == 0 || d > 1) {
            return std::make_pair(from, count);
        }
        from = *boost::adjacent_vertices(from, graph_).first;
        ++count;
    }
    return std::make_pair(from, count);
}

template <typename V, typename G>
auto count_out_weight(const V& v, const G& g)
{
    using T = decltype(g[typename boost::graph_traits<G>::edge_descriptor()].weight);
    const auto p = boost::out_edges(v, g);
    return std::accumulate(p.first, p.second, T {0},
                           [&g] (const auto curr, const auto& e) {
                               return curr + g[e].weight;
                           });
}

template <typename R, typename T>
auto compute_transition_score(const T edge_weight, const T total_out_weight) noexcept
{
    if (total_out_weight == 0) {
        return R {0};
    } else if (edge_weight == 0) {
        return -10 * std::log10(R {1} / total_out_weight);
    } else if (edge_weight == total_out_weight) {
        return R {1} / edge_weight;
    } else {
        return -10 * std::log10(static_cast<R>(edge_weight) / total_out_weight);
    }
}

void Assembler::set_out_edge_transition_scores(const Vertex v)
{
    const auto total_out_weight = count_out_weight(v, graph_);
    const auto p = boost::out_edges(v, graph_);
    using R = GraphEdge::ScoreType;
    std::for_each(p.first, p.second, [this, total_out_weight] (const Edge& e) {
        graph_[e].transition_score = compute_transition_score<R>(graph_[e].weight, total_out_weight);
    });
}

template <typename R, typename G, typename PropertyMap>
struct TransitionScorer : public boost::default_dfs_visitor
{
    using Vertex = typename boost::graph_traits<G>::vertex_descriptor;
    
    PropertyMap map;
    
    TransitionScorer(PropertyMap m) : map {m} {}
    
    void discover_vertex(Vertex v, const G& g)
    {
        const auto total_out_weight = count_out_weight(v, g);
        const auto p = boost::out_edges(v, g);
        std::for_each(p.first, p.second, [this, &g, total_out_weight] (const auto& e) {
            boost::put(map, e, compute_transition_score<R>(g[e].weight, total_out_weight));
        });
    }
};

void Assembler::set_all_edge_transition_scores_from(const Vertex src)
{
    auto score_map = boost::get(&GraphEdge::transition_score, graph_);
    TransitionScorer<GraphEdge::ScoreType, KmerGraph, decltype(score_map)> vis {score_map};
    boost::depth_first_search(graph_, boost::visitor(vis)
                              .vertex_index_map(boost::get(&GraphNode::index, graph_)));
    for (auto p = boost::edges(graph_); p.first != p.second; ++p.first) {
        graph_[*p.first].transition_score = score_map[*p.first];
    }
}

void Assembler::set_all_in_edge_transition_scores(const Vertex v, const GraphEdge::ScoreType score)
{
    const auto p = boost::in_edges(v, graph_);
    std::for_each(p.first, p.second, [this, score] (const Edge e) {
        graph_[e].transition_score = score;
    });
}

void Assembler::block_all_in_edges(const Vertex v)
{
    const auto total_out_weight = count_out_weight(v, graph_);
    set_all_in_edge_transition_scores(v, compute_transition_score<GraphEdge::ScoreType>(0u, total_out_weight));
}

Assembler::PredecessorMap Assembler::find_shortest_scoring_paths(const Vertex from, const bool use_weights) const
{
    assert(from != null_vertex());
    std::unordered_map<Vertex, Vertex> result {};
    result.reserve(boost::num_vertices(graph_));
    if (use_weights) {
        boost::dag_shortest_paths(graph_, from,
                                  boost::weight_map(boost::get(&GraphEdge::weight, graph_))
                                  .predecessor_map(boost::make_assoc_property_map(result))
                                  .vertex_index_map(boost::get(&GraphNode::index, graph_)));
    } else {
        boost::dag_shortest_paths(graph_, from,
                                  boost::weight_map(boost::get(&GraphEdge::transition_score, graph_))
                                  .predecessor_map(boost::make_assoc_property_map(result))
                                  .vertex_index_map(boost::get(&GraphNode::index, graph_)));
    }
    return result;
}

bool Assembler::is_on_path(const Vertex v, const PredecessorMap& predecessors, const Vertex from) const
{
    if (v == from) return true;
    assert(predecessors.count(from) == 1);
    auto itr = predecessors.find(from);
    while (itr != std::end(predecessors) && itr->first != itr->second) {
        if (itr->second == v) return true;
        itr = predecessors.find(itr->second);
    }
    return false;
}

bool Assembler::is_on_path(const Edge e, const PredecessorMap& predecessors, const Vertex from) const
{
    assert(predecessors.count(from) == 1);
    const auto last = std::cend(predecessors);
    auto itr1 = predecessors.find(from);
    assert(itr1 != last);
    auto itr2 = predecessors.find(itr1->second);
    Edge path_edge; bool good;
    while (itr2 != last && itr1 != itr2) {
        std::tie(path_edge, good) = boost::edge(itr2->second, itr1->second, graph_);
        assert(good);
        if (path_edge == e) {
            return true;
        }
        itr1 = itr2;
        itr2 = predecessors.find(itr1->second);
    }
    return false;
}

Assembler::Path Assembler::extract_full_path(const PredecessorMap& predecessors, const Vertex from) const
{
    assert(predecessors.count(from) == 1);
    Path result {from};
    auto itr = predecessors.find(from);
    while (itr != std::end(predecessors) && itr->first != itr->second) {
        result.push_front(itr->second);
        itr = predecessors.find(itr->second);
    }
    return result;
}

std::tuple<Assembler::Vertex, Assembler::Vertex, unsigned>
Assembler::backtrack_until_nonreference(const PredecessorMap& predecessors, Vertex from) const
{
    assert(predecessors.count(from) == 1);
    auto v = predecessors.at(from);
    unsigned count {1};
    const auto head = reference_head();
    while (v != head) {
        assert(from != v); // was not reachable from source
        const auto p = boost::edge(v, from, graph_);
        assert(p.second);
        if (!is_reference(p.first)) break;
        from = v;
        assert(predecessors.count(from) == 1);
        v = predecessors.at(from);
        ++count;
    }
    return std::make_tuple(v, from, count);
}

Assembler::Path Assembler::extract_nonreference_path(const PredecessorMap& predecessors, Vertex from) const
{
    Path result {from};
    from = predecessors.at(from);
    while (!is_reference(from)) {
        result.push_front(from);
        from = predecessors.at(from);
    }
    return result;
}

template <typename Map>
auto count_unreachables(const Map& predecessors)
{
    return std::count_if(std::cbegin(predecessors), std::cend(predecessors),
                         [] (const auto& p) { return p.first == p.second; });
}

template <typename Path, typename Map>
void erase_all(const Path& path, Map& dominator_tree)
{
    for (const auto& v : path) dominator_tree.erase(v);
}

template <typename V, typename BidirectionalIt, typename Map>
bool is_dominated_by_path(const V& vertex, const BidirectionalIt first, const BidirectionalIt last,
                          const Map& dominator_tree)
{
    const auto& dominator = dominator_tree.at(vertex);
    const auto rfirst = std::make_reverse_iterator(last);
    const auto rlast  = std::make_reverse_iterator(first);
    // reverse because more likely to be a closer vertex
    return std::find(rfirst, rlast, dominator) != rlast;
}

Assembler::Edge Assembler::head_edge(const Path& path) const
{
    assert(path.size() > 1);
    Edge e; bool good;
    std::tie(e, good) = boost::edge(path[0], path[1], graph_);
    assert(good);
    return e;
}

int Assembler::head_mean_base_quality(const Path& path) const
{
    const auto& fork_edge = graph_[head_edge(path)];
    return fork_edge.weight > 0 ? fork_edge.base_quality_sum / fork_edge.weight : 0;
}
int Assembler::tail_mean_base_quality(const Path& path) const
{
    if (path.size() < 3) return 0;
    int base_quality_sum {0};
    const auto add_edge = [&] (const auto& u, const auto& v) {
        Edge e; bool good;
        std::tie(e, good) = boost::edge(u, v, graph_);
        assert(good);
        base_quality_sum += graph_[e].base_quality_sum;
        return graph_[e].weight;
    };
    auto total_weight = std::inner_product(std::next(std::cbegin(path)), std::prev(std::cend(path)),
                       std::next(std::cbegin(path), 2), GraphEdge::WeightType {0},
                       std::plus<> {}, add_edge);
    return base_quality_sum / total_weight;
}

double Assembler::get_min_bubble_score(Vertex ref_head, Vertex ref_tail, BubbleScoreSetter min_bubble_scorer) const
{
    const auto ref_head_itr = std::find(std::cbegin(reference_vertices_), std::cend(reference_vertices_), ref_head);
    assert(ref_head_itr != std::cend(reference_vertices_));
    const auto ref_head_idx = static_cast<std::size_t>(std::distance(std::cbegin(reference_vertices_), ref_head_itr));
    const auto ref_tail_itr = std::find(ref_head_itr, std::cend(reference_vertices_), ref_tail);
    assert(ref_tail_itr != std::cend(reference_vertices_));
    const auto ref_tail_idx = static_cast<std::size_t>(std::distance(std::cbegin(reference_vertices_), ref_tail_itr));
    return min_bubble_scorer(ref_head_idx, ref_tail_idx);
}

namespace {

double base_quality_probability(const int base_quality)
{
    // If the given base quality is zero then there were no observations so just report 1
    return base_quality > 0 ? maths::phred_to_probability<>(base_quality) : 1.0;
}

} // namespace

double Assembler::bubble_score(const Path& path, const EdgeSet& exclude) const
{
    if (path.size() < 2) return 0;
    const auto weight_stats = compute_weight_stats(path, exclude);
    auto result = static_cast<double>(weight_stats.max);
    if (params_.use_strand_bias) {
        GraphEdge::WeightType context_forward_weight {0}, context_reverse_weight {0};
        const auto add_strand_weights = [&] (const auto& p) {
            std::for_each(p.first, p.second, [&] (const Edge e) {
                const auto forward_weight = graph_[e].forward_strand_weight;
                context_forward_weight += forward_weight;
                context_reverse_weight += (graph_[e].weight - forward_weight);
            });
        };
        if (boost::in_degree(path.front(), graph_) > 0) {
            add_strand_weights(boost::in_edges(path.front(), graph_));
        } else {
            add_strand_weights(boost::out_edges(path.front(), graph_));
        }
        // path.front() is reference, path.back() is not
        const auto p = boost::adjacent_vertices(path.back(), graph_);
        std::for_each(p.first, p.second, [&] (Vertex v) {
            if (is_reference(v)) {
                if (boost::out_degree(v, graph_) > 0) {
                    add_strand_weights(boost::out_edges(v, graph_));
                } else {
                    add_strand_weights(boost::in_edges(v, graph_));
                }
            }
        });
        const auto context_total_weight = context_forward_weight + context_reverse_weight;
        const auto context_avg_weight = std::max(context_total_weight / 2, 1u);
        const auto approx_vaf = std::min(static_cast<double>(weight_stats.mean) / context_avg_weight, 1.0);
        if (approx_vaf < 0.9) {
            // don't apply to putative homozygous variants, otherwise we could exclude
            // due to misalligned non-variant reads.
            auto pval = maths::fisher_exact_test(weight_stats.total_forward, context_forward_weight,
                                                 weight_stats.total_reverse, context_reverse_weight);
            if (pval < 0.001) result *= pval;
        }
    }
    if (params_.use_base_qualities) {
        result *= base_quality_probability(head_mean_base_quality(path));
        if (path.size() > 2) {
            result *= base_quality_probability(tail_mean_base_quality(path));
        }
    }
    return result;
}

std::vector<Assembler::EdgePath> Assembler::extract_k_shortest_paths(Vertex src, Vertex dst, unsigned k) const
{
    auto weights = boost::get(&GraphEdge::transition_score, graph_);
    auto indices = boost::get(&GraphNode::index, graph_);
    const auto ksps = boost::yen_ksp(graph_, src, dst, std::move(weights), std::move(indices), k);
    std::vector<EdgePath> result {};
    result.reserve(k);
    for (const auto& p : ksps) {
        result.emplace_back(std::cbegin(p.second), std::cend(p.second));
    }
    return result;
}

struct BfsSearcherSuccess {};

template <typename Vertex>
struct BfsSearcher : public boost::default_bfs_visitor
{
    BfsSearcher(Vertex v) : v_ {v} {}
    
    template <typename Graph>
    void discover_vertex(Vertex v, const Graph& g) const
    {
        if (v == v_) throw BfsSearcherSuccess {};
    }
private:
    Vertex v_;
};

template <typename Vertex>
auto make_bfs_searcher(Vertex v)
{
    return BfsSearcher<Vertex> {v};
}

std::deque<Assembler::Variant>
Assembler::extract_bubble_paths(unsigned k, const BubbleScoreSetter min_bubble_scorer)
{
    auto num_remaining_alt_kmers = num_kmers() - num_reference_kmers();
    std::deque<Variant> result {};
    boost::optional<DominatorMap> dominator_tree {};
    bool use_weights {false};
    EdgeSet scored_edges {boost::num_edges(graph_), EdgeHash {*this}};
    while (k > 0 && num_remaining_alt_kmers > 0) {
        auto predecessors = find_shortest_scoring_paths(reference_head(), use_weights);
        assert(count_unreachables(predecessors) == 1);
        Vertex ref, alt; unsigned rhs_kmer_count;
        std::tie(alt, ref, rhs_kmer_count) = backtrack_until_nonreference(predecessors, reference_tail());
        if (alt == reference_head()) {
            // complete reference path is shortest path
            if (dominator_tree) {
                if (use_weights) {
                    utils::append(extract_bubble_paths_with_ksp(k, min_bubble_scorer, scored_edges), result);
                    return result;
                } else {
                    use_weights = true;
                    continue;
                }
            } else {
                dominator_tree = build_dominator_tree(reference_head());
                const auto nondominant_reference = extract_nondominant_reference(*dominator_tree);
                for (Vertex v : nondominant_reference) {
                    block_all_in_edges(v);
                }
                continue;
            }
        }
        bool removed_bubble {false};
        while (alt != reference_head()) {
            auto alt_path = extract_nonreference_path(predecessors, alt);
            assert(!alt_path.empty());
            assert(predecessors.count(alt_path.front()) == 1);
            const auto ref_before_bubble = predecessors.at(alt_path.front());
            auto ref_seq = make_reference(ref_before_bubble, ref);
            alt_path.push_front(ref_before_bubble);
            const auto min_bubble_score = get_min_bubble_score(ref_before_bubble, ref, min_bubble_scorer);
            const auto score = bubble_score(alt_path, scored_edges);
            const bool is_extractable {score >= min_bubble_score};
            auto alt_seq = make_sequence(alt_path);
            alt_path.pop_front();
            rhs_kmer_count += count_kmers(ref_seq, kmer_size());
            if (is_extractable) {
                const auto pos = reference_head_position_ + reference_size() - sequence_length(rhs_kmer_count, kmer_size());
                result.emplace_front(pos, std::move(ref_seq), std::move(alt_seq));
            }
            for_each_edge(alt_path, [&] (const Edge e) { scored_edges.insert(e); });
            --rhs_kmer_count; // because we padded one reference kmer to make ref_seq
            Edge edge_to_alt; bool good;
            std::tie(edge_to_alt, good) = boost::edge(alt, ref, graph_);
            assert(good);
            if (alt_path.size() == 1 && is_simple_deletion(edge_to_alt)) {
                remove_edge(alt_path.front(), ref);
                set_out_edge_transition_scores(alt_path.front());
            } else {
                auto vertex_before_bridge = ref_before_bubble;
                assert(is_reference(vertex_before_bridge));
                const auto bifurication_point_itr = is_bridge_until(alt_path);
                if (bifurication_point_itr == std::cend(alt_path)) {
                    /*
                           -> ref -> ref -> ->
                          /                    \
                       ref                     ref*
                          \                    /
                            -> alt -> alt -> alt
                      
                       The entire alt path can be removed as it never needs to be explored again;
                       the reference path can always be taken.
                    */
                    remove_path(alt_path);
                    regenerate_vertex_indices();
                    set_out_edge_transition_scores(vertex_before_bridge);
                    num_remaining_alt_kmers -= alt_path.size();
                    alt_path.clear();
                    removed_bubble = true;
                } else if (joins_reference_only(bifurication_point_itr, std::cend(alt_path))) {
                    /*
                            -> ref -> ref -> ref -> ref ->
                           /                              \
                        ref -> alt -> -> alt              ref
                           \                 \            /
                            -> alt -> alt -> alt* -> alt ->
                      
                       The alt path can be removed up until the bifurication point as one of the other
                       alt paths can be taken in future paths.
                    */
                    if (bifurication_point_itr == std::cbegin(alt_path)) {
                        remove_edge(ref_before_bubble, alt_path.front());
                        alt_path.clear();
                    } else {
                        alt_path.erase(bifurication_point_itr, std::cend(alt_path));
                        remove_path(alt_path);
                    }
                    regenerate_vertex_indices();
                    set_out_edge_transition_scores(vertex_before_bridge);
                    num_remaining_alt_kmers -= alt_path.size();
                    removed_bubble = true;
                } else if (boost::in_degree(*bifurication_point_itr, graph_) == 1) {
                    const auto next_bifurication_point_itr = is_bridge_until(std::next(bifurication_point_itr), std::cend(alt_path));
                    if (next_bifurication_point_itr != std::cend(alt_path)) {
                        if (boost::out_degree(*next_bifurication_point_itr, graph_) == 1) {
                            if (joins_reference_only(next_bifurication_point_itr, std::cend(alt_path))) {
                                const auto p = boost::adjacent_vertices(*bifurication_point_itr, graph_);
                                auto is_simple_bubble = std::all_of(p.first, p.second, [&] (Vertex v) {
                                    if (v == *std::next(bifurication_point_itr)) {
                                        return true;
                                    } else {
                                        try {
                                            boost::breadth_first_search(graph_, v,
                                                                        boost::visitor(make_bfs_searcher(*next_bifurication_point_itr)).
                                                                        vertex_index_map(boost::get(&GraphNode::index, graph_)));
                                        } catch (const BfsSearcherSuccess&) {
                                            return true;
                                        }
                                        return false;
                                    }
                                });
                                if (is_simple_bubble) {
                                    /*
                                            -> ref -> ref -> ref -> ref ->
                                           /                               \
                                        ref         -> alt -> -> alt        ref
                                           \      /                  \     /
                                            -> alt* -> alt -> alt -> alt ->
                                    */
                                    const auto bifurication_point = *bifurication_point_itr;
                                    if (std::next(bifurication_point_itr) == next_bifurication_point_itr) {
                                        remove_edge(*bifurication_point_itr, *std::next(bifurication_point_itr));
                                        alt_path.clear();
                                    } else {
                                        alt_path.erase(std::cbegin(alt_path), std::next(bifurication_point_itr));
                                        alt_path.erase(next_bifurication_point_itr, std::cend(alt_path));
                                        remove_path(alt_path);
                                    }
                                    regenerate_vertex_indices();
                                    set_all_edge_transition_scores_from(bifurication_point);
                                    num_remaining_alt_kmers -= alt_path.size();
                                    removed_bubble = true;
                                }
                            }
                        }
                    } else {
                        /*
                                -> ref -> ref -> ref -> ref ->
                               /                        /      \
                            ref         -> alt -> -> alt        ref
                               \      /                        /
                                -> alt* -> alt -> alt -> alt ->
                        */
                        const auto bifurication_point = *bifurication_point_itr;
                        alt_path.erase(std::cbegin(alt_path), std::next(bifurication_point_itr));
                        if (alt_path.empty()) {
                            remove_edge(bifurication_point, ref);
                        } else {
                            remove_path(alt_path);
                        }
                        regenerate_vertex_indices();
                        set_all_edge_transition_scores_from(bifurication_point);
                        num_remaining_alt_kmers -= alt_path.size();
                        removed_bubble = true;
                    }
                }
            }
            unsigned kmer_count_to_alt;
            std::tie(alt, ref, kmer_count_to_alt) = backtrack_until_nonreference(predecessors, ref_before_bubble);
            rhs_kmer_count += kmer_count_to_alt;
            if (!use_weights && k > 0) --k;
        }
        if (!removed_bubble && use_weights) {
            if (can_prune_reference_flanks()) {
                prune_reference_flanks();
                regenerate_vertex_indices();
            }
            utils::append(extract_bubble_paths_with_ksp(k, min_bubble_scorer, scored_edges), result);
            return result;
        } else if (!removed_bubble) {
            use_weights = true;
        }
        assert(boost::out_degree(reference_head(), graph_) > 0);
        assert(boost::in_degree(reference_tail(), graph_) > 0);
        if (can_prune_reference_flanks()) {
            prune_reference_flanks();
            regenerate_vertex_indices();
        }
    }
    return result;
}

std::deque<Assembler::SubGraph> Assembler::find_independent_subgraphs() const
{
    assert(!reference_vertices_.empty());
    const auto diverges  = [this] (const Vertex& v) { return boost::out_degree(v, graph_) > 1; };
    const auto coalesces = [this] (const Vertex& v) { return boost::in_degree(v, graph_) > 1; };
    auto subgraph_head_itr = std::find_if(std::cbegin(reference_vertices_), std::cend(reference_vertices_), diverges);
    if (subgraph_head_itr == std::cend(reference_vertices_)) {
        return {{reference_head(), reference_tail(), 0}};
    }
    auto candidate_subgraph_tail_itr = std::find_if(subgraph_head_itr, std::cend(reference_vertices_), coalesces);
    if (candidate_subgraph_tail_itr == std::cend(reference_vertices_)) {
        return {{reference_head(), reference_tail(), 0}};
    } else {
        std::deque<Vertex> coalescent_points {*candidate_subgraph_tail_itr};
        std::copy_if(std::next(candidate_subgraph_tail_itr), std::cend(reference_vertices_),
                     std::front_inserter(coalescent_points), coalesces);
        const auto dominator = build_dominator_tree(reference_head());
        std::deque<SubGraph> result {};
        while (subgraph_head_itr != std::cend(reference_vertices_)) {
            auto subbgraph_end = std::find_if(std::cbegin(coalescent_points), std::cend(coalescent_points),
                                              [&] (const Vertex& v) { return dominator.at(v) == *subgraph_head_itr; });
            assert(subbgraph_end != std::cend(coalescent_points));
            auto subgraph_offset = static_cast<std::size_t>(std::distance(std::cbegin(reference_vertices_), subgraph_head_itr));
            result.push_back({*subgraph_head_itr, *subbgraph_end, subgraph_offset});
            subgraph_head_itr = std::find_if(std::find(subgraph_head_itr, std::cend(reference_vertices_), *subbgraph_end),
                                             std::cend(reference_vertices_), diverges);
            coalescent_points.erase(subbgraph_end, std::cend(coalescent_points));
        }
        return result;
    }
}

std::deque<Assembler::Variant> 
Assembler::extract_bubble_paths_with_ksp(const unsigned k, const BubbleScoreSetter min_bubble_scorer, EdgeSet& scored_edges)
{
    const auto subgraphs = find_independent_subgraphs();
    std::deque<Variant> result {};
    for (const auto& subgraph : subgraphs) {
        auto shortest_paths = extract_k_shortest_paths(subgraph.head, subgraph.tail, k);
        for (const auto& path : shortest_paths) {
            assert(!path.empty());
            const auto is_alt_edge = [this] (const Edge e) { return !is_reference(e); };
            auto alt_head_itr = std::find_if(std::cbegin(path), std::cend(path), is_alt_edge);
            auto lhs_kmer_count = std::distance(std::cbegin(path), alt_head_itr);
            while (alt_head_itr != std::cend(path)) {
                const auto ref_before_bubble = boost::source(*alt_head_itr, graph_);
                assert(is_reference(ref_before_bubble));
                const auto alt_tail_itr = std::find_if(alt_head_itr, std::cend(path), [this] (Edge e) { return is_target_reference(e); });
                assert(alt_tail_itr != std::cend(path));
                const auto ref_after_bubble = boost::target(*alt_tail_itr, graph_);
                assert(!is_reference(*alt_tail_itr));
                assert(is_reference(ref_after_bubble));
                auto ref_seq = make_reference(ref_before_bubble, ref_after_bubble);
                Path alt_path {};
                std::transform(alt_head_itr, std::next(alt_tail_itr), std::back_inserter(alt_path),
                               [this] (Edge e) { return boost::source(e, graph_); });
                const auto num_ref_kmers = count_kmers(ref_seq, kmer_size());
                const auto min_bubble_score = get_min_bubble_score(ref_before_bubble, ref_after_bubble, min_bubble_scorer);
                const auto score = bubble_score(alt_path, scored_edges);
                if (score >= min_bubble_score) {
                    auto alt_seq = make_sequence(alt_path);
                    const auto pos = reference_head_position_ + subgraph.reference_offset + lhs_kmer_count;
                    result.emplace_front(pos, std::move(ref_seq), std::move(alt_seq));
                }
                for_each_edge(alt_path, [&] (const Edge e) { scored_edges.insert(e); });
                alt_head_itr = std::find_if(std::next(alt_tail_itr), std::cend(path), is_alt_edge);
                lhs_kmer_count += num_ref_kmers + std::distance(alt_tail_itr, alt_head_itr) - 1;
            }
        }
    }
    return result;
}

// debug

std::ostream& operator<<(std::ostream& os, const Assembler::Kmer& kmer)
{
    std::copy(std::cbegin(kmer), std::cend(kmer), std::ostreambuf_iterator<char> {os});
    return os;
}

void Assembler::print_reference_head() const
{
    std::cout << "reference head is " << kmer_of(reference_head()) << std::endl;
}

void Assembler::print_reference_tail() const
{
    std::cout << "reference tail is " << kmer_of(reference_tail()) << std::endl;
}

void Assembler::print_reference_path() const
{
    print(reference_vertices_);
}

void Assembler::print(const Edge e) const
{
    std::cout << kmer_of(boost::source(e, graph_)) << "->" << kmer_of(boost::target(e, graph_));
}

void Assembler::print(const Path& path) const
{
    assert(!path.empty());
    std::transform(std::cbegin(path), std::prev(std::cend(path)), std::ostream_iterator<Kmer> {std::cout, "->"},
                   [this] (const Vertex v) { return kmer_of(v); });
    std::cout << kmer_of(path.back());
}

void Assembler::print_verbose(const Path& path) const
{
    if (path.size() < 2) return;
    std::transform(std::cbegin(path), std::prev(std::cend(path)), std::next(std::cbegin(path)),
                   std::ostream_iterator<std::string> {std::cout, "->"},
                   [this] (const auto& u, const auto& v) {
                       Edge e; bool good;
                       std::tie(e, good) = boost::edge(u, v, graph_);
                       assert(good);
                       auto result = static_cast<std::string>(this->kmer_of(v));
                       result += "(";
                       result += std::to_string(graph_[e].weight);
                       result += ", " + std::to_string(graph_[e].forward_strand_weight);
                       result += ", " + std::to_string(graph_[e].base_quality_sum);
                       result += ")";
                       return result;
                   });
    std::cout << kmer_of(path.back());
}

void Assembler::print_dominator_tree() const
{
    const auto dom_tree = build_dominator_tree(reference_head());
    for (const auto& p : dom_tree) {
        std::cout << kmer_of(p.first) << " dominated by " << kmer_of(p.second) << std::endl;
    }
}

// non-member methods

bool operator==(const Assembler::Variant& lhs, const Assembler::Variant& rhs) noexcept
{
    return lhs.begin_pos == rhs.begin_pos && lhs.ref.size() == rhs.ref.size() && lhs.alt == rhs.alt;
}

} // namespace coretools
} // namespace octopus
