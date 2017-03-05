// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "assembler.hpp"

#include <iterator>
#include <algorithm>
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

#include "ksp/yen_ksp.hpp"

#include "utils/sequence_utils.hpp"
#include "utils/append.hpp"

#include "timers.hpp"

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

Assembler::BadReferenceSequence::BadReferenceSequence(NucleotideSequence reference_sequence)
: std::invalid_argument {"bad reference sequence"}
, reference_sequence_ {std::move(reference_sequence)}
{}

const char* Assembler::BadReferenceSequence::what() const noexcept
{
    return std::invalid_argument::what();
}

Assembler::Assembler(const unsigned kmer_size)
: k_ {kmer_size}
, reference_kmers_ {}
, reference_head_position_ {0}
, reference_vertices_ {}
{}

Assembler::Assembler(const unsigned kmer_size, const NucleotideSequence& reference)
: k_ {kmer_size}
, reference_kmers_ {}
, reference_head_position_ {0}
, reference_vertices_ {}
{
    insert_reference_into_empty_graph(reference);
}

unsigned Assembler::kmer_size() const noexcept
{
    return k_;
}

void Assembler::insert_reference(const NucleotideSequence& sequence)
{
    if (sequence.size() >= k_) {
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

void Assembler::insert_read(const NucleotideSequence& sequence)
{
    if (sequence.size() >= k_) {
        auto kmer_begin = std::cbegin(sequence);
        auto kmer_end   = std::next(kmer_begin, k_);
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
                   ++next_kmer_begin, ++next_kmer_end, ++ref_kmer_itr, ++ref_vertex_itr, ++ref_edge_itr) {
                if (std::equal(next_kmer_begin, next_kmer_end, std::cbegin(*ref_kmer_itr))) {
                    assert(ref_edge_itr != std::cend(reference_edges_));
                    increment_weight(*ref_edge_itr);
                } else {
                    break;
                }
            }
            if (next_kmer_end > std::cend(sequence)) {
                return;
            }
            if (next_kmer_begin != std::next(kmer_begin)) {
                kmer_begin = std::prev(next_kmer_begin);
                kmer_end   = std::prev(next_kmer_end);
                if (kmer_end <= std::cend(sequence)) {
                    prev_kmer = Kmer {kmer_begin, kmer_end};
                }
            }
        } else {
            ++kmer_begin;
            ++kmer_end;
        }
        for (; kmer_end <= std::cend(sequence); ++kmer_begin, ++kmer_end) {
            Kmer kmer {kmer_begin, kmer_end};
            const auto kmer_itr = vertex_cache_.find(kmer);
            if (kmer_itr == std::cend(vertex_cache_)) {
                const auto v = add_vertex(kmer);
                if (v) {
                    if (prev_kmer_good) {
                        const auto u = vertex_cache_.at(prev_kmer);
                        add_edge(u, *v, 1);
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
                        increment_weight(e);
                    } else {
                        add_edge(u, v, 1);
                    }
                } else {
                    if (is_reference(kmer_itr->second)) {
                        ref_kmer_itr = std::find(ref_kmer_itr, std::cend(reference_kmers_), kmer);
                        assert(ref_kmer_itr != std::cend(reference_kmers_));
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
                                increment_weight(*ref_edge_itr);
                            } else {
                                break;
                            }
                        }
                        if (next_kmer_end > std::cend(sequence)) {
                            return;
                        }
                        if (next_kmer_begin != std::next(kmer_begin)) {
                            kmer_begin = std::prev(next_kmer_begin);
                            kmer_end   = std::prev(next_kmer_end);
                            if (kmer_end <= std::cend(sequence)) {
                                kmer = Kmer {kmer_begin, kmer_end};
                            }
                        }
                    }
                    prev_kmer_good = true;
                }
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

template <typename G>
struct CycleDetector : public boost::default_dfs_visitor
{
    using Edge = typename boost::graph_traits<G>::edge_descriptor;
    CycleDetector(bool& is_acyclic) : is_acyclic_ {is_acyclic} {}
    void back_edge(Edge e, const G& g)
    {
        if (boost::source(e, g) != boost::target(e, g)) {
            is_acyclic_ = false;
        }
    }
protected:
    bool& is_acyclic_;
};

bool Assembler::is_acyclic() const
{
    if (graph_has_trivial_cycle()) return false;
    bool is_acyclic {true};
    CycleDetector<KmerGraph> vis {is_acyclic};
    const auto index_map = boost::get(&GraphNode::index, graph_);
    boost::depth_first_search(graph_, boost::visitor(vis).vertex_index_map(index_map));
    return is_acyclic;
}

bool Assembler::is_all_reference() const
{
    const auto p = boost::edges(graph_);
    return std::all_of(p.first, p.second, [this] (const Edge& e) { return is_reference(e); });
}

bool Assembler::prune(const unsigned min_weight)
{
    if (!is_reference_unique_path()) {
        clear();
        return false;
    }
    auto old_size = boost::num_vertices(graph_);
    if (old_size < 2) return true;
    
    remove_trivial_nonreference_cycles();
    auto new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) {
            return true;
        }
        old_size = new_size;
    }
    assert(is_reference_unique_path());
    remove_low_weight_edges(min_weight);
    remove_disconnected_vertices();
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) {
            return true;
        }
        old_size = new_size;
    }
    assert(is_reference_unique_path());
    remove_vertices_that_cant_be_reached_from(reference_head());
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) return true;
        old_size = new_size;
    }
    assert(is_reference_unique_path());
    remove_vertices_past(reference_tail());
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) return true;
        old_size = new_size;
    }
    assert(is_reference_unique_path());
    remove_vertices_that_cant_reach(reference_tail());
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) return true;
        old_size = new_size;
    }
    assert(is_reference_unique_path());
    try {
        if (can_prune_reference_flanks()) {
            prune_reference_flanks();
        }
    } catch (const boost::not_a_dag& e) {
        clear();
        return false;
    }
    assert(is_reference_unique_path());
    if (is_reference_empty()) {
        clear();
        return true;
    }
    if (can_prune_reference_flanks()) {
        // something is wrong, have seen cases, bug?
        clear();
        return false;
    }
    new_size = boost::num_vertices(graph_);
    assert(new_size != 0);
    assert(!(boost::num_edges(graph_) == 0 && new_size > 1));
    assert(is_reference_unique_path());
    if (new_size != old_size) {
        regenerate_vertex_indices();
        old_size = new_size;
    }
    
    return true;
}

void Assembler::clear()
{
    graph_.clear();
    vertex_cache_.clear();
    reference_kmers_.clear();
    reference_kmers_.shrink_to_fit();
    reference_vertices_.clear();
    reference_kmers_.shrink_to_fit();
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
Assembler::extract_variants(const unsigned max_bubbles, const double min_bubble_score)
{
    if (is_empty() || is_all_reference()) return {};
    set_all_edge_transition_scores_from(reference_head());
    auto result = extract_bubble_paths(max_bubbles, min_bubble_score);
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

void Assembler::write_dot(std::ostream& out) const
{
    const auto vertex_writer = [this] (std::ostream& out, Vertex v) {
        if (is_reference(v)) {
            out << " [color=blue]" << std::endl;
        } else {
            out << " [color=red]" << std::endl;
        }
        out << " [label=\"" << kmer_of(v) << "\"]" << std::endl;
    };
    const auto edge_writer = [this] (std::ostream& out, Edge e) {
        if (is_reference(e)) {
            out << " [color=blue]" << std::endl;
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

// Assembler private types
Assembler::GraphEdge::GraphEdge(const GraphEdge::WeightType weight, const bool is_reference)
: weight {weight}
, is_reference {is_reference}
{}

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
void Assembler::insert_reference_into_empty_graph(const NucleotideSequence& sequence)
{
    assert(sequence.size() >= k_);
    vertex_cache_.reserve(sequence.size() + std::pow(4, 5));
    auto kmer_begin = std::cbegin(sequence);
    auto kmer_end   = std::next(kmer_begin, k_);
    reference_kmers_.emplace_back(kmer_begin, kmer_end);
    if (!contains_kmer(reference_kmers_.back())) {
        const auto u = add_vertex(reference_kmers_.back(), true);
        if (!u) {
            throw BadReferenceSequence {sequence};
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
                throw BadReferenceSequence {sequence};
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
    assert(sequence.size() >= k_);
    assert(reference_kmers_.empty());
    vertex_cache_.reserve(vertex_cache_.size() + sequence.size() + std::pow(4, 5));
    auto kmer_begin = std::cbegin(sequence);
    auto kmer_end   = std::next(kmer_begin, k_);
    reference_kmers_.emplace_back(kmer_begin, kmer_end);
    if (!contains_kmer(reference_kmers_.back())) {
        const auto u = add_vertex(reference_kmers_.back(), true);
        if (!u) {
            throw BadReferenceSequence {sequence};
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
                throw BadReferenceSequence {sequence};
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
    return sequence_length(reference_kmers_.size(), k_);
}

void Assembler::regenerate_vertex_indices()
{
    const auto p = boost::vertices(graph_);
    unsigned idx {0};
    std::for_each(p.first, p.second, [this, &idx] (Vertex v) { graph_[v].index = idx++; });
}

bool Assembler::is_reference_unique_path() const
{
    auto u = reference_head();
    const auto tail = reference_tail();
    const auto is_reference_edge = [this] (const Edge e) { return is_reference(e); };
    
    while (u != tail) {
        const auto p = boost::out_edges(u, graph_);
        const auto it = std::find_if(p.first, p.second, is_reference_edge);
        assert(it != p.second);
        if (std::any_of(boost::next(it), p.second, is_reference_edge)) {
            return false;
        }
        u = boost::target(*it, graph_);
    }
    
    const auto p = boost::out_edges(tail, graph_);
    return std::none_of(p.first, p.second, is_reference_edge);
}

Assembler::Vertex Assembler::null_vertex() const
{
    return boost::graph_traits<KmerGraph>::null_vertex();
}

boost::optional<Assembler::Vertex> Assembler::add_vertex(const Kmer& kmer, const bool is_reference)
{
    if (!utils::is_canonical_dna(kmer)) return boost::none;
    const auto u = boost::add_vertex(GraphNode {boost::num_vertices(graph_), kmer, is_reference}, graph_);
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

Assembler::Edge Assembler::add_edge(const Vertex u, const Vertex v,
                                    const GraphEdge::WeightType weight,
                                    const bool is_reference)
{
    return boost::add_edge(u, v, GraphEdge {weight, is_reference}, graph_).first;
}

Assembler::Edge Assembler::add_reference_edge(const Vertex u, const Vertex v)
{
    return add_edge(u, v, 0, true);
}

void Assembler::remove_edge(const Vertex u, const Vertex v)
{
    boost::remove_edge(u, v, graph_);
}

void Assembler::remove_edge(const Edge e)
{
    boost::remove_edge(e, graph_);
}

void Assembler::increment_weight(const Edge e)
{
    ++graph_[e].weight;
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

bool Assembler::is_reference_empty() const noexcept
{
    return reference_kmers_.empty();
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

Assembler::NucleotideSequence Assembler::make_sequence(const Path& path) const
{
    assert(!path.empty());
    NucleotideSequence result(k_ + path.size() - 1, 'N');
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
    result.reserve(2 * k_);
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

bool Assembler::joins_reference_only(const Vertex v) const
{
    if (boost::out_degree(v, graph_) != 1) {
        return false;
    }
    return is_reference(*boost::out_edges(v, graph_).first);
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
    return e == *boost::in_edges(path.front(), graph_).first
           ||  e == *boost::out_edges(path.back(), graph_).first;
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

void Assembler::remove_trivial_nonreference_cycles()
{
    boost::remove_edge_if([this] (const Edge e) {
        return !is_reference(e) && is_trivial_cycle(e);
    }, graph_);
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

bool Assembler::is_low_weight(const Edge e, const unsigned min_weight) const
{
    if (is_reference(e)) return false;
    const auto edge_weight = graph_[e].weight;
    if (edge_weight >= min_weight) return false;
    const auto source_weight = sum_source_in_edge_weight(e);
    if (source_weight < min_weight) return true;
    const auto target_weight = sum_target_out_edge_weight(e);
    return (source_weight + edge_weight + target_weight) < 3 * min_weight;
}

void Assembler::remove_low_weight_edges(const unsigned min_weight)
{
    boost::remove_edge_if([this, min_weight] (const Edge& e) {
        return is_low_weight(e, min_weight);
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
                                boost::visitor(vis)
                                .vertex_index_map(boost::get(&GraphNode::index, graph_)));
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
    if (is_reference_empty()) return;
    
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
    reference_edges_.pop_back();
}

void Assembler::prune_reference_flanks()
{
    if (is_reference_empty()) return;
    auto new_head_itr = std::cbegin(reference_vertices_);
    const auto is_bridge_vertex = [this] (const Vertex v) { return is_bridge(v); };
    if (boost::out_degree(reference_head(), graph_) == 1) {
        new_head_itr = std::find_if_not(std::next(new_head_itr), std::cend(reference_vertices_), is_bridge_vertex);
        std::for_each(std::cbegin(reference_vertices_), new_head_itr, [this] (const Vertex u) {
            remove_edge(u, *boost::adjacent_vertices(u, graph_).first);
            remove_vertex(u);
            pop_reference_head();
        });
    }
    if (new_head_itr != std::cend(reference_vertices_) && boost::in_degree(reference_tail(), graph_) == 1) {
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
auto compute_transition_score(const T edge_weight, const T total_out_weight, const R max_score = 100) noexcept
{
    if (edge_weight == 0) {
        return R {0};
    } else if (edge_weight == 0) {
        return max_score;
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

Assembler::PredecessorMap Assembler::find_shortest_scoring_paths(const Vertex from) const
{
    assert(from != null_vertex());
    std::unordered_map<Vertex, Vertex> result {};
    result.reserve(boost::num_vertices(graph_));
    boost::dag_shortest_paths(graph_, from,
                              boost::weight_map(boost::get(&GraphEdge::transition_score, graph_))
                              .predecessor_map(boost::make_assoc_property_map(result))
                              .vertex_index_map(boost::get(&GraphNode::index, graph_)));
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

double Assembler::bubble_score(const Path& path) const
{
    const auto path_weight = weight(path);
    const auto num_low_weights = count_low_weights(path, 1);
    const auto num_low_weight_flanks = count_low_weight_flanks(path, 1);
    auto score = static_cast<double>(path_weight) / path.size();
    score /= std::max((num_low_weights - num_low_weight_flanks) + 2 * num_low_weight_flanks, 1u);
    return score;
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

std::deque<Assembler::Variant>
Assembler::extract_bubble_paths(unsigned k, const double min_bubble_score)
{
    auto num_remaining_alt_kmers = num_kmers() - num_reference_kmers();
    std::deque<Variant> result {};
    while (k > 0 && num_remaining_alt_kmers > 0) {
        auto predecessors = find_shortest_scoring_paths(reference_head());
        assert(count_unreachables(predecessors) == 1);
        
        Vertex ref, alt; unsigned rhs_kmer_count;
        std::tie(alt, ref, rhs_kmer_count) = backtrack_until_nonreference(predecessors, reference_tail());
        
        if (alt == reference_head()) {
            // complete reference path is shortest path
            utils::append(extract_bubble_paths_with_ksp(k, min_bubble_score), result);
            return result;
        }
        bool removed_bubble {false};
        while (alt != reference_head()) {
            auto alt_path = extract_nonreference_path(predecessors, alt);
            assert(!alt_path.empty());
            assert(predecessors.count(alt_path.front()) == 1);
            const auto ref_before_bubble = predecessors.at(alt_path.front());
            auto ref_seq = make_reference(ref_before_bubble, ref);
            alt_path.push_front(ref_before_bubble);
            const auto extractable = bubble_score(alt_path) >= min_bubble_score;
            auto alt_seq = make_sequence(alt_path);
            alt_path.pop_front();
            rhs_kmer_count += count_kmers(ref_seq, k_);
            if (extractable) {
                const auto pos = reference_head_position_ + reference_size() - sequence_length(rhs_kmer_count, k_);
                result.emplace_front(pos, std::move(ref_seq), std::move(alt_seq));
            }
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
                const auto bifurication_point = is_bridge_until(alt_path);
                if (bifurication_point == std::cend(alt_path)) {
                    remove_path(alt_path);
                    regenerate_vertex_indices();
                    set_out_edge_transition_scores(vertex_before_bridge);
                    num_remaining_alt_kmers -= alt_path.size();
                    alt_path.clear();
                    removed_bubble = true;
                } else if (joins_reference_only(*bifurication_point)) {
                    alt_path.erase(bifurication_point, std::cend(alt_path));
                    remove_path(alt_path);
                    regenerate_vertex_indices();
                    set_out_edge_transition_scores(vertex_before_bridge);
                    num_remaining_alt_kmers -= alt_path.size();
                    removed_bubble = true;
                }
            }
            unsigned kmer_count_to_alt;
            std::tie(alt, ref, kmer_count_to_alt) = backtrack_until_nonreference(predecessors, ref_before_bubble);
            rhs_kmer_count += kmer_count_to_alt;
            if (k > 0) --k;
        }
        if (!removed_bubble) {
            utils::append(extract_bubble_paths_with_ksp(k, min_bubble_score), result);
            return result;
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

std::deque<Assembler::Variant> Assembler::extract_bubble_paths_with_ksp(const unsigned k, const double min_bubble_score)
{
    auto shortest_paths = extract_k_shortest_paths(reference_head(), reference_tail(), k);
    std::deque<Variant> result {};
    for (const auto& path : shortest_paths) {
        assert(!path.empty());
        assert(boost::source(path.front(), graph_) == reference_head());
        assert(boost::target(path.back(), graph_) == reference_tail());
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
            const auto num_ref_kmers = count_kmers(ref_seq, k_);
            if (bubble_score(alt_path) >= min_bubble_score) {
                auto alt_seq = make_sequence(alt_path);
                const auto pos = reference_head_position_ + lhs_kmer_count;
                result.emplace_front(pos, std::move(ref_seq), std::move(alt_seq));
            }
            alt_head_itr = std::find_if(std::next(alt_tail_itr), std::cend(path), is_alt_edge);
            lhs_kmer_count += num_ref_kmers + std::distance(alt_tail_itr, alt_head_itr) - 1;
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

void Assembler::print_weighted(const Path& path) const
{
    if (path.size() < 2) return;
    std::transform(std::cbegin(path), std::prev(std::cend(path)), std::next(std::cbegin(path)),
                   std::ostream_iterator<std::string> {std::cout, "->"},
                   [this] (const auto& u, const auto& v) {
                       Edge e; bool good;
                       std::tie(e, good) = boost::edge(u, v, graph_);
                       assert(good);
                       return static_cast<std::string>(kmer_of(v)) + "(" + std::to_string(graph_[e].weight) + ")";
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
    return lhs.begin_pos == rhs.begin_pos && lhs.alt == rhs.alt;
}

namespace debug {

void print(const Assembler& assembler)
{
    for (auto ep = boost::edges(assembler.graph_); ep.first != ep.second; ++ep.first) {
        const auto e = *ep.first;
        assembler.print(e);
        std::cout << " weight = " << assembler.graph_[e].weight << " ";
        if (assembler.is_reference(e)) {
            std::cout << "ref";
        } else {
            std::cout << "alt";
        }
        std::cout << " (";
        if (assembler.is_source_reference(e)) {
            std::cout << "ref";
        } else {
            std::cout << "alt";
        }
        std::cout << ",";
        if (assembler.is_target_reference(e)) {
            std::cout << "ref";
        } else {
            std::cout << "alt";
        }
        std::cout << ")" << '\n';
    }
}

}

} // namespace coretools
} // namespace octopus
