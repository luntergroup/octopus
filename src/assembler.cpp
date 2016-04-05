//
//  assembler.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "assembler.hpp"

#include <iterator>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <limits>
#include <cassert>

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

#include <iostream>

namespace
{
    template <typename I>
    auto sequence_length(const I num_kmers, const unsigned kmer_size)
    {
        return num_kmers + kmer_size - 1;
    }
    
    template <typename T>
    auto count_kmers(const T& sequence, const unsigned kmer_size)
    {
        if (sequence.size() < kmer_size) {
            return std::size_t {0};
        } else {
            return sequence.size() - kmer_size + 1;
        }
    }
}

// public methods

Assembler::Variant::Variant(std::size_t pos, SequenceType ref, SequenceType alt)
:
begin_pos {pos},
ref {std::move(ref)},
alt {std::move(alt)}
{}

Assembler::Assembler(const unsigned kmer_size)
:
k_ {kmer_size},
reference_kmers_ {},
reference_head_position_ {0}
{}

Assembler::Assembler(const unsigned kmer_size,
                     const SequenceType& reference)
:
k_ {kmer_size},
reference_kmers_ {},
reference_head_position_ {0}
{
    if (!reference_kmers_.empty()) {
        throw std::runtime_error {"Assembler: only one reference sequence can be inserted into the graph"};
    }
    
    if (reference.size() < k_) {
        throw std::runtime_error {"Assembler:: reference length must >= kmer_size"};
    }
    
    vertex_cache_.reserve(reference.size() + std::pow(4, 5));
    
    auto it1 = std::cbegin(reference);
    auto it2 = std::next(it1, k_);
    
    reference_kmers_.emplace_back(it1, it2);
    
    if (!contains_kmer(reference_kmers_.back())) {
        add_vertex(reference_kmers_.back(), true);
    }
    
    ++it1;
    ++it2;
    
    for (; it2 <= std::cend(reference); ++it1, ++it2) {
        reference_kmers_.emplace_back(it1, it2);
        
        if (!contains_kmer(reference_kmers_.back())) {
            const auto v = add_vertex(reference_kmers_.back(), true);
            if (v) {
                const auto u = vertex_cache_.at(reference_kmers_.crbegin()[1]);
                add_edge(u, *v, 0, true);
            }
        } else {
            const auto u = vertex_cache_.at(reference_kmers_.crbegin()[1]);
            const auto v = vertex_cache_.at(reference_kmers_.back());
            add_edge(u, v, 0, true);
        }
    }
    
    reference_kmers_.shrink_to_fit();
}

void Assembler::insert_reference(const SequenceType& sequence)
{
    if (!reference_kmers_.empty()) {
        throw std::runtime_error {"Assembler: only one reference sequence can be inserted into the graph"};
    }
    
    if (sequence.size() < k_) {
        throw std::runtime_error {"Assembler:: reference length must >= kmer_size"};
    }
    
    vertex_cache_.reserve(vertex_cache_.size() + sequence.size() + std::pow(4, 5));
    
    auto it1 = std::cbegin(sequence);
    auto it2 = std::next(it1, k_);
    
    reference_kmers_.emplace_back(it1, it2);
    
    if (!contains_kmer(reference_kmers_.back())) {
        const auto u = add_vertex(reference_kmers_.back(), true);
        if (!u) {
            throw std::runtime_error {"Assembler: bad reference"};
        }
    } else {
        set_vertex_reference(reference_kmers_.back());
    }
    
    ++it1;
    ++it2;
    
    for (; it2 <= std::cend(sequence); ++it1, ++it2) {
        reference_kmers_.emplace_back(it1, it2);
        
        if (!contains_kmer(reference_kmers_.back())) {
            const auto v = add_vertex(reference_kmers_.back(), true);
            if (v) {
                const auto u = vertex_cache_.at(reference_kmers_.crbegin()[1]);
                add_edge(u, *v, 1, true);
            } else {
                throw std::runtime_error {"Assembler: bad reference"};
            }
        } else {
            set_vertex_reference(reference_kmers_.back());
            const auto u = vertex_cache_.at(reference_kmers_.crbegin()[1]);
            const auto v = vertex_cache_.at(reference_kmers_.back());
            
            const auto p = boost::edge(u, v, graph_);
            
            if (!p.second) {
                add_edge(u, v, 0, true);
            } else {
                set_edge_reference(p.first);
            }
        }
    }
    
    vertex_cache_.rehash(vertex_cache_.size());
    reference_kmers_.shrink_to_fit();
    
    regenerate_vertex_indices();
    
    reference_head_position_ = 0;
}

void Assembler::insert_read(const SequenceType& sequence)
{
    if (sequence.size() < k_) {
        return;
    }
    
    auto it1 = std::cbegin(sequence);
    auto it2 = std::next(it1, k_);
    
    Kmer prev_kmer {it1, it2};
    
    bool prev_kmer_good {true};
    
    if (!contains_kmer(prev_kmer)) {
        const auto u = add_vertex(prev_kmer);
        if (!u) {
            prev_kmer_good = false;
        }
    }
    
    ++it1;
    ++it2;
    
    SequenceType kmer;
    kmer.resize(k_);
    
    for (; it2 <= std::cend(sequence); ++it1, ++it2) {
        std::copy(it1, it2, std::begin(kmer));
        
        if (!contains_kmer(kmer)) {
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
        } else if (prev_kmer_good) {
            const auto u = vertex_cache_.at(prev_kmer);
            const auto v = vertex_cache_.at(kmer);
            
            Edge e;
            bool found;
            std::tie(e, found) = boost::edge(u, v, graph_);
            
            if (found) {
                increment_weight(e);
            } else {
                add_edge(u, v, 1);
            }
        } else {
            prev_kmer_good = true;
        }
        
        prev_kmer = kmer;
    }
}

std::size_t Assembler::num_kmers() const noexcept
{
    return vertex_cache_.size();
}

bool Assembler::empty() const noexcept
{
    return vertex_cache_.empty();
}

template <typename G>
struct CycleDetector : public boost::default_dfs_visitor
{
    using Edge = typename boost::graph_traits<G>::edge_descriptor;
    explicit CycleDetector(bool& is_acyclic) : is_acyclic_ {is_acyclic} {}
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
    if (graph_has_trivial_cycle()) {
        return false;
    }
    
    bool is_acyclic {true};
    
    CycleDetector<KmerGraph> vis {is_acyclic};
    
    const auto index_map = boost::get(&GraphNode::index, graph_);
    
    boost::depth_first_search(graph_, boost::visitor(vis).vertex_index_map(index_map));
    
    return is_acyclic;
}

bool Assembler::is_all_reference() const
{
    const auto p = boost::edges(graph_);
    return std::all_of(p.first, p.second,
                       [this] (const Edge& e) { return is_reference(e); });
}

void Assembler::remove_trivial_nonreference_cycles()
{
    boost::remove_edge_if([this] (const Edge e) {
        return !is_reference(e) && is_trivial_cycle(e);
    }, graph_);
}

bool Assembler::prune(const unsigned min_weight)
{
    auto old_size = boost::num_vertices(graph_);
    
    remove_low_weight_edges(min_weight);
    
    remove_disconnected_vertices();
    
    auto new_size = boost::num_vertices(graph_);
    
    if (new_size != old_size) {
        regenerate_vertex_indices();
        old_size = new_size;
    }
    
    remove_vertices_that_cant_be_reached_from(reference_head());
    
    new_size = boost::num_vertices(graph_);
    
    if (new_size != old_size) {
        regenerate_vertex_indices();
        old_size = new_size;
    }
    
    remove_vertices_past_reference_tail();
    
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        old_size = new_size;
    }
    
    remove_vertices_that_cant_reach(reference_tail());
    
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        old_size = new_size;
    }
    
    try {
        prune_reference_flanks();
    } catch (boost::not_a_dag& e) {
        clear();
        return false;
    }
    
    if (is_reference_empty()) {
        clear();
        return true;
    }
    
    new_size = boost::num_vertices(graph_);
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
}

std::deque<Assembler::Variant> Assembler::extract_variants(unsigned max)
{
    std::deque<Variant> result {};
    
    if (empty() || is_all_reference()) {
        return result;
    }
    
    set_all_edge_log_probabilities_from(reference_head());
    
    extract_highest_probability_bubbles(result);
    
    // TODO
//    while (max > 0 && !is_all_reference()) {
//        const auto old_size = result.size();
//        extract_highest_probability_bubbles(result);
//        if (result.size() == old_size) break;
//        --max;
//    }
    
    return result;
}

// private methods

Assembler::GraphEdge::GraphEdge(const unsigned weight, bool is_reference)
:
weight {weight},
is_reference {is_reference}
{}

bool Assembler::contains_kmer(const Kmer& kmer) const noexcept
{
    return vertex_cache_.count(kmer) > 0;
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
    unsigned i {0};
    std::for_each(p.first, p.second, [this, &i] (Vertex v) { graph_[v].index = i++; });
}

template <typename T>
bool is_dna(const T& sequence)
{
    return std::all_of(std::cbegin(sequence), std::cend(sequence),
                       [] (const char base) {
                           return base == 'A' || base == 'C' || base == 'G' || base == 'T';
                       });
}

Assembler::Vertex Assembler::null_vertex() const
{
    return boost::graph_traits<KmerGraph>::null_vertex();
}

boost::optional<Assembler::Vertex> Assembler::add_vertex(const Kmer& kmer, const bool is_reference)
{
    if (!is_dna(kmer)) {
        return boost::none;
    }
    const auto next_index = static_cast<unsigned>(boost::num_vertices(graph_));
    const auto u = boost::add_vertex(GraphNode {next_index, kmer, is_reference}, graph_);
    vertex_cache_.emplace(kmer, u);
    return u;
}

void Assembler::remove_vertex(const Vertex v)
{
    const auto c = vertex_cache_.erase(kmer_of(v));
    assert(c == 1);
    boost::remove_vertex(v, graph_);
}

void Assembler::clear_and_remove_vertex(const Vertex v)
{
    //std::cout << "clearing vertex " << graph_[v].kmer << std::endl;
    const auto c = vertex_cache_.erase(kmer_of(v));
    assert(c == 1);
    boost::clear_vertex(v, graph_);
    boost::remove_vertex(v, graph_);
}

void Assembler::add_edge(const Vertex u, const Vertex v,
                         const unsigned weight, const bool is_reference)
{
    boost::add_edge(u, v, GraphEdge {weight, is_reference}, graph_);
}

void Assembler::remove_edge(Vertex u, Vertex v)
{
    boost::remove_edge(u, v, graph_);
}

void Assembler::remove_edge(Edge e)
{
    boost::remove_edge(e, graph_);
}

void Assembler::increment_weight(Edge e)
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
    return vertex_cache_.at(reference_kmers_.front());
}

Assembler::Vertex Assembler::reference_tail() const
{
    return vertex_cache_.at(reference_kmers_.back());
}

Assembler::Vertex Assembler::next_reference(const Vertex u) const
{
    const auto p = boost::out_edges(u, graph_);
    const auto it = std::find_if(p.first, p.second,
                                 [this] (const Edge e) { return is_reference(e); });
    assert(it != p.second);
    return boost::target(*it, graph_);
}

Assembler::Vertex Assembler::prev_reference(const Vertex v) const
{
    const auto p = boost::in_edges(v, graph_);
    const auto it = std::find_if(p.first, p.second,
                                 [this] (const Edge e) { return is_reference(e); });
    assert(it != p.second);
    return boost::source(*it, graph_);
}

Assembler::SequenceType Assembler::make_reference(Vertex from, const Vertex to) const
{
    SequenceType result {};
    
    const auto null = null_vertex();
    
    if (from == to || from == null) {
        return result;
    }
    
    auto last = to;
    
    if (last == null) {
        if (from == reference_tail()) {
            return kmer_of(from);
        }
        last = reference_tail();
    }
    
    result.reserve(2 * k_);
    
    result.insert(std::end(result), std::cbegin(kmer_of(from)), std::cend(kmer_of(from)));
    
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

bool Assembler::is_trivial_cycle(const Edge e) const
{
    return boost::source(e, graph_) == boost::target(e, graph_);
}

bool Assembler::graph_has_trivial_cycle() const
{
    const auto p = boost::edges(graph_);
    return std::any_of(p.first, p.second,
                       [this] (const Edge& e) {
                           return is_trivial_cycle(e);
                       });
}

bool Assembler::is_bridge(const Vertex v) const
{
    return boost::in_degree(v, graph_) == 1 && boost::out_degree(v, graph_) == 1;
}

void Assembler::remove_low_weight_edges(const unsigned min_weight)
{
    boost::remove_edge_if([this, min_weight] (const Edge& e) {
        return !is_reference(e) && graph_[e].weight < min_weight;
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

void Assembler::remove_vertices_that_cant_be_reached_from(const Vertex v)
{
    const auto reachables = find_reachable_kmers(v);
    
    VertexIterator vi, vi_end, vi_next;
    
    std::tie(vi, vi_end) = boost::vertices(graph_);
    for (vi_next = vi; vi != vi_end; vi = vi_next) {
        ++vi_next;
        if (reachables.count(*vi) == 0) {
            clear_and_remove_vertex(*vi);
        }
    }
}

void Assembler::remove_vertices_that_cant_reach(const Vertex v)
{
    if (is_reference_empty()) {
        return;
    }
    
    const auto transpose = boost::make_reverse_graph(graph_);
    
    const auto index_map = boost::get(&GraphNode::index, transpose);
    
    std::unordered_set<Vertex> reachables {};
    
    auto vis = boost::make_bfs_visitor(
                       boost::write_property(boost::typed_identity_property_map<Vertex>(),
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

void Assembler::remove_vertices_past_reference_tail()
{
    auto reachables = find_reachable_kmers(reference_tail());
    
    reachables.erase(reference_tail());
    
    for (const Vertex v : reachables) {
        clear_and_remove_vertex(v);
    }
}

void Assembler::prune_reference_flanks()
{
    if (is_reference_empty()) {
        return;
    }
    
    std::deque<Vertex> sorted_vertices {};
    
    boost::topological_sort(graph_, std::front_inserter(sorted_vertices),
                            boost::vertex_index_map(boost::get(&GraphNode::index, graph_)));
    
    if (!is_reference(sorted_vertices.back())) {
        std::cout << reference_kmers_.back() << " " << kmer_of(sorted_vertices.back()) << std::endl;
    }
    
    assert(is_reference(sorted_vertices.front()));
    assert(is_reference(sorted_vertices.back()));
    
    const auto it = std::find_if_not(std::cbegin(sorted_vertices),
                                     std::cend(sorted_vertices),
                                     [this] (const Vertex v) {
                                         return boost::out_degree(v, graph_) == 1
                                         && is_reference(*boost::out_edges(v, graph_).first);
                                     });
    
    std::for_each(std::cbegin(sorted_vertices), it,
                  [this] (const Vertex u) {
                      const auto v = *boost::adjacent_vertices(u, graph_).first;
                      remove_edge(u, v);
                      remove_vertex(u);
                      reference_kmers_.pop_front();
                      ++reference_head_position_;
                  });
    
    const auto it2 = std::find_if_not(std::crbegin(sorted_vertices), std::make_reverse_iterator(it),
                                      [this] (const Vertex v) {
                                          return boost::in_degree(v, graph_) == 1
                                          && is_reference(*boost::in_edges(v, graph_).first);
                                      });
    
    std::for_each(std::crbegin(sorted_vertices), it2,
                  [this] (const Vertex u) {
                      const auto v = *boost::inv_adjacent_vertices(u, graph_).first;
                      remove_edge(v, u);
                      remove_vertex(u);
                      reference_kmers_.pop_back();
                  });
}

std::unordered_map<Assembler::Vertex, Assembler::Vertex>
Assembler::build_dominator_tree(const Vertex from) const
{
    std::unordered_map<Vertex, Vertex> dom_tree;
    dom_tree.reserve(boost::num_vertices(graph_));
    
    auto dom_tree_pred_map = boost::make_assoc_property_map(dom_tree);
    
    boost::lengauer_tarjan_dominator_tree(graph_, from, dom_tree_pred_map);
    
    auto it = std::begin(dom_tree);
    
    for (; it != std::end(dom_tree);) {
        if (it->second == null_vertex()) {
            it = dom_tree.erase(it);
        } else {
            ++it;
        }
    }
    
    dom_tree.rehash(dom_tree.size());
    
    return dom_tree;
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

std::unordered_set<Assembler::Vertex>
Assembler::extract_nondominants_on_path(const Path& path) const
{
    assert(!path.empty());
    
    std::cout << "reachable kmers from " << kmer_of(path.front()) << std::endl;
    for (auto v : find_reachable_kmers(path.front())) {
        std::cout << kmer_of(v) << " ";
    }
    std::cout << std::endl;
    
    auto dt = build_dominator_tree(path.front());
    std::cout << dt.size() << std::endl;
    std::cout << "dominator tree from " << kmer_of(path.front()) << std::endl;
    for (auto p : dt) {
        std::cout << kmer_of(p.first) << " dominated by " << kmer_of(p.second) << std::endl;
    }
    
    const auto nondominants_from_path_head = extract_nondominants(path.front());
    
    std::cout << "nondominants from " << kmer_of(path.front()) << std::endl;
    for (auto v : nondominants_from_path_head) {
        std::cout << kmer_of(v) << " ";
    }
    std::cout << std::endl;
    exit(0);
    
    std::unordered_set<Vertex> result {};
    result.reserve(path.size());
    
    std::copy_if(std::cbegin(path), std::cend(path),
                 std::inserter(result, std::begin(result)),
                 [&] (const Vertex v) {
                     return nondominants_from_path_head.count(v) == 1;
                 });
    
    return result;
}

std::deque<Assembler::Vertex> Assembler::find_nondominant_reference(const Vertex from) const
{
    const auto dom_tree = build_dominator_tree(from);
    
    std::unordered_set<Vertex> dominators {};
    dominators.reserve(dom_tree.size());
    
    for (const auto& p : dom_tree) {
        dominators.emplace(p.second);
    }
    
    std::deque<Vertex> result {};
    
    for (const auto& p : dom_tree) {
        if (is_reference(p.first) && dominators.count(p.first) == 0) {
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
    const auto p = boost::out_edges(v, g);
    return std::accumulate(p.first, p.second, 0,
                           [&g] (const auto curr, const auto& e) {
                               return curr + g[e].weight;
                           });
}

void Assembler::set_out_edge_log_probabilities(const Vertex v)
{
    const auto total_out_weight = count_out_weight(v, graph_);
    
    const auto p = boost::out_edges(v, graph_);
    
    std::for_each(p.first, p.second, [this, total_out_weight] (const Edge& e) {
        graph_[e].neg_log_probability =  std::abs(std::log(static_cast<double>(graph_[e].weight)
                                                           / total_out_weight));
    });
}

template <typename G, typename PropertyMap>
struct ProbabilitySetter : public boost::default_dfs_visitor
{
    using Vertex = typename boost::graph_traits<G>::vertex_descriptor;
    
    PropertyMap map;
    
    explicit ProbabilitySetter(PropertyMap m) : map {m} {}
    
    void discover_vertex(Vertex v, const G& g)
    {
        const auto total_out_weight = count_out_weight(v, g);
        
        const auto p = boost::out_edges(v, g);
        
        std::for_each(p.first, p.second, [this, &g, total_out_weight] (const auto& e) {
            boost::put(map, e, std::abs(std::log(static_cast<double>(g[e].weight)
                                                 / total_out_weight)));
        });
    }
};

void Assembler::set_all_edge_log_probabilities_from(const Vertex src)
{
    const auto index_map = boost::get(&GraphNode::index, graph_);
    
    auto prob_map = boost::get(&GraphEdge::neg_log_probability, graph_);
    
    ProbabilitySetter<KmerGraph, decltype(prob_map)> vis {prob_map};
    
    boost::depth_first_search(graph_, boost::visitor(vis).vertex_index_map(index_map));
    
    for (auto p = boost::edges(graph_); p.first != p.second; ++p.first) {
        graph_[*p.first].neg_log_probability = prob_map[*p.first];
    }
}

void Assembler::set_nondominant_reference_paths_impossible(const Vertex from)
{
    const auto nondominant_ref = find_nondominant_reference(from);
    
    Edge prev;
    
    for (const auto v : nondominant_ref) {
        const auto p = boost::in_edges(v, graph_);
        std::for_each(p.first, p.second, [this, &prev] (const Edge e) {
            if (is_source_reference(e)) {
                prev = e;
            }
            graph_[e].neg_log_probability = 100;
        });
        const auto p2 = boost::out_edges(v, graph_);
        std::for_each(p2.first, p2.second, [this] (const Edge e) {
            graph_[e].neg_log_probability = 100;
        });
        
        // while (true) {
        //
        // }
    }
}

void Assembler::clear_and_remove_all(const std::unordered_set<Vertex>& vertices)
{
    for (const Vertex v : vertices) {
        clear_and_remove_vertex(v);
    }
}

Assembler::SequenceType Assembler::make_sequence(const Path& path) const
{
    SequenceType result {};
    result.reserve(k_ + path.size() - 1);
    
    result.insert(std::end(result),
                  std::cbegin(kmer_of(path.front())),
                  std::cend(kmer_of(path.front())));
    
    std::for_each(std::next(std::cbegin(path)), std::cend(path),
                  [this, &result] (const Vertex v) {
                      result.push_back(back_base_of(v));
                  });
    
    return result;
}

bool Assembler::is_bridge(const Path& path) const
{
    return std::all_of(std::cbegin(path), std::cend(path),
                       [this] (const Vertex v) { return is_bridge(v); });
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

Assembler::PredecessorMap Assembler::find_shortest_paths(const Vertex from) const
{
    assert(from != null_vertex());
    
    std::unordered_map<Vertex, Vertex> preds {};
    preds.reserve(boost::num_vertices(graph_));
    
    boost::dag_shortest_paths(graph_, from,
                              boost::weight_map(boost::get(&GraphEdge::neg_log_probability, graph_))
                              .predecessor_map(boost::make_assoc_property_map(preds))
                              .vertex_index_map(boost::get(&GraphNode::index, graph_)));
    
    return preds;
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
        if (is_reference(p.first)) {
            break;
        }
        from = v;
        assert(predecessors.count(from) == 1);
        v = predecessors.at(from);
        ++count;
    }
    
    return std::make_tuple(v, from, count);
}

Assembler::Path
Assembler::extract_nonreference_path(const PredecessorMap& predecessors, Vertex from) const
{
    assert(!is_reference(from));
    
    Path result {};
    
    result.push_front(from);
    
    from = predecessors.at(from);
    
    while (!is_reference(from)) {
        result.push_front(from);
        from = predecessors.at(from);
    }
    
    return result;
}

void Assembler::extract_highest_probability_bubbles(std::deque<Variant>& result)
{
    std::cout << "reference = " << make_reference(reference_head(), null_vertex()) << std::endl;
    
    std::cout << "head = " << kmer_of(reference_head()) << std::endl;
    
    auto u = reference_head();
    while (boost::out_degree(u, graph_) > 0) {
        u = next_reference(u);
    }
    std::cout << kmer_of(u) << std::endl;
    exit(0);
    
    const auto predecessors = find_shortest_paths(reference_head());
    
    std::cout << std::count_if(std::cbegin(predecessors), std::cend(predecessors),
                              [] (const auto& p) {
                                  return p.first == p.second;
                              }) << std::endl;
    exit(0);
    
    Vertex ref, alt;
    unsigned rhs_kmer_count;
    
    std::tie(alt, ref, rhs_kmer_count) = backtrack_until_nonreference(predecessors, reference_tail());
    
    std::cout << "ref = " << kmer_of(ref) << std::endl;
    std::cout << "alt = " << kmer_of(alt) << std::endl;
    std::cout << "ref->reference_tail = " << make_reference(ref, null_vertex()) << std::endl;
    std::cout << "rhs_kmer_count = " << rhs_kmer_count << std::endl;
    
    if (alt == reference_head()) {
        const auto p = boost::edges(graph_);
        bool b {false};
        std::for_each(p.first, p.second, [this, &b] (const Edge e) {
            if (is_reference(e)) {
                if (graph_[e].neg_log_probability < 1000) b = true;
                graph_[e].neg_log_probability = 1000;
            }
        });
        if (b) {
            extract_highest_probability_bubbles(result);
        }
        
        return;
    }
    
    while (alt != reference_head()) {
        auto alt_path = extract_nonreference_path(predecessors, alt);
        
//        std::cout << "alt_path = ";
//        print(alt_path);
//        std::cout << std::endl;
        
        const auto ref_before_bubble = predecessors.at(alt_path.front());
        
        std::cout << "ref_before_bubble = " << kmer_of(ref_before_bubble) << std::endl;
        
        auto ref_seq = make_reference(next_reference(ref_before_bubble), ref);
        auto alt_seq = make_sequence(alt_path);
        
        rhs_kmer_count += count_kmers(ref_seq, k_);
        
        std::cout << "ref_seq = " << ref_seq << std::endl;
        std::cout << "ref_seq.size() = " << ref_seq.size() << std::endl;
        std::cout << "alt_seq = " << alt_seq << std::endl;
        std::cout << "alt_seq.size() = " << alt_seq.size() << std::endl;
        
//        std::cout << "ref = " << ref_seq << std::endl;
//        std::cout << "alt = " << alt_seq << std::endl;
        
        if (is_bridge(alt_path)) {
            remove_path(alt_path);
            regenerate_vertex_indices();
            set_out_edge_log_probabilities(alt_path.front());
        } else {
            // TEST
            //add_edge(vertex_cache_.at("ATA"), vertex_cache_.at("AAG"), 1);
            auto nondominants = extract_nondominants_on_path(alt_path);
            exit(0);
            clear_and_remove_all(nondominants);
        }
        
        const auto pos = reference_head_position_ + reference_size() - sequence_length(rhs_kmer_count, k_);
        std::cout << "pos = " << pos << std::endl;
        result.emplace_front(pos, std::move(ref_seq), std::move(alt_seq));
        
        unsigned kmer_count_to_alt;
        std::tie(alt, ref, kmer_count_to_alt) = backtrack_until_nonreference(predecessors, ref_before_bubble);
        
        rhs_kmer_count += kmer_count_to_alt;
        
        std::cout << "ref = " << kmer_of(ref) << std::endl;
        std::cout << "alt = " << kmer_of(alt) << std::endl;
        std::cout << "ref->ref_before_bubble = " << make_reference(ref, next_reference(ref_before_bubble)) << std::endl;
        std::cout << "kmer_count_to_alt = " << kmer_count_to_alt << std::endl;
        std::cout << "rhs_kmer_count = " << rhs_kmer_count << std::endl;
    }
    
    const auto old_size = boost::num_vertices(graph_);
    prune_reference_flanks();
    const auto new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
    }
}

void Assembler::print(const Edge e) const
{
    std::cout << kmer_of(boost::source(e, graph_)) << "->" << kmer_of(boost::target(e, graph_));
}

void Assembler::print(const Path& path) const
{
    assert(!path.empty());
    std::transform(std::cbegin(path), std::prev(std::cend(path)),
                   std::ostream_iterator<Kmer> {std::cout, "->"},
                   [this] (const Vertex v) {
                       return kmer_of(v);
                   });
    std::cout << kmer_of(path.back());
}

// non-member methods

namespace debug
{
    void print_edges(const Assembler& assembler)
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
    
    void print_vertices(const Assembler& assembler)
    {
        
    }
}
