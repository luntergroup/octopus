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
#include <tuple>
#include <cassert>
#include <numeric>
#include <limits>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/dominator_tree.hpp>

#include <iostream>

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
reference_kmers_ {}
{}

Assembler::Assembler(const unsigned kmer_size,
                     const SequenceType& reference)
:
k_ {kmer_size},
reference_kmers_ {}
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
                add_edge(u, *v, 0);
            }
        } else {
            const auto u = vertex_cache_.at(reference_kmers_.crbegin()[1]);
            const auto v = vertex_cache_.at(reference_kmers_.back());
            add_edge(u, v, 0);
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
    
    bool prev_kmer_good {true};
    
    if (!contains_kmer(reference_kmers_.back())) {
        const auto u = add_vertex(reference_kmers_.back(), true);
        if (!u) {
            prev_kmer_good = false;
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
                if (prev_kmer_good) {
                    const auto u = vertex_cache_.at(reference_kmers_.crbegin()[1]);
                    add_edge(u, *v, 1);
                }
                prev_kmer_good = true;
            } else {
                prev_kmer_good = false;
            }
        } else if (prev_kmer_good) {
            set_vertex_reference(reference_kmers_.back());
            const auto u = vertex_cache_.at(reference_kmers_.crbegin()[1]);
            const auto v = vertex_cache_.at(reference_kmers_.back());
            if (!boost::edge(u, v, graph_).second) {
                add_edge(u, v, 0);
            }
        } else {
            prev_kmer_good = true;
        }
    }
    
    vertex_cache_.rehash(vertex_cache_.size());
    reference_kmers_.shrink_to_fit();
    
    regenerate_vertex_indices();
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
class DfsAcyclicVisitor : public boost::default_dfs_visitor
{
public:
    using Edge = typename boost::graph_traits<G>::edge_descriptor;
    explicit DfsAcyclicVisitor(bool& is_acyclic) : is_acyclic_ {is_acyclic} {}
    void back_edge(Edge e, const G& g)
    {
        if (boost::source(e, g) != boost::target(e, g)) {
            is_acyclic_ = false;
        }
    }
private:
    bool& is_acyclic_;
};

template <typename G>
bool has_trivial_cycle(const G& graph)
{
    const auto p = boost::edges(graph);
    return std::any_of(p.first, p.second,
                       [&graph] (const auto& e) {
                           return boost::source(e, graph) == boost::target(e, graph);
                       });
}

bool Assembler::is_acyclic() const
{
    if (has_trivial_cycle(graph_)) {
        return false;
    }
    
    bool is_acyclic {true};
    
    DfsAcyclicVisitor<KmerGraph> vis {is_acyclic};
    
    const auto index_map = boost::get(&GraphNode::index, graph_);
    
    boost::depth_first_search(graph_, boost::visitor(vis).vertex_index_map(index_map));
    
    return is_acyclic;
}

bool Assembler::is_all_reference() const
{
    const auto p = boost::vertices(graph_);
    return std::all_of(p.first, p.second,
                       [this] (const Vertex& v) {
                           return is_reference(v);
                       });
}

void Assembler::prune(const unsigned min_weight)
{
    while (!empty()) {
        const auto old_size = boost::num_vertices(graph_);
        
        remove_low_weight_edges(min_weight);
        remove_disconnected_vertices();
        prune_reference_head();
        prune_reference_tail();
        
        if (reference_head() == reference_tail()) {
            clear();
            return;
        }
        
        const auto new_size1 = boost::num_vertices(graph_);
        
        if (new_size1 != old_size) {
            regenerate_vertex_indices();
        }
        
        prune_disconnected_subgraphs();
        prune_dangling_paths();
        cleanup_reference();
        prune_reference_head();
        prune_reference_tail();
        
        const auto new_size2 = boost::num_vertices(graph_);
        
        if (new_size1 != new_size2) {
            regenerate_vertex_indices();
        } else {
            break;
        }
        
        if (new_size2 == old_size) {
            break;
        }
    }
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
    
    if (empty()) {
        return result;
    }
    
    set_all_edge_log_probabilities_from(reference_head());
    
    while (max > 0 && !is_all_reference()) {
        extract_highest_probability_bubbles(result);
        --max;
    }
    
    return result;
}

// private methods

Assembler::GraphEdge::GraphEdge(unsigned weight)
:
weight {weight}
{}

bool Assembler::contains_kmer(const Kmer& kmer) const noexcept
{
    return vertex_cache_.count(kmer) > 0;
}

std::size_t Assembler::count_kmer(const Kmer& kmer) const noexcept
{
    return vertex_cache_.count(kmer);
}

void Assembler::set_vertex_reference(const Kmer& kmer)
{
    graph_[vertex_cache_.at(kmer)].is_reference = true;
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

void Assembler::remove_vertex(Vertex v)
{
    const auto c = vertex_cache_.erase(graph_[v].kmer);
    assert(c == 1);
    boost::remove_vertex(v, graph_);
}

void Assembler::clear_and_remove_vertex(Vertex v)
{
    //std::cout << "clearing vertex " << graph_[v].kmer << std::endl;
    const auto c = vertex_cache_.erase(graph_[v].kmer);
    assert(c == 1);
    boost::clear_vertex(v, graph_);
    boost::remove_vertex(v, graph_);
}

void Assembler::add_edge(const Vertex u, const Vertex v, const unsigned weight)
{
    boost::add_edge(u, v, GraphEdge {weight}, graph_);
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
    return is_source_reference(e) && is_target_reference(e);
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
    const auto p = boost::adjacent_vertices(u, graph_);
    return *std::find_if(p.first, p.second,
                         [this] (const Vertex v) {
                             return is_reference(v);
                         });
}

Assembler::Vertex Assembler::prev_reference(Vertex v) const
{
    const auto p = boost::inv_adjacent_vertices(v, graph_);
    return *std::find_if(p.first, p.second,
                         [this] (const Vertex u) {
                             return is_reference(u);
                         });
}

Assembler::SequenceType Assembler::make_reference(Vertex from, const Vertex to) const
{
    SequenceType result {};
    result.reserve(10 * k_);
    
    result.insert(std::end(result), std::cbegin(graph_[from].kmer), std::cend(graph_[from].kmer));
    
    from = next_reference(from);
    
    while (from != to) {
        result.push_back(graph_[from].kmer.back());
        from = next_reference(from);
    }
    
    result.shrink_to_fit();
    
    return result;
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

void Assembler::cleanup_reference()
{
    const auto it = std::find_if_not(std::begin(reference_kmers_),
                                     std::end(reference_kmers_),
                                     [this] (const Kmer& kmer) {
                                         return vertex_cache_.count(kmer) == 0;
                                     });
    
    reference_kmers_.erase(std::begin(reference_kmers_), it);
    
    if (reference_kmers_.empty()) {
        return;
    }
    
    const auto it2 = std::find_if_not(std::make_reverse_iterator(std::end(reference_kmers_)),
                                      std::make_reverse_iterator(it),
                                      [this] (const Kmer& kmer) {
                                          return vertex_cache_.count(kmer) == 0;
                                      }).base();
    
    reference_kmers_.erase(it2, std::end(reference_kmers_));
}

void Assembler::prune_reference_head()
{
    auto u = reference_head();
    
    if (boost::in_degree(u, graph_) > 0) {
        return;
    }
    
    const auto last = reference_tail();
    
    while (u != last && boost::in_degree(u, graph_) == 0
           && boost::out_degree(u, graph_) == 1) {
        const auto v = *boost::adjacent_vertices(u, graph_).first;
        remove_edge(u, v);
        remove_vertex(u);
        reference_kmers_.pop_front();
        u = v;
    }
}

void Assembler::prune_reference_tail()
{
    auto v = reference_tail();
    
    if (boost::out_degree(v, graph_) > 0) {
        return;
    }
    
    const auto first = reference_head();
    
    while (v != first && boost::in_degree(v, graph_) == 1
           && boost::out_degree(v, graph_) == 0) {
        const auto u = *boost::inv_adjacent_vertices(v, graph_).first;
        remove_edge(u, v);
        remove_vertex(v);
        reference_kmers_.pop_back();
        v = u;
    }
}

std::unordered_set<Assembler::Vertex> Assembler::find_reachable_kmers(Vertex src) const
{
    const auto index_map = boost::get(&GraphNode::index, graph_);
    
    std::unordered_set<Vertex> result {};
    
    auto vis = boost::make_bfs_visitor(
                                       boost::write_property(
                                                             boost::typed_identity_property_map<Vertex>(),
                                                             std::inserter(result, std::begin(result)),
                                                             boost::on_discover_vertex()));
    
    boost::breadth_first_search(graph_, src, boost::visitor(vis).vertex_index_map(index_map));
    
    return result;
}

void Assembler::prune_disconnected_subgraphs()
{
    const auto reference_reachable = find_reachable_kmers(reference_head());
    
    VertexIterator vi, vi_end, vi_next;
    
    std::tie(vi, vi_end) = boost::vertices(graph_);
    for (vi_next = vi; vi != vi_end; vi = vi_next) {
        ++vi_next;
        if (reference_reachable.count(*vi) == 0) {
            clear_and_remove_vertex(*vi);
        }
    }
}

std::deque<Assembler::Vertex> Assembler::find_dangling_tails() const
{
    std::deque<Vertex> result {};
    
    for (auto p = boost::vertices(graph_); p.first != p.second; ++p.first) {
        const auto v = *p.first;
        if (boost::out_degree(v, graph_) == 0  && v != reference_tail()) {
            result.push_back(v);
        }
    }
    
    return result;
}

void Assembler::prune_backwards_until_bifurcation(Vertex v)
{
    auto u = *boost::inv_adjacent_vertices(v, graph_).first;
    
    while (boost::out_degree(u, graph_) == 1) {
        remove_edge(u, v);
        remove_vertex(v);
        v = u;
        u = *boost::inv_adjacent_vertices(v, graph_).first;
    }
}

void Assembler::prune_dangling_paths()
{
    const auto tails = find_dangling_tails();
    for (const auto v : tails) {
        prune_backwards_until_bifurcation(v);
    }
}

std::unordered_map<Assembler::Vertex, Assembler::Vertex>
Assembler::build_dominator_tree(const Vertex from) const
{
    std::unordered_map<Vertex, Vertex> dom_tree;
    dom_tree.reserve(boost::num_vertices(graph_));
    
    auto dom_tree_pred_map = boost::make_assoc_property_map(dom_tree);
    
    boost::lengauer_tarjan_dominator_tree(graph_, from, dom_tree_pred_map);
    
    return dom_tree;
}

std::deque<Assembler::Vertex> Assembler::find_nondominants(const Vertex from) const
{
    const auto dom_tree = build_dominator_tree(from);
    
    std::unordered_set<Vertex> dominators {};
    dominators.reserve(dom_tree.size());
    
    for (const auto& p : dom_tree) {
        dominators.emplace(p.second);
    }
    
    std::deque<Vertex> result {};
    
    for (const auto& p : dom_tree) {
        if (dominators.count(p.first) == 0) {
            result.push_back(p.first);
        }
    }
    
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
struct DfsProbabilityVisitor : public boost::default_dfs_visitor
{
    using Vertex = typename boost::graph_traits<G>::vertex_descriptor;
    
    PropertyMap map;
    
    explicit DfsProbabilityVisitor(PropertyMap m) : map {m} {}
    
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
    
    DfsProbabilityVisitor<KmerGraph, decltype(prob_map)> vis {prob_map};
    
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

void Assembler::extract_highest_probability_bubbles(std::deque<Variant>& result)
{
    const auto wm = boost::get(&GraphEdge::neg_log_probability, graph_);
    
    std::unordered_map<Vertex, Vertex> pmi {};
    pmi.reserve(num_kmers());
    
    auto pm = boost::make_assoc_property_map(pmi);
    
    const auto index_map = boost::get(&GraphNode::index, graph_);
    
    std::cout << "calculating shortest paths..." << std::endl;
    
    boost::dag_shortest_paths(graph_, reference_head(),
                              boost::weight_map(wm).predecessor_map(pm).vertex_index_map(index_map));
    
    std::cout << "done calculating shortest paths..." << std::endl;
    
    auto v1 = reference_tail();
    std::cout << graph_[v1].kmer << std::endl;
    auto v2 = pm[v1];
    
    while (v2 != reference_head() && is_reference(v2) && v1 != v2) {
        std::cout << graph_[v2].kmer << std::endl;
        v1 = v2;
        v2 = pm[v1];
    }
    
    std::cout << graph_[v2].kmer << std::endl;
    
    if (v2 == reference_head()) {
        const auto nondominant_ref = find_nondominant_reference(v2);
        
        std::cout << "current graph" << std::endl;
        std::cout << "nondominant ref" << std::endl;
        for (auto v : nondominant_ref) {
            const auto p = boost::in_edges(v, graph_);
            std::for_each(p.first, p.second, [this] (const Edge e) {
                graph_[e].neg_log_probability = 100;
            });
            const auto p2 = boost::out_edges(v, graph_);
            std::for_each(p2.first, p2.second, [this] (const Edge e) {
                graph_[e].neg_log_probability = 100;
            });
            std::cout << graph_[v].kmer << std::endl;
        }
        
        for (auto ep = boost::edges(graph_); ep.first != ep.second; ++ep.first) {
            const auto e = *ep.first;
            print(e);
            std::cout << " weight = " << graph_[e].weight << " "
            << "(";
            if (is_source_reference(e)) {
                std::cout << "ref";
            } else {
                std::cout << "alt";
            }
            std::cout << ",";
            if (is_target_reference(e)) {
                std::cout << "ref";
            } else {
                std::cout << "alt";
            }
            std::cout << ") " << wm[e] << '\n';
        }
        
        //exit(0);
        
        return;
    }
    
    auto u1 = v2;
    auto u2 = pm[u1];
    
    while (u2 != reference_head()) {
        SequenceType alt {};
        alt.reserve(k_ * 10);
        
        alt.insert(std::end(alt),
                   std::make_reverse_iterator(std::cend(graph_[u1].kmer)),
                   std::make_reverse_iterator(std::cbegin(graph_[u1].kmer)));
        
        while (!is_reference(u2)) {
            alt.push_back(graph_[u2].kmer.front());
            u1 = u2;
            u2 = pm[u1];
        }
        
        alt.shrink_to_fit();
        
        std::reverse(std::begin(alt), std::end(alt));
        
        auto ref = make_reference(next_reference(u2), v1);
        
        std::cout << "ref = " << ref << std::endl;
        std::cout << "alt = " << alt << std::endl;
        
        result.emplace_back(0, std::move(ref), std::move(alt));
        
        auto b = find_bifurcation(u1, v1).first;
        
        if (b == v1) {
            remove_edge(u2, u1);
            while (u1 != v1) {
                auto tmp = *boost::adjacent_vertices(u1, graph_).first;
                remove_edge(u1, tmp);
                if (boost::in_degree(tmp, graph_) > 0) {
                    break;
                }
                remove_vertex(u1);
                u1 = tmp;
            }
            regenerate_vertex_indices();
            set_out_edge_log_probabilities(u2);
        } else {
            std::cout << "TODO" << std::endl;
            return;
        }
        
        while (u2 != reference_head() && is_reference(u2)) {
            u1 = u2;
            u2 = pm[u1];
        }
        
        const auto old_size = boost::num_vertices(graph_);
        prune_reference_tail();
        const auto new_size = boost::num_vertices(graph_);
        if (new_size != old_size) {
            regenerate_vertex_indices();
        }
    }
}

void Assembler::print(Edge e) const
{
    std::cout << graph_[boost::source(e, graph_)].kmer
    << "->"
    << graph_[boost::target(e, graph_)].kmer;
}

// non-member methods

namespace debug
{
    void print_edges(const Assembler& assembler)
    {
        for (auto ep = boost::edges(assembler.graph_); ep.first != ep.second; ++ep.first) {
            const auto e = *ep.first;
            assembler.print(e);
            std::cout << " weight = " << assembler.graph_[e].weight << " "
            << "(";
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
