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

#include "timers.hpp"

namespace octopus { namespace coretools {

namespace {
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
} // namespace

// public methods

Assembler::BadReferenceSequence::BadReferenceSequence(NucleotideSequence reference_sequence)
:
std::invalid_argument {"bad reference sequence"},
reference_sequence_ {std::move(reference_sequence)}
{}

const char* Assembler::BadReferenceSequence::what() const noexcept
{
    return std::invalid_argument::what();
}

Assembler::Assembler(const unsigned kmer_size)
:
k_ {kmer_size},
reference_kmers_ {},
reference_head_position_ {0}
{}

Assembler::Assembler(const unsigned kmer_size,
                     const NucleotideSequence& reference)
:
k_ {kmer_size},
reference_kmers_ {},
reference_head_position_ {0}
{
    insert_reference_into_empty_graph(reference);
}

unsigned Assembler::kmer_size() const noexcept
{
    return k_;
}

void Assembler::insert_reference(const NucleotideSequence& sequence)
{
    if (is_empty()) {
        insert_reference_into_empty_graph(sequence);
    } else {
        insert_reference_into_populated_graph(sequence);
    }
}

void Assembler::insert_read(const NucleotideSequence& sequence)
{
    if (sequence.size() < k_) return;
    
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
    
    for (; it2 <= std::cend(sequence); ++it1, ++it2) {
        Kmer kmer {it1, it2};
        
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
            bool e_in_graph;
            std::tie(e, e_in_graph) = boost::edge(u, v, graph_);
            
            if (e_in_graph) {
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
    
    remove_vertices_that_cant_be_reached_from(reference_head());
    
    new_size = boost::num_vertices(graph_);
    
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) return true;
        old_size = new_size;
    }
    
    remove_vertices_past(reference_tail());
    
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) return true;
        old_size = new_size;
    }
    
    remove_vertices_that_cant_reach(reference_tail());
    
    new_size = boost::num_vertices(graph_);
    if (new_size != old_size) {
        regenerate_vertex_indices();
        if (new_size < 2) return true;
        old_size = new_size;
    }
    
    try {
        if (can_prune_reference_flanks()) {
            prune_reference_flanks();
        }
    } catch (boost::not_a_dag& e) {
        clear();
        return false;
    }
    
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
}

std::deque<Assembler::Variant> Assembler::extract_variants(const unsigned max)
{
    if (is_empty() || is_all_reference()) {
        return std::deque<Variant> {};
    }
    
    set_all_edge_transition_scores_from(reference_head());
    
    return extract_k_highest_scoring_bubble_paths(max);
}

// Assembler private types

Assembler::GraphEdge::GraphEdge(const GraphEdge::WeightType weight, const bool is_reference)
:
weight {weight},
is_reference {is_reference}
{}

// Kmer

Assembler::Kmer::Kmer(SequenceIterator first, SequenceIterator last) noexcept
:
first_ {first}, last_ {last},
hash_ {boost::hash_range(first_, last_)}
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

// Assembler private methods

void Assembler::insert_reference_into_empty_graph(const NucleotideSequence& sequence)
{
    if (sequence.size() < k_) {
        throw std::runtime_error {"Assembler:: reference length must >= kmer_size"};
    }
    
    vertex_cache_.reserve(sequence.size() + std::pow(4, 5));
    
    auto it1 = std::cbegin(sequence);
    auto it2 = std::next(it1, k_);
    
    reference_kmers_.emplace_back(it1, it2);
    
    if (!contains_kmer(reference_kmers_.back())) {
        const auto u = add_vertex(reference_kmers_.back(), true);
        if (!u) {
            throw BadReferenceSequence {sequence};
        }
    }
    
    ++it1;
    ++it2;
    
    for (; it2 <= std::cend(sequence); ++it1, ++it2) {
        reference_kmers_.emplace_back(it1, it2);
        
        const auto& kmer = reference_kmers_.back();
        
        if (!contains_kmer(kmer)) {
            const auto v = add_vertex(kmer, true);
            if (v) {
                const auto u = vertex_cache_.at(std::crbegin(reference_kmers_)[1]);
                add_reference_edge(u, *v);
            } else {
                throw BadReferenceSequence {sequence};
            }
        } else {
            const auto u = vertex_cache_.at(std::crbegin(reference_kmers_)[1]);
            const auto v = vertex_cache_.at(kmer);
            add_reference_edge(u, v);
        }
    }
    
    reference_kmers_.shrink_to_fit();
}

void Assembler::insert_reference_into_populated_graph(const NucleotideSequence& sequence)
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
            throw BadReferenceSequence {sequence};
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
                const auto u = vertex_cache_.at(std::crbegin(reference_kmers_)[1]);
                add_reference_edge(u, *v);
            } else {
                throw BadReferenceSequence {sequence};
            }
        } else {
            const auto u = vertex_cache_.at(std::crbegin(reference_kmers_)[1]);
            const auto v = vertex_cache_.at(reference_kmers_.back());
            
            set_vertex_reference(v);
            
            Edge e;
            bool e_in_graph;
            std::tie(e, e_in_graph) = boost::edge(u, v, graph_);
            
            if (e_in_graph) {
                set_edge_reference(e);
            } else {
                add_reference_edge(u, v);
            }
        }
    }
    
    vertex_cache_.rehash(vertex_cache_.size());
    reference_kmers_.shrink_to_fit();
    
    regenerate_vertex_indices();
    
    reference_head_position_ = 0;
}

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

namespace
{
    template <typename T>
    bool is_dna(const T& sequence)
    {
        return std::all_of(std::cbegin(sequence), std::cend(sequence),
                           [] (const char base) {
                               return base == 'A' || base == 'C' || base == 'G' || base == 'T';
                           });
    }
} // namespace

Assembler::Vertex Assembler::null_vertex() const
{
    return boost::graph_traits<KmerGraph>::null_vertex();
}

boost::optional<Assembler::Vertex> Assembler::add_vertex(const Kmer& kmer, const bool is_reference)
{
    if (!is_dna(kmer)) return boost::none;
    const auto u = boost::add_vertex(GraphNode {boost::num_vertices(graph_), kmer, is_reference},
                                     graph_);
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
    const auto c = vertex_cache_.erase(kmer_of(v));
    assert(c == 1);
    boost::clear_vertex(v, graph_);
    boost::remove_vertex(v, graph_);
}

void Assembler::clear_and_remove_all(const std::unordered_set<Vertex>& vertices)
{
    for (const Vertex v : vertices) {
        clear_and_remove_vertex(v);
    }
}

void Assembler::add_edge(const Vertex u, const Vertex v,
                         const GraphEdge::WeightType weight,
                         const bool is_reference)
{
    boost::add_edge(u, v, GraphEdge {weight, is_reference}, graph_);
}

void Assembler::add_reference_edge(const Vertex u, const Vertex v)
{
    add_edge(u, v, 0, true);
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
    return vertex_cache_.at(reference_kmers_.front());
}

Assembler::Vertex Assembler::reference_tail() const
{
    return vertex_cache_.at(reference_kmers_.back());
}

Assembler::Vertex Assembler::next_reference(const Vertex u) const
{
    const auto p = boost::out_edges(u, graph_);
    const auto it = std::find_if(p.first, p.second, [this] (const Edge e) { return is_reference(e); });
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

std::size_t Assembler::num_reference_kmers() const
{
    const auto p = boost::vertices(graph_);
    return std::count_if(p.first, p.second, [this] (const Vertex& v) { return is_reference(v); });
}

Assembler::NucleotideSequence Assembler::make_sequence(const Path& path) const
{
    NucleotideSequence result {};
    result.reserve(k_ + path.size() - 1);
    
    const auto& first_kmer = kmer_of(path.front());
    result.insert(std::end(result), std::cbegin(first_kmer), std::cend(first_kmer));
    
    std::for_each(std::next(std::cbegin(path)), std::cend(path),
                  [this, &result] (const Vertex v) {
                      result.push_back(back_base_of(v));
                  });
    
    return result;
}

Assembler::NucleotideSequence Assembler::make_reference(Vertex from, const Vertex to) const
{
    NucleotideSequence result {};
    
    const auto null = null_vertex();
    
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
    result.insert(std::end(result), std::cbegin(first_kmer), std::cend(first_kmer));
    
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
    return std::any_of(p.first, p.second,
                       [this] (const Edge& e) {
                           return is_trivial_cycle(e);
                       });
}

bool Assembler::is_simple_deletion(Edge e) const
{
    return !is_reference(e) && is_source_reference(e) && is_target_reference(e);
}

bool Assembler::is_on_path(const Edge e, const Path& path) const
{
    if (path.size() < 2) return false;
    
    auto it1 = std::cbegin(path);
    auto it2 = std::next(it1);
    
    const auto last = std::cend(path);
    
    Edge path_edge;
    bool good;
    
    for (; it2 != last; ++it1, ++it2) {
        std::tie(path_edge, good) = boost::edge(*it1, *it2, graph_);
        assert(good);
        if (path_edge == e) return true;
    }
    
    return false;
}

void Assembler::remove_trivial_nonreference_cycles()
{
    boost::remove_edge_if([this] (const Edge e) {
        return !is_reference(e) && is_trivial_cycle(e);
    }, graph_);
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
    
    for (const Vertex u : reachables) {
        clear_and_remove_vertex(u);
    }
}

bool Assembler::can_prune_reference_flanks() const
{
    return boost::out_degree(reference_head(), graph_) == 1 || boost::in_degree(reference_tail(), graph_) == 1;
}

void Assembler::prune_reference_flanks()
{
    if (is_reference_empty()) return;
    
    // NB: I don't think this topological_sort is really needed (just iterate from reference_head
    // and reference_tail). Leaving it in for now as it's helping to uncover bugs!
    
    std::deque<Vertex> sorted_vertices {};
    
    boost::topological_sort(graph_, std::front_inserter(sorted_vertices),
                            boost::vertex_index_map(boost::get(&GraphNode::index, graph_)));
    
    assert(sorted_vertices.front() == reference_head() && sorted_vertices.back() == reference_tail());
    
    const auto it = std::find_if_not(std::cbegin(sorted_vertices), std::cend(sorted_vertices),
                                     [this] (const Vertex v) {
                                         return boost::out_degree(v, graph_) == 1
                                         && is_reference(*boost::out_edges(v, graph_).first);
                                     });
    
    std::for_each(std::cbegin(sorted_vertices), it,
                  [this] (const Vertex u) {
                      remove_edge(u, *boost::adjacent_vertices(u, graph_).first);
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
                      remove_edge(*boost::inv_adjacent_vertices(u, graph_).first, u);
                      remove_vertex(u);
                      reference_kmers_.pop_back();
                  });
}

Assembler::DominatorMap
Assembler::build_dominator_tree(const Vertex from) const
{
    DominatorMap result;
    result.reserve(boost::num_vertices(graph_));
    
    boost::lengauer_tarjan_dominator_tree(graph_, from,  boost::make_assoc_property_map(result));
    
    auto it = std::begin(result);
    
    for (; it != std::end(result);) {
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
auto compute_transition_score(const T edge_weight, const T total_out_weight, const R max_score = 100)
{
    if (total_out_weight == 0) {
        return R {0};
    }
    
    if (edge_weight == 0) {
        return max_score;
    }
    
    return std::abs(std::log(static_cast<R>(edge_weight) / total_out_weight));
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

bool Assembler::is_blocked(Edge e) const
{
    return graph_[e].transition_score >= BlockedScore;
}

void Assembler::block_edge(const Edge e)
{
    graph_[e].transition_score = BlockedScore;
}

void Assembler::block_all_in_edges(const Vertex v)
{
    set_all_in_edge_transition_scores(v, BlockedScore);
}

bool Assembler::all_in_edges_are_blocked(const Vertex v) const
{
    const auto p = boost::in_edges(v, graph_);
    return std::all_of(p.first, p.second, [this] (const Edge e) {
        return graph_[e].transition_score == BlockedScore;
    });
}

void Assembler::block_all_vertices(const std::deque<Vertex>& vertices)
{
    for (const Vertex& v : vertices) {
        block_all_in_edges(v);
    }
}

bool Assembler::all_vertices_are_blocked(const std::deque<Vertex>& vertices) const
{
    return std::all_of(std::cbegin(vertices), std::cend(vertices),
                       [this] (const Vertex v) { return all_in_edges_are_blocked(v); });
}

Assembler::PredecessorMap Assembler::find_shortest_scoring_paths(const Vertex from) const
{
    assert(from != null_vertex());
    
    std::unordered_map<Vertex, Vertex> preds {};
    preds.reserve(boost::num_vertices(graph_));
    
    boost::dag_shortest_paths(graph_, from,
                              boost::weight_map(boost::get(&GraphEdge::transition_score, graph_))
                              .predecessor_map(boost::make_assoc_property_map(preds))
                              .vertex_index_map(boost::get(&GraphNode::index, graph_)));
    
    return preds;
}

bool Assembler::is_on_path(const Vertex v, const PredecessorMap& predecessors, const Vertex from) const
{
    if (v == from) return true;
    
    assert(predecessors.count(from) == 1);
    
    auto it = predecessors.find(from);
    
    while (it != std::end(predecessors) && it->first != it->second) {
        if (it->second == v) return true;
        it = predecessors.find(it->second);
    }
    
    return false;
}

bool Assembler::is_on_path(const Edge e, const PredecessorMap& predecessors, const Vertex from) const
{
    assert(predecessors.count(from) == 1);
    
    auto it1 = predecessors.find(from);
    auto it2 = predecessors.find(it1->second);
    
    const auto last = std::cend(predecessors);
    
    Edge path_edge;
    bool good;
    
    while (it2 != last && it1 != it2) {
        std::tie(path_edge, good) = boost::edge(it2->second, it1->second, graph_);
        assert(good);
        if (path_edge == e) {
            return true;
        }
        it1 = it2;
        it2 = predecessors.find(it1->second);
    }
    
    return false;
}

Assembler::Path Assembler::extract_full_path(const PredecessorMap& predecessors, const Vertex from) const
{
    assert(predecessors.count(from) == 1);
    
    Path result {from};
    
    auto it = predecessors.find(from);
    
    while (it != std::end(predecessors) && it->first != it->second) {
        result.push_front(it->second);
        it = predecessors.find(it->second);
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

Assembler::Path
Assembler::extract_nonreference_path(const PredecessorMap& predecessors, Vertex from) const
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
    // reverse becasue more likely to be a closer vertex
    return std::find(rfirst, rlast, dominator) != rlast;
}

std::deque<Assembler::Variant> Assembler::extract_k_highest_scoring_bubble_paths(unsigned k)
{
    // TODO: should implement Eppstein's algorithm here
    
    auto dominator_tree = build_dominator_tree(reference_head());
    
    auto num_remaining_alt_kmers = num_kmers() - num_reference_kmers();
    
    Vertex ref, alt;
    unsigned rhs_kmer_count;
    
    boost::optional<Edge> blocked_edge {};
    
    std::deque<Variant> result {};
    
    unsigned max_blockings {50}; // HACK
    
    while (k > 0 && num_remaining_alt_kmers > 0) {
        auto predecessors = find_shortest_scoring_paths(reference_head());
        
        if (blocked_edge) {
            if (max_blockings == 0) { // HACK
                return result; // HACK
            } // HACK
            
            --max_blockings; // HACK
            
            // TODO: This is almost certainly not optimal and is it even guaranteed to terminate?
            if (!is_on_path(boost::target(*blocked_edge, graph_), predecessors, reference_tail())) {
                set_out_edge_transition_scores(boost::source(*blocked_edge, graph_));
                blocked_edge = boost::none;
            } else {
                const auto p = boost::out_edges(boost::target(*blocked_edge, graph_), graph_);
                if (std::all_of(p.first, p.second, [this] (const Edge e) { return is_blocked(e); })) {
                    return result; // Othersie might not terminate?
                }
            }
        }
        
        assert(count_unreachables(predecessors) == 1);
        
        std::tie(alt, ref, rhs_kmer_count) = backtrack_until_nonreference(predecessors, reference_tail());
        
        if (alt == reference_head()) {
            // complete reference path is shortest path
            const auto nondominant_reference = extract_nondominant_reference(dominator_tree);
            
            if (all_vertices_are_blocked(nondominant_reference)) {
                return result; // nothing more we can do?
            }
            
            block_all_vertices(nondominant_reference);
            
            continue;
        }
        
        while (alt != reference_head()) {
            auto alt_path = extract_nonreference_path(predecessors, alt);
            
            assert(!alt_path.empty());
            
            const auto ref_before_bubble = predecessors.at(alt_path.front());
            
            auto ref_seq = make_reference(ref_before_bubble, ref);
            alt_path.push_front(ref_before_bubble);
            auto alt_seq = make_sequence(alt_path);
            alt_path.pop_front();
            
            rhs_kmer_count += count_kmers(ref_seq, k_);
            
            const auto pos = reference_head_position_ + reference_size() - sequence_length(rhs_kmer_count, k_);
            
            result.emplace_front(pos, std::move(ref_seq), std::move(alt_seq));
            
            --rhs_kmer_count; // because we padded one reference kmer to make ref_seq
            
            Edge edge_to_alt;
            bool good;
            std::tie(edge_to_alt, good) = boost::edge(alt, ref, graph_);
            
            assert(good);
            
            if (alt_path.size() == 1 && is_simple_deletion(edge_to_alt)) {
                remove_edge(alt_path.front(), ref);
                set_out_edge_transition_scores(alt_path.front());
            } else {
                auto vertex_before_bridge = ref_before_bubble;
                
                while (!alt_path.empty()) {
                    const auto it = is_bridge_until(alt_path);
                    
                    if (it == std::cend(alt_path)) {
                        remove_path(alt_path);
                        regenerate_vertex_indices();
                        set_out_edge_transition_scores(vertex_before_bridge);
                        erase_all(alt_path, dominator_tree);
                        num_remaining_alt_kmers -= alt_path.size();
                        alt_path.clear();
                    } else if (joins_reference_only(*it)) {
                        alt_path.erase(it, std::end(alt_path));
                        remove_path(alt_path);
                        regenerate_vertex_indices();
                        set_out_edge_transition_scores(vertex_before_bridge);
                        erase_all(alt_path, dominator_tree);
                        num_remaining_alt_kmers -= alt_path.size();
                        break;
                    } else if (is_dominated_by_path(*it, std::cbegin(alt_path), it, dominator_tree)) {
                        vertex_before_bridge = *it;
                        alt_path.erase(std::begin(alt_path), std::next(it));
                    } else {
                        // TODO: This is almost certainly not optimal... fortunatly it seems to
                        // be a rare case.
                        Edge e;
                        bool good;
                        std::tie(e, good) = boost::edge(*std::prev(it), *it, graph_);
                        assert(good);
                        block_edge(e);
                        blocked_edge = e;
                        break;
                    }
                }
            }
            
            unsigned kmer_count_to_alt;
            std::tie(alt, ref, kmer_count_to_alt) = backtrack_until_nonreference(predecessors, ref_before_bubble);
            
            rhs_kmer_count += kmer_count_to_alt;
            
            if (k > 0) --k;
        }
        
        assert(boost::out_degree(reference_head(), graph_) > 0);
        assert(boost::in_degree(reference_tail(), graph_) > 0);
        
        if (can_prune_reference_flanks()) {
            prune_reference_flanks();
            regenerate_vertex_indices();
            dominator_tree = build_dominator_tree(reference_head());
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

void Assembler::print_dominator_tree() const
{
    const auto dom_tree = build_dominator_tree(reference_head());
    
    for (const auto& p : dom_tree) {
        std::cout << kmer_of(p.first) << " dominated by " << kmer_of(p.second) << std::endl;
    }
}

// non-member methods

bool operator==(const Assembler::Variant& lhs, const Assembler::Variant& rhs)
{
    return lhs.begin_pos == rhs.begin_pos && lhs.alt == rhs.alt;
}

bool operator!=(const Assembler::Variant& lhs, const Assembler::Variant& rhs)
{
    return !operator==(lhs, rhs);
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
