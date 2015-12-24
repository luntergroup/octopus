//
//  kmer_graph.hpp
//  Octopus
//
//  Created by Daniel Cooke on 22/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_kmer_graph_hpp
#define Octopus_kmer_graph_hpp

#include <string>
#include <vector>
#include <unordered_set>
#include <functional>
#include <tuple>
#include <list>
#include <iterator>
#include <stdexcept>
#include <cstdint>
#include <unordered_map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/functional/hash.hpp>

#include "string_utils.hpp"
#include "hash_functions.hpp"

#include <iostream>

//namespace Octopus {

template <typename ColourType, typename StringStoragePolicy>
class KmerGraph : StringStoragePolicy
{
public:
    using typename StringStoragePolicy::InputType;
    
    using SizeType = std::uint_fast32_t;
    
    KmerGraph() = delete;
    explicit KmerGraph(unsigned k);
    //explicit KmerGraph(unsigned k, std::function<int(ColourType)> f_colour_weight);
    ~KmerGraph() = default;
    
    KmerGraph(const KmerGraph&)            = default;
    KmerGraph& operator=(const KmerGraph&) = default;
    KmerGraph(KmerGraph&&)                 = default;
    KmerGraph& operator=(KmerGraph&&)      = default;
    
    void add_sequence(InputType sequence, SizeType position, ColourType colour);
    
    std::vector<std::string> get_contigs(unsigned max_num_paths);
    
    bool is_acyclic() const;
    unsigned num_kmers() const noexcept;
    unsigned num_connections() const noexcept;
    void clear();
    
    void print_kmers() const;
    void print_kmers(SizeType position) const;
    void print_kmers(SizeType from, ColourType colour) const;
    
private:
    using typename StringStoragePolicy::ReferenceType;
    using StringStoragePolicy::store;
    
    struct KmerEdge
    {
        ReferenceType kmer;
        std::unordered_map<ColourType, unsigned> colours;
        int weight;
        std::vector<SizeType> positions;
    };
    
    using Graph = boost::adjacency_list<
        boost::listS, boost::listS, boost::directedS, boost::no_property, KmerEdge
    >;
    using Vertex             = typename boost::graph_traits<Graph>::vertex_descriptor;
    using Edge               = typename boost::graph_traits<Graph>::edge_descriptor;
    using VertexIterator     = typename boost::graph_traits<Graph>::vertex_iterator;
    using EdgeIterator       = typename boost::graph_traits<Graph>::edge_iterator;
    using OutEdgeIterator    = typename boost::graph_traits<Graph>::out_edge_iterator;
    using VertexPair         = std::pair<VertexIterator, VertexIterator>;
    using EdgePair           = std::pair<EdgeIterator, EdgeIterator>;
    using OutEdgePair        = std::pair<OutEdgeIterator, OutEdgeIterator>;
    using VertexIndexMapImpl = std::unordered_map<Vertex, unsigned>;
    using VertexIndexMap     = boost::associative_property_map<VertexIndexMapImpl>;
    
    using Bubble  = std::pair<std::vector<Edge>, std::vector<Edge>>;
    using Bubbles = std::vector<Bubble>;
    
    class DfsVisitor : public boost::default_dfs_visitor
    {
    public:
        explicit DfsVisitor(bool& is_acyclic) : is_acyclic_ {is_acyclic} {}
        void back_edge(Edge e, const Graph& g);
    private:
        bool& is_acyclic_;
    };
    
    struct EdgeHash
    {
        size_t operator()(const Edge& edge) const
        {
            size_t seed {};
            //boost::hash_combine(seed, boost::hash<Vertex>()(boost::source(edge, graph_)));
            //boost::hash_combine(seed, boost::hash<Vertex>()(boost::target(edge, graph_)));
            return seed;
        }
    };
    
    using PositionMap = std::unordered_map<SizeType, std::unordered_set<Edge, EdgeHash>>;
    
    const unsigned k_;
    Graph graph_;
    std::unordered_map<ReferenceType, Vertex> kmer_vertex_map_;
    std::function<int(ColourType)> f_colour_weight_;
    VertexIndexMapImpl vertex_indices_impl_;
    VertexIndexMap vertex_indices_;
    PositionMap positions_;
    
    void add_kmer(ReferenceType kmer, SizeType index, ColourType colour);
    std::pair<Vertex, Vertex> get_vertices(ReferenceType kmer);
    std::pair<Vertex, bool> get_vertex(ReferenceType kmer_prefix_or_suffix) const;
    Vertex add_vertex(ReferenceType k_minus_1_mer);
    void add_edge(Vertex source, Vertex target, ReferenceType kmer, SizeType index, ColourType colour);
    
    bool is_in_graph(ReferenceType k_minus_1_mer) const;
    unsigned get_next_index() const;
    
    Bubbles find_bubbles(Vertex source, SizeType min_index, SizeType max_index, ColourType to_follow) const;
    std::vector<std::string> get_all_euler_paths(Vertex source, unsigned max_num_paths);
    template <typename ForwardIterator>
    std::string convert_path_to_string(ForwardIterator begin, ForwardIterator end) const;
    
    ReferenceType get_prefix(ReferenceType kmer) const;
    ReferenceType get_suffix(ReferenceType kmer) const;
};

template <typename ColourType, typename T>
KmerGraph<ColourType, T>::KmerGraph(unsigned k)
:
k_ {k},
graph_ {},
kmer_vertex_map_ {},
f_colour_weight_ {[] (ColourType c) { return 1; }},
vertex_indices_impl_ {},
vertex_indices_ {vertex_indices_impl_}
{}

template <typename ColourType, typename T>
void KmerGraph<ColourType, T>::add_sequence(InputType sequence, SizeType position, ColourType colour)
{
    if (sequence.size() < k_) return;
    
    ReferenceType sequence_ref = store(sequence);
    
    auto num_kmers = static_cast<unsigned>(sequence_ref.size()) - (k_ - 1);
    
    for (unsigned i = 0; i < num_kmers; ++i, ++position) {
        add_kmer(sequence_ref.substr(i, k_), position, colour);
    }
}

template <typename C, typename T>
std::vector<std::string> KmerGraph<C, T>::get_contigs(unsigned max_num_paths)
{
    return get_all_euler_paths(boost::vertex(0, graph_), max_num_paths);
}

template <typename C, typename T>
void KmerGraph<C, T>::DfsVisitor::back_edge(Edge e, const Graph &g)
{
    if (boost::source(e, g) != boost::target(e, g)) is_acyclic_ = false;
}

template <typename C, typename T>
bool KmerGraph<C, T>::is_acyclic() const
{
    bool is_acyclic {true};
    DfsVisitor visitor {is_acyclic};
    boost::depth_first_search(graph_, visitor(visitor).vertex_index_map(vertex_indices_));
    return is_acyclic;
}

template <typename C, typename T>
unsigned KmerGraph<C, T>::num_kmers() const noexcept
{
    return static_cast<unsigned>(boost::num_edges(graph_));
}

template <typename C, typename T>
unsigned KmerGraph<C, T>::num_connections() const noexcept
{
    return static_cast<unsigned>(boost::num_vertices(graph_));
}

template <typename C, typename T>
void KmerGraph<C, T>::clear()
{
    graph_.clear();
    kmer_vertex_map_.clear();
    vertex_indices_impl_.clear();
}

template <typename ColourType, typename T>
void KmerGraph<ColourType, T>::add_kmer(ReferenceType kmer, SizeType index, ColourType colour)
{
    Vertex source, target;
    std::tie(source, target) = get_vertices(kmer);
    add_edge(source, target, kmer, index, colour);
}

template <typename C, typename T>
std::pair<typename KmerGraph<C, T>::Vertex, typename KmerGraph<C, T>::Vertex>
KmerGraph<C, T>::get_vertices(ReferenceType kmer)
{
    auto kmer_prefix = get_prefix(kmer);
    auto kmer_suffix = get_suffix(kmer);
    
    Vertex source, target;
    if (kmer_prefix == kmer_suffix) {
        source = (is_in_graph(kmer_prefix)) ? get_vertex(kmer_prefix).first : add_vertex(kmer_prefix);
        target = source;
    } else {
        source = (is_in_graph(kmer_prefix)) ? get_vertex(kmer_prefix).first : add_vertex(kmer_prefix);
        target = (is_in_graph(kmer_suffix)) ? get_vertex(kmer_suffix).first : add_vertex(kmer_suffix);
    }
    
    return {source, target};
}

template <typename C, typename T>
std::pair<typename KmerGraph<C, T>::Vertex, bool>
KmerGraph<C, T>::get_vertex(ReferenceType kmer_prefix_or_suffix) const
{
    EdgePair ep;
    
    for (ep = boost::edges(graph_); ep.first != ep.second; ++ep.first) {
        if (Octopus::is_prefix(kmer_prefix_or_suffix, graph_[*ep.first].kmer)) {
            return {boost::source(*ep.first, graph_), true};
        }
        if (Octopus::is_suffix(kmer_prefix_or_suffix, graph_[*ep.first].kmer)) {
            return {boost::target(*ep.first, graph_), true};
        }
    }
    
    return {boost::source(*ep.first, graph_), false};
}

template <typename C, typename T>
typename KmerGraph<C, T>::Vertex KmerGraph<C, T>::add_vertex(ReferenceType k_minus_1_mer)
{
    auto result = boost::add_vertex(graph_);
    kmer_vertex_map_.emplace(k_minus_1_mer, result);
    boost::put(vertex_indices_, result, get_next_index());
    return result;
}

template <typename ColourType, typename T>
void KmerGraph<ColourType, T>::add_edge(Vertex source, Vertex target, ReferenceType kmer,
                                        SizeType position, ColourType colour)
{
    auto existing_edge = boost::edge(source, target, graph_);
    
    if (existing_edge.second) {
        ++graph_[existing_edge.first].colours[colour];
        graph_[existing_edge.first].weight += f_colour_weight_(colour);
        graph_[existing_edge.first].positions.emplace_back(position);
        positions_[position].emplace(existing_edge.first);
    } else {
        auto new_edge = boost::add_edge(source, target, graph_).first;
        
        graph_[new_edge].kmer = kmer;
        ++graph_[new_edge].colours[colour];
        graph_[new_edge].weight = f_colour_weight_(colour);
        graph_[new_edge].positions.emplace_back(position);
        positions_[position].emplace(new_edge);
    }
}

template <typename C, typename T>
bool KmerGraph<C, T>::is_in_graph(ReferenceType k_minus_1_mer) const
{
    return kmer_vertex_map_.count(k_minus_1_mer) > 0;
}

template <typename C, typename T>
unsigned KmerGraph<C, T>::get_next_index() const
{
    return static_cast<unsigned>(boost::num_vertices(graph_)) - 1;
}

template <typename C, typename T>
std::vector<std::string> KmerGraph<C, T>::get_all_euler_paths(Vertex source, unsigned max_num_paths)
{
    std::vector<std::string> euler_paths {};
    
//    std::list<Vertex> an_euler_path {};
//    std::list<Edge> path {};
//    auto all_verticies = boost::vertices(graph_);
//    std::unordered_set<Vertex> unvisited_verticies {all_verticies.first, all_verticies.second};
//    Vertex current_vertex {the_source};
//    unsigned num_edges_to_vist {num_kmers()};
//    OutEdgeIterator out_edge_begin, out_edge_end;
//    
//    while (num_edges_to_vist > 0) {
//        std::tie(out_edge_begin, out_edge_end) = boost::out_edges(current_vertex, graph_);
//        for (; out_edge_begin != out_edge_end; ++out_edge_begin) {
//            if (edge_counts[*out_edge_begin] < graph_[*out_edge_begin].colours.size()) {
//                ++edge_counts[*out_edge_begin];
//                an_euler_path.emplace_back(current_vertex);
//                path.emplace_back(*out_edge_begin);
//                current_vertex = boost::target(*out_edge_begin, graph_);
//                --num_edges_to_vist;
//                break;
//            }
//        }
//        if (out_edge_begin == out_edge_end) {
//            std::cout << "There are " << num_edges_to_vist << " edges remaining" << std::endl;
//            break;
//        }
//    }
//    
//    euler_paths.emplace_back(convert_path_to_string(path.cbegin(), path.cend()));
    
    return euler_paths;
}

template <typename ColourType, typename T>
typename KmerGraph<ColourType, T>::Bubbles
KmerGraph<ColourType, T>::find_bubbles(Vertex source, SizeType min_index, SizeType max_index,
                                       ColourType to_follow) const
{
    Bubbles result {};
    
    return result;
}

template <typename C, typename T>
template <typename ForwardIterator>
std::string KmerGraph<C, T>::convert_path_to_string(ForwardIterator begin, ForwardIterator end) const
{
    auto first_kmer_prefix = get_prefix(graph_[*begin].kmer);
    
    std::string result {first_kmer_prefix.cbegin(), first_kmer_prefix.cend()};
    result.reserve(k_ + std::distance(begin, end));
    
    for (; begin != end; ++begin) {
        result.push_back(graph_[*begin].kmer.back());
    }
    
    return result;
}

template <typename C, typename T>
typename KmerGraph<C, T>::ReferenceType KmerGraph<C, T>::get_prefix(ReferenceType kmer) const
{
    return kmer.substr(0, k_ - 1);
}

template <typename C, typename T>
typename KmerGraph<C, T>::ReferenceType KmerGraph<C, T>::get_suffix(ReferenceType kmer) const
{
    return kmer.substr(1, k_ - 1);
}

template <typename C, typename T>
void KmerGraph<C, T>::print_kmers() const
{
    auto kmer_map = boost::get(&KmerEdge::kmer, graph_);
    for (auto ep = edges(graph_); ep.first != ep.second; ++ep.first) {
        Edge e = *ep.first;
        std::cout << kmer_map[e] << "(" << graph_[e].colours.size() << ") ";
    }
    std::cout << std::endl;
}

template <typename C, typename T>
void KmerGraph<C, T>::print_kmers(SizeType position) const
{
    auto kmer_map = boost::get(&KmerEdge::kmer, graph_);
    
    for (const auto& edge : positions_.at(position)) {
        std::cout << kmer_map[edge] << "(" << graph_[edge].colours.size() << ") ";
    }
}

template <typename ColourType, typename T>
void KmerGraph<ColourType, T>::print_kmers(SizeType from, ColourType colour) const
{
    auto kmer_map = boost::get(&KmerEdge::kmer, graph_);
    auto edges = positions_.at(from);
    
    auto edge = *std::find_if(std::cbegin(edges), std::cend(edges),
                              [this, colour] (const auto& edge) {
                                  return graph_[edge].colours.count(colour) > 0;
                              });
    
    std::cout <<  kmer_map[edge] << std::endl;
}

//} // namespace Octopus

#endif
