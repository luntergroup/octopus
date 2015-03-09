//
//  kmer_graph.h
//  Octopus
//
//  Created by Daniel Cooke on 22/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_kmer_graph_h
#define Octopus_kmer_graph_h

#include <string>
#include <vector>
#include <unordered_set>
#include <functional> // std::function
#include <tuple>      // std::tie
#include <list>
#include <iterator>   // std::cbegin etc
#include <stdexcept>
#include <cstdint>
#include <unordered_map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>

#include "utils.h"

#include <iostream>

template <typename ColourType, typename StringStoragePolicy>
class KmerGraph : StringStoragePolicy
{
public:
    using SizeType = std::uint_fast32_t;
    using typename StringStoragePolicy::InputType;
    
    KmerGraph() = delete;
    explicit KmerGraph(unsigned k);
    ~KmerGraph() = default;
    
    KmerGraph(const KmerGraph&)            = default;
    KmerGraph& operator=(const KmerGraph&) = default;
    KmerGraph(KmerGraph&&)                 = default;
    KmerGraph& operator=(KmerGraph&&)      = default;
    
    void add_sequence(InputType a_sequence, SizeType the_index, ColourType the_colour);
    void set_colour_weight_map(std::function<int(ColourType)> f_colour_weight);
    
    std::vector<std::string> get_contigs(unsigned max_num_paths);
    
    bool is_acyclic() const;
    unsigned get_num_kmers() const noexcept;
    unsigned get_num_connections() const noexcept;
    void clear();
    
    void print_kmers() const;
    
private:
    using typename StringStoragePolicy::ReferenceType;
    using StringStoragePolicy::store;
    
    struct KmerEdge
    {
        ReferenceType the_kmer;
        std::unordered_map<ColourType, unsigned> the_colours;
        int weight;
        std::vector<SizeType> the_indices;
    };
    
    using Graph = boost::adjacency_list<
        boost::listS, boost::listS, boost::directedS, boost::no_property, KmerEdge
    >;
    using Vertex                 = typename boost::graph_traits<Graph>::vertex_descriptor;
    using Edge                   = typename boost::graph_traits<Graph>::edge_descriptor;
    using VertexIterator         = typename boost::graph_traits<Graph>::vertex_iterator;
    using EdgeIterator           = typename boost::graph_traits<Graph>::edge_iterator;
    using OutEdgeIterator        = typename boost::graph_traits<Graph>::out_edge_iterator;
    using VertexPair             = std::pair<VertexIterator, VertexIterator>;
    using EdgePair               = std::pair<EdgeIterator, EdgeIterator>;
    using OutEdgePair            = std::pair<OutEdgeIterator, OutEdgeIterator>;
    using VertexIndexMapImpl     = std::unordered_map<Vertex, unsigned>;
    using VertexIndexMap         = boost::associative_property_map<VertexIndexMapImpl>;
    
    using Bubble  = std::pair<std::vector<Edge>, std::vector<Edge>>;
    using Bubbles = std::vector<Bubble>;
    
    class DfsVisitor : public boost::default_dfs_visitor
    {
    public:
        DfsVisitor(bool& is_acyclic) : is_acyclic_ {is_acyclic} {}
        void back_edge(Edge e, const Graph& g);
    private:
        bool& is_acyclic_;
    };
    
    const unsigned k_;
    Graph the_graph_;
    std::unordered_map<ReferenceType, Vertex> kmer_vertex_map_;
    std::function<int(ColourType)> f_colour_weight_;
    VertexIndexMapImpl vertex_indices_impl_;
    VertexIndexMap vertex_indices_;
    
    void add_kmer(ReferenceType the_kmer, SizeType the_index, ColourType the_colour);
    std::pair<Vertex, Vertex> get_vertices(ReferenceType a_kmer);
    std::pair<Vertex, bool> get_vertex(ReferenceType a_kmer_prefix_or_suffix) const;
    Vertex add_vertex(ReferenceType a_k_minus_1_mer);
    void add_edge(Vertex source, Vertex target, ReferenceType the_kmer, SizeType the_index,
                  ColourType the_colour);
    
    bool is_in_graph(ReferenceType a_k_minus_1_mer) const;
    unsigned get_next_index() const;
    
    Bubbles find_bubbles(Vertex the_source, SizeType min_index, SizeType max_index,
                         ColourType to_follow) const;
    std::vector<std::string> get_all_euler_paths(Vertex the_source, unsigned max_num_paths);
    template <typename ForwardIterator>
    std::string convert_path_to_string(ForwardIterator begin, ForwardIterator end) const;
    
    ReferenceType get_prefix(ReferenceType a_kmer) const;
    ReferenceType get_suffix(ReferenceType a_kmer) const;
};

template <typename ColourType, typename T>
KmerGraph<ColourType, T>::KmerGraph(unsigned k)
:k_ {k},
the_graph_ {},
kmer_vertex_map_ {},
f_colour_weight_ {[] (ColourType c) { return 1; }},
vertex_indices_impl_ {},
vertex_indices_ {vertex_indices_impl_}
{}

template <typename ColourType, typename T>
void KmerGraph<ColourType, T>::set_colour_weight_map(std::function<int(ColourType)> f_colour_weight)
{
    f_colour_weight_ = f_colour_weight;
}

template <typename ColourType, typename T>
void KmerGraph<ColourType, T>::add_sequence(InputType a_sequence, SizeType the_index,
                                                ColourType the_colour)
{
    if (a_sequence.size() < k_) return;
    ReferenceType sequence_ref = store(a_sequence);
    unsigned num_kmers = static_cast<unsigned>(sequence_ref.size()) - (k_ - 1);
    for (unsigned i = 0; i < num_kmers; ++i) {
        add_kmer(sequence_ref.substr(i, k_), the_index, the_colour);
    }
}

template <typename C, typename T>
std::vector<std::string> KmerGraph<C, T>::get_contigs(unsigned max_num_paths)
{
    return get_all_euler_paths(boost::vertex(0, the_graph_), max_num_paths);
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
    DfsVisitor a_visitor {is_acyclic};
    boost::depth_first_search(the_graph_, visitor(a_visitor).vertex_index_map(vertex_indices_));
    return is_acyclic;
}

template <typename C, typename T>
unsigned KmerGraph<C, T>::get_num_kmers() const noexcept
{
    return static_cast<unsigned>(boost::num_edges(the_graph_));
}

template <typename C, typename T>
unsigned KmerGraph<C, T>::get_num_connections() const noexcept
{
    return static_cast<unsigned>(boost::num_vertices(the_graph_));
}

template <typename C, typename T>
void KmerGraph<C, T>::clear()
{
    the_graph_.clear();
    kmer_vertex_map_.clear();
    vertex_indices_impl_.clear();
}

template <typename ColourType, typename T>
void KmerGraph<ColourType, T>::add_kmer(ReferenceType the_kmer, SizeType the_index,
                                            ColourType the_colour)
{
    Vertex source, target;
    std::tie(source, target) = get_vertices(the_kmer);
    add_edge(source, target, the_kmer, the_index, the_colour);
}

template <typename C, typename T>
std::pair<typename KmerGraph<C, T>::Vertex, typename KmerGraph<C, T>::Vertex>
KmerGraph<C, T>::get_vertices(ReferenceType a_kmer)
{
    auto kmer_prefix = get_prefix(a_kmer);
    auto kmer_suffix = get_suffix(a_kmer);
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
KmerGraph<C, T>::get_vertex(ReferenceType a_kmer_prefix_or_suffix) const
{
    EdgePair ep;
    for (ep = boost::edges(the_graph_); ep.first != ep.second; ++ep.first) {
        if (is_prefix(a_kmer_prefix_or_suffix, the_graph_[*ep.first].the_kmer)) {
            return {boost::source(*ep.first, the_graph_), true};
        }
        if (is_suffix(a_kmer_prefix_or_suffix, the_graph_[*ep.first].the_kmer)) {
            return {boost::target(*ep.first, the_graph_), true};
        }
    }
    return {boost::source(*ep.first, the_graph_), false};
}

template <typename C, typename T>
typename KmerGraph<C, T>::Vertex KmerGraph<C, T>::add_vertex(ReferenceType a_k_minus_1_mer)
{
    auto the_new_vertex = boost::add_vertex(the_graph_);
    kmer_vertex_map_.emplace(a_k_minus_1_mer, the_new_vertex);
    boost::put(vertex_indices_, the_new_vertex, get_next_index());
    return the_new_vertex;
}

template <typename ColourType, typename T>
void KmerGraph<ColourType, T>::add_edge(Vertex source, Vertex target, ReferenceType the_kmer,
                                            SizeType the_index, ColourType the_colour)
{
    auto the_existing_edge = boost::edge(source, target, the_graph_);
    if (the_existing_edge.second) {
        ++the_graph_[the_existing_edge.first].the_colours[the_colour];
        the_graph_[the_existing_edge.first].weight += f_colour_weight_(the_colour);
        the_graph_[the_existing_edge.first].the_indices.emplace_back(the_index);
    } else {
        auto a_new_edge = boost::add_edge(source, target, the_graph_).first;
        the_graph_[a_new_edge].the_kmer = the_kmer;
        ++the_graph_[a_new_edge].the_colours[the_colour];
        the_graph_[a_new_edge].weight = f_colour_weight_(the_colour);
        the_graph_[a_new_edge].the_indices.emplace_back(the_index);
    }
}

template <typename C, typename T>
bool KmerGraph<C, T>::is_in_graph(ReferenceType a_k_minus_1_mer) const
{
    return kmer_vertex_map_.count(a_k_minus_1_mer) > 0;
}

template <typename C, typename T>
unsigned KmerGraph<C, T>::get_next_index() const
{
    return static_cast<unsigned>(boost::num_vertices(the_graph_)) - 1;
}

template <typename C, typename T>
std::vector<std::string> KmerGraph<C, T>::get_all_euler_paths(Vertex the_source,
                                                                  unsigned max_num_paths)
{
    std::vector<std::string> euler_paths {};
    
//    std::list<Vertex> an_euler_path {};
//    std::list<Edge> path {};
//    auto all_verticies = boost::vertices(the_graph_);
//    std::unordered_set<Vertex> unvisited_verticies {all_verticies.first, all_verticies.second};
//    Vertex current_vertex {the_source};
//    unsigned num_edges_to_vist {get_num_kmers()};
//    OutEdgeIterator out_edge_begin, out_edge_end;
//    
//    while (num_edges_to_vist > 0) {
//        std::tie(out_edge_begin, out_edge_end) = boost::out_edges(current_vertex, the_graph_);
//        for (; out_edge_begin != out_edge_end; ++out_edge_begin) {
//            if (edge_counts[*out_edge_begin] < the_graph_[*out_edge_begin].the_colours.size()) {
//                ++edge_counts[*out_edge_begin];
//                an_euler_path.emplace_back(current_vertex);
//                path.emplace_back(*out_edge_begin);
//                current_vertex = boost::target(*out_edge_begin, the_graph_);
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
KmerGraph<ColourType, T>::find_bubbles(Vertex the_source, SizeType min_index, SizeType max_index,
                                       ColourType to_follow) const
{
    Bubbles result {};
    
    
    
    return result;
}

template <typename C, typename T>
template <typename ForwardIterator>
std::string KmerGraph<C, T>::convert_path_to_string(ForwardIterator begin, ForwardIterator end) const
{
    auto the_first_kmer_prefix = get_prefix(the_graph_[*begin].the_kmer);
    std::string result {the_first_kmer_prefix.cbegin(), the_first_kmer_prefix.cend()};
    result.reserve(k_ + std::distance(begin, end));
    for (; begin != end; ++begin) {
        result.push_back(the_graph_[*begin].the_kmer.back());
    }
    return result;
}

template <typename C, typename T>
typename KmerGraph<C, T>::ReferenceType KmerGraph<C, T>::get_prefix(ReferenceType a_kmer) const
{
    return a_kmer.substr(0, k_ - 1);
}

template <typename C, typename T>
typename KmerGraph<C, T>::ReferenceType KmerGraph<C, T>::get_suffix(ReferenceType a_kmer) const
{
    return a_kmer.substr(1, k_ - 1);
}

template <typename C, typename T>
void KmerGraph<C, T>::print_kmers() const
{
    auto kmer_map = boost::get(&KmerEdge::the_kmer, the_graph_);
    for (auto ep = edges(the_graph_); ep.first != ep.second; ++ep.first) {
        Edge e = *ep.first;
        std::cout << kmer_map[e] << "(" << the_graph_[e].the_colours.size() << ") ";
    }
    std::cout << std::endl;
}

#endif
