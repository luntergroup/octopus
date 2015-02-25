//
//  kmer_assembler.h
//  Octopus
//
//  Created by Daniel Cooke on 22/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_kmer_assembler_h
#define Octopus_kmer_assembler_h

#include <string>
#include <vector>
#include <unordered_set>
#include <map>
#include <functional>
#include <tuple>
#include <list>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include "variant.h"

template <typename ColourType, typename StringStoragePolicy>
class KmerAssembler : public StringStoragePolicy
{
public:
    using typename StringStoragePolicy::InputType;
    
    KmerAssembler() = delete;
    KmerAssembler(unsigned k);
    ~KmerAssembler() = default;
    
    KmerAssembler(const KmerAssembler&)            = default;
    KmerAssembler& operator=(const KmerAssembler&) = default;
    KmerAssembler(KmerAssembler&&)                 = default;
    KmerAssembler& operator=(KmerAssembler&&)      = default;
    
    void add_sequence(InputType a_sequence, ColourType the_colour);
    
    std::vector<std::string> get_contigs();
    std::vector<Variant> get_variants(ColourType the_reference_colour);
    
    void set_colour_weight_map(std::function<int(ColourType)> f_colour_weight);
    void clear();
    
    unsigned get_num_kmers() const noexcept;
    
    void print_kmers() const;
    
private:
    using typename StringStoragePolicy::ReferenceType;
    using StringStoragePolicy::store;
    
    struct Kmer
    {
        ReferenceType the_kmer;
        std::vector<ColourType> the_colours;
        int weight;
    };
    
    using Graph_t = boost::adjacency_list<
        boost::listS, boost::listS, boost::bidirectionalS, boost::no_property, Kmer
    >;
    using Vertex                 = typename boost::graph_traits<Graph_t>::vertex_descriptor;
    using Edge                   = typename boost::graph_traits<Graph_t>::edge_descriptor;
    using VertexIterator         = typename boost::graph_traits<Graph_t>::vertex_iterator;
    using EdgeIterator           = typename boost::graph_traits<Graph_t>::edge_iterator;
    using OutEdgeIterator        = typename boost::graph_traits<Graph_t>::out_edge_iterator;
    using VertexPair             = std::pair<VertexIterator, VertexIterator>;
    using EdgePair               = std::pair<EdgeIterator, EdgeIterator>;
    using OutEdgePair            = std::pair<OutEdgeIterator, OutEdgeIterator>;
    using VertexIndexMap         = std::map<Vertex, unsigned>;
    using VertexIndexPropertyMap = boost::associative_property_map<VertexIndexMap>;
    using EdgeCountMap           = std::map<Edge, unsigned>;
    
    const unsigned k_;
    Graph_t the_graph_;
    std::unordered_set<ReferenceType> added_kmer_suffixes_and_prefixes_;
    std::function<int(ColourType)> f_colour_weight_;
    
    void add_kmer(ReferenceType the_kmer, ColourType the_colour);
    std::pair<Vertex, Vertex> get_vertices(ReferenceType a_kmer);
    std::pair<Vertex, bool> get_vertex(ReferenceType a_kmer_prefix_or_suffix) const;
    Vertex add_vertex(ReferenceType a_k_minus_1_mer);
    void add_edge(Vertex source, Vertex target, ReferenceType the_kmer, ColourType the_colour);
    bool is_in_graph(EdgeIterator an_edge) const;
    bool is_in_graph(ReferenceType a_k_minus_1_mer) const;
    template <typename C> std::string convert_path_to_string(const C& path) const;
    ReferenceType get_prefix(ReferenceType a_kmer) const;
    ReferenceType get_suffix(ReferenceType a_kmer) const;
    bool is_prefix(ReferenceType lhs, ReferenceType rhs) const noexcept;
    bool is_suffix(ReferenceType lhs, ReferenceType rhs) const noexcept;
    
    std::vector<std::string> get_all_euler_paths(Vertex the_ColourType);
};

template <typename ColourType, typename T>
KmerAssembler<ColourType, T>::KmerAssembler(unsigned k)
:k_ {k},
the_graph_ {},
added_kmer_suffixes_and_prefixes_ {},
f_colour_weight_ {[] (ColourType c) { return 1; }}
{}

template <typename ColourType, typename T>
void KmerAssembler<ColourType, T>::set_colour_weight_map(std::function<int(ColourType)> f_colour_weight)
{
    f_colour_weight_ = f_colour_weight;
}

template <typename ColourType, typename T>
void KmerAssembler<ColourType, T>::add_sequence(InputType a_sequence, ColourType the_colour)
{
    if (a_sequence.size() < k_) return;
    ReferenceType sequence_ref = store(a_sequence);
    unsigned num_kmers = static_cast<unsigned>(sequence_ref.size()) - (k_ - 1);
    for (unsigned i = 0; i < num_kmers; ++i) {
        add_kmer(sequence_ref.substr(i, k_), the_colour);
    }
}

template <typename ColourType, typename T>
void KmerAssembler<ColourType, T>::add_kmer(ReferenceType the_kmer, ColourType the_colour)
{
    Vertex source, target;
    std::tie(source, target) = get_vertices(the_kmer);
    add_edge(source, target, the_kmer, the_colour);
}

template <typename C, typename T>
std::pair<typename KmerAssembler<C, T>::Vertex, typename KmerAssembler<C, T>::Vertex>
KmerAssembler<C, T>::get_vertices(ReferenceType a_kmer)
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
std::pair<typename KmerAssembler<C, T>::Vertex, bool>
KmerAssembler<C, T>::get_vertex(ReferenceType a_kmer_prefix_or_suffix) const
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
typename KmerAssembler<C, T>::Vertex KmerAssembler<C, T>::add_vertex(ReferenceType a_k_minus_1_mer)
{
    added_kmer_suffixes_and_prefixes_.emplace(a_k_minus_1_mer);
    return boost::add_vertex(the_graph_);
}

template <typename ColourType, typename T>
void KmerAssembler<ColourType, T>::add_edge(Vertex source, Vertex target, ReferenceType the_kmer,
                                            ColourType the_colour)
{
    auto the_existing_edge = boost::edge(source, target, the_graph_);
    if (the_existing_edge.second) {
        the_graph_[the_existing_edge.first].the_colours.emplace_back(the_colour);
        the_graph_[the_existing_edge.first].weight += f_colour_weight_(the_colour);
    } else {
        auto a_new_edge = boost::add_edge(source, target, the_graph_).first;
        the_graph_[a_new_edge].the_kmer = the_kmer;
        the_graph_[a_new_edge].the_colours.emplace_back(the_colour);
        the_graph_[a_new_edge].weight = f_colour_weight_(the_colour);
    }
}

template <typename C, typename T>
bool KmerAssembler<C, T>::is_in_graph(EdgeIterator an_edge) const
{
    return an_edge != boost::edges(the_graph_).second;
}

template <typename C, typename T>
bool KmerAssembler<C, T>::is_in_graph(ReferenceType a_k_minus_1_mer) const
{
    return added_kmer_suffixes_and_prefixes_.count(a_k_minus_1_mer) > 0;
}

template <typename C, typename T>
std::vector<std::string> KmerAssembler<C, T>::get_contigs()
{
    return get_all_euler_paths(boost::vertex(0, the_graph_));
}

template <typename C, typename T>
std::vector<std::string> KmerAssembler<C, T>::get_all_euler_paths(Vertex the_source)
{
    std::vector<std::string> euler_paths {};
    
    std::list<Vertex> an_euler_path {};
    std::list<Edge> path {};
    auto all_verticies = boost::vertices(the_graph_);
    std::unordered_set<Vertex> unvisited_verticies {all_verticies.first, all_verticies.second};
    Vertex current_vertex {the_source};
    unsigned num_edges_to_vist {get_num_kmers()};
    EdgeCountMap edge_counts {};
    OutEdgeIterator out_edge_begin, out_edge_end;
    
    while (num_edges_to_vist > 0) {
        std::tie(out_edge_begin, out_edge_end) = boost::out_edges(current_vertex, the_graph_);
        for (; out_edge_begin != out_edge_end; ++out_edge_begin) {
            if (edge_counts[*out_edge_begin] < the_graph_[*out_edge_begin].the_colours.size()) {
                ++edge_counts[*out_edge_begin];
                an_euler_path.emplace_back(current_vertex);
                path.emplace_back(*out_edge_begin);
                current_vertex = boost::target(*out_edge_begin, the_graph_);
                --num_edges_to_vist;
                break;
            }
        }
        if (out_edge_begin == out_edge_end) {
            std::cout << "There are " << num_edges_to_vist << " edges remaining" << std::endl;
            break;
        }
    }
    
    euler_paths.emplace_back(convert_path_to_string(path));
    
    return euler_paths;
}

template <typename ColourType, typename T>
std::vector<Variant> KmerAssembler<ColourType, T>::get_variants(ColourType the_reference_colour)
{
    std::vector<Variant> result {};
    
    return result;
}

template <typename C, typename T>
template <typename Container>
std::string KmerAssembler<C, T>::convert_path_to_string(const Container& path) const
{
    std::string result {};
    result.reserve(k_ + path.size() - 1);
    for (const auto& edge : path) {
        auto the_edge_kmer = the_graph_[edge].the_kmer;
        if (result.empty()) {
            result += std::string {the_edge_kmer.cbegin(), the_edge_kmer.cend()};
        } else {
            result.push_back(the_edge_kmer.back());
        }
    }
    return result;
}

template <typename C, typename T>
void KmerAssembler<C, T>::clear()
{
    the_graph_.clear();
    added_kmer_suffixes_and_prefixes_.clear();
}

template <typename C, typename T>
unsigned KmerAssembler<C, T>::get_num_kmers() const noexcept
{
    return static_cast<unsigned>(boost::num_edges(the_graph_));
}

template <typename C, typename T>
typename KmerAssembler<C, T>::ReferenceType KmerAssembler<C, T>::get_prefix(ReferenceType a_kmer) const
{
    return a_kmer.substr(0, k_ - 1);
}

template <typename C, typename T>
typename KmerAssembler<C, T>::ReferenceType KmerAssembler<C, T>::get_suffix(ReferenceType a_kmer) const
{
    return a_kmer.substr(1, k_ - 1);
}

template <typename C, typename T>
bool KmerAssembler<C, T>::is_prefix(ReferenceType lhs, ReferenceType rhs) const noexcept
{
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}

template <typename C, typename T>
bool KmerAssembler<C, T>::is_suffix(ReferenceType lhs, ReferenceType rhs) const noexcept
{
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::next(std::cbegin(rhs)));
}

template <typename C, typename T>
void KmerAssembler<C, T>::print_kmers() const
{
    auto kmer_map = boost::get(&Kmer::the_kmer, the_graph_);
    for (auto ep = edges(the_graph_); ep.first != ep.second; ++ep.first) {
        Edge e = *ep.first;
        std::cout << kmer_map[e] << "(" << the_graph_[e].the_colours.size() << ") ";
    }
    std::cout << std::endl;
}

#endif
