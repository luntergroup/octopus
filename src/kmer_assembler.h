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
#include <boost/graph/adjacency_list.hpp>

template <typename T, typename C>
class KmerAssembler
{
    KmerAssembler() = delete;
    KmerAssembler(unsigned k);
    KmerAssembler() = default;
    
    KmerAssembler(const KmerAssembler&)            = default;
    KmerAssembler& operator=(const KmerAssembler&) = default;
    KmerAssembler(KmerAssembler&&)                 = default;
    KmerAssembler& operator=(KmerAssembler&&)      = default;
    
    template <typename T_>
    void add_sequence(T_&& the_sequence, C the_colour);
    void clear();

private:
    struct Kmer
    {
        T the_kmer;
        C the_colour;
        int weight;
    };
    
    using Graph_t = adjacency_list<
    boost::multisetS, boost::listS, boost::bidirectionalS, boost::no_property, Kmer
    >;
    using Vertex            = graph_traits<Graph_t>::vertex_descriptor;
    using Edge              = graph_traits<Graph_t>::edge_descriptor;
    using VertexIterator    = graph_traits<Graph_t>::vertex_iterator;
    using EdgeIterator      = graph_traits<Graph_t>::edge_iterator;
    using OutEdgeIterator   = graph_traits<Graph_t>::out_edge_iterator;
    using VertexPair        = std::pair<VertexIterator, VertexIterator>;
    using EdgePair          = std::pair<EdgeIterator, EdgeIterator>;
    using OutEdgePair       = std::pair<OutEdgeIterator, OutEdgeIterator>;
    
    const unsigned k_;
    Graph_t the_graph_;
    std::unordered_set<T> added_kmer_suffixes_and_prefixes_;
    
    template <typename T_>
    void add_sequence(T_&& the_sequence, C the_colour);
    template <typename T_>
    void add_kmer(T_&& the_kmer, C the_colour);
    std::pair<Vertex, Vertex> get_vertices(const T& a_kmer);
    std::pair<Assembler::Vertex, bool> get_vertex(const T& a_kmer_prefix_or_suffix) const;
    Vertex add_vertex(const T& a_k_minus_1_mer);
    bool merge_with_existing_edge(Vertex& source, Vertex& target, C the_colour);
    template <typename T_>
    void add_edge(Vertex& source, Vertex& target, T_&& the_kmer, C the_colour);
    bool is_from_another_source(Source lhs, Source rhs) const noexcept;
    bool is_in_graph(EdgeIterator an_edge) const;
    bool is_in_graph(const T& a_k_minus_1_mer) const;
    unsigned num_parallel_edges(const OutEdgePair& an_edge_range) const noexcept;
    bool is_single_edge(const OutEdgePair& an_edge_range) const noexcept;
    
    std::string get_prefix(const T& a_kmer) const;
    std::string get_suffix(const T& a_kmer) const;
    bool is_prefix(const T& lhs, const T& rhs) const noexcept;
    bool is_suffix(const T& lhs, const T& rhs) const noexcept;
};

template <typename T, typename C>
KmerAssembler::KmerAssembler(unsigned k)
:k_ {k},
the_graph_ {},
added_kmer_suffixes_and_prefixes_ {}
{}

template <typename T, typename C>
void KmerAssembler::add_reference_contig(const std::string& a_contig)
{
    add_sequence(a_contig, Source::Reference);
}

template <typename T, typename C>

template <typename T, typename C>

template <typename T, typename C>

template <typename T, typename C>

template <typename T, typename C>

template <typename T, typename C>

template <typename T, typename C>

template <typename T, typename C>

template <typename T, typename C>

#endif
