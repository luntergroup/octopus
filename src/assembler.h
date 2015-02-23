//
//  assembler.h
//  Octopus
//
//  Created by Daniel Cooke on 20/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__assembler__
#define __Octopus__assembler__

#include <string>
#include <vector>
#include <unordered_set>
#include <map>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>

#include "variant.h"
#include "genomic_region.h"
#include "aligned_read.h"

using std::cbegin;
using std::cend;
using std::uint_fast32_t;

using boost::adjacency_list;
using boost::graph_traits;
using boost::vertices;
using boost::edges;

class Assembler
{
public:
    Assembler() = delete;
    Assembler(unsigned k);
    ~Assembler() = default;
    
    Assembler(const Assembler&)            = default;
    Assembler& operator=(const Assembler&) = default;
    Assembler(Assembler&&)                 = default;
    Assembler& operator=(Assembler&&)      = default;
    
    void add_reference_contig(GenomicRegion the_region, const std::string& a_contig);
    void add_read(const AlignedRead& a_read);
    
    std::vector<std::string> get_contigs();
    std::vector<Variant> get_variants();
    void clear();
    
    // TEST: delete after testing
    unsigned get_num_verticies() const noexcept;
    unsigned get_num_edges() const noexcept;
    void print_vertices() const;
    void print_edges() const;
    
private:
    
    enum class Source
    {
        Reference, Read, ReferenceAndRead
    };
    
    struct Kmer
    {
        std::string the_kmer;
        Source the_source;
        GenomicRegion the_region;
        int weight;
    };
    
    using Graph_t = adjacency_list<
        boost::multisetS, boost::listS, boost::bidirectionalS, boost::no_property, Kmer
    >;
    using Vertex                 = graph_traits<Graph_t>::vertex_descriptor;
    using Edge                   = graph_traits<Graph_t>::edge_descriptor;
    using VertexIterator         = graph_traits<Graph_t>::vertex_iterator;
    using EdgeIterator           = graph_traits<Graph_t>::edge_iterator;
    using OutEdgeIterator        = graph_traits<Graph_t>::out_edge_iterator;
    using VertexPair             = std::pair<VertexIterator, VertexIterator>;
    using EdgePair               = std::pair<EdgeIterator, EdgeIterator>;
    using OutEdgePair            = std::pair<OutEdgeIterator, OutEdgeIterator>;
    using VertexIndexMap         = std::map<Vertex, unsigned>;
    using VertexIndexPropertyMap = boost::associative_property_map<VertexIndexMap>;
    
    const unsigned k_;
    Graph_t the_graph_;
    std::unordered_set<std::string> added_kmer_suffixes_and_prefixes_;
    
    void add_sequence(const std::string& the_sequence, Source the_source);
    void add_kmer(std::string the_kmer, Source the_source);
    std::pair<Vertex, Vertex> get_vertices(const std::string& a_kmer);
    std::pair<Assembler::Vertex, bool> get_vertex(const std::string& a_kmer_prefix_or_suffix) const;
    Vertex add_vertex(const std::string& a_k_minus_1_mer);
    bool merge_with_existing_edge(Vertex& source, Vertex& target, Source the_source);
    void add_edge(Vertex& source, Vertex& target, std::string&& the_kmer, Source the_source);
    bool is_from_another_source(Source lhs, Source rhs) const noexcept;
    bool is_in_graph(EdgeIterator an_edge) const;
    bool is_in_graph(const std::string& a_k_minus_1_mer) const;
    unsigned num_parallel_edges(const OutEdgePair& an_edge_range) const noexcept;
    bool is_single_edge(const OutEdgePair& an_edge_range) const noexcept;
    unsigned get_next_vertex_index() const;
    
    std::string get_prefix(const std::string& a_kmer) const;
    std::string get_suffix(const std::string& a_kmer) const;
    bool is_prefix(const std::string& lhs, const std::string& rhs) const noexcept;
    bool is_suffix(const std::string& lhs, const std::string& rhs) const noexcept;
    
    std::vector<std::string> get_all_euler_paths(Vertex the_source);
};

inline unsigned Assembler::get_num_verticies() const noexcept
{
    return static_cast<unsigned>(boost::num_vertices(the_graph_));
}

inline unsigned Assembler::get_num_edges() const noexcept
{
    return static_cast<unsigned>(boost::num_edges(the_graph_));
}

inline void Assembler::print_vertices() const
{
    for (auto vp = vertices(the_graph_); vp.first != vp.second; ++vp.first) {
        Vertex v = *vp.first;
        if (boost::out_degree(v, the_graph_) > 0) {
            Edge an_out_edge = *boost::out_edges(v, the_graph_).first;
            auto prefix = the_graph_[an_out_edge].the_kmer.substr(0, k_ - 1);
            std::cout << prefix << " ";
        } else {
            Edge an_in_edge = *boost::in_edges(v, the_graph_).first;
            auto suffix = the_graph_[an_in_edge].the_kmer.substr(1, k_ - 1);
            std::cout << suffix << " ";
        }
    }
    std::cout << std::endl;
}

inline void Assembler::print_edges() const
{
    auto kmer_map = boost::get(&Kmer::the_kmer, the_graph_);
    for (auto ep = edges(the_graph_); ep.first != ep.second; ++ep.first) {
        Edge e = *ep.first;
        std::cout << kmer_map[e] << " ";
    }
    std::cout << std::endl;
}

#endif /* defined(__Octopus__assembler__) */
