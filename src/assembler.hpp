//
//  assembler.hpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef assembler_hpp
#define assembler_hpp

#include <vector>
#include <deque>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cstddef>
#include <utility>

#include <boost/graph/adjacency_list.hpp>
#include <boost/optional.hpp>

class Assembler;

namespace boost
{
    template <typename G>
    static decltype(auto) get(vertex_index_t, G& g);
    template <typename G>
    static decltype(auto) get(vertex_index_t, const G& g);
}

namespace debug
{
    void print_edges(const Assembler& assembler);
    void print_vertices(const Assembler& assembler);
}

class Assembler
{
public:
    using SequenceType = std::string;
    
    struct Variant
    {
        Variant() = default;
        explicit Variant(std::size_t pos, SequenceType ref, SequenceType alt);
        std::size_t begin_pos;
        SequenceType ref, alt;
    };
    
    Assembler() = delete;
    
    explicit Assembler(unsigned kmer_size);
    explicit Assembler(unsigned kmer_size, const SequenceType& reference);
    
    ~Assembler() = default;
    
    Assembler(const Assembler&)            = default;
    Assembler& operator=(const Assembler&) = default;
    Assembler(Assembler&&)                 = default;
    Assembler& operator=(Assembler&&)      = default;
    
    void insert_reference(const SequenceType& sequence);
    void insert_read(const SequenceType& sequence);
    
    std::size_t num_kmers() const noexcept;
    bool empty() const noexcept;
    bool is_acyclic() const;
    bool is_all_reference() const;
    
    void prune(unsigned min_weight);
    void clear();
    
    std::deque<Variant> extract_variants(unsigned max = 100);
    
    friend void debug::print_edges(const Assembler& assembler);
    friend void debug::print_vertices(const Assembler& assembler);
    
private:
    using Kmer = SequenceType;
    
    struct GraphEdge
    {
        GraphEdge() = default;
        GraphEdge(unsigned weight);
        
        unsigned weight;
        double neg_log_probability;
    };
    
    struct GraphNode
    {
        GraphNode() = default;
        template <typename T>
        explicit GraphNode(std::size_t index, T&& kmer, bool is_reference = false)
        :
        index {index},
        kmer {std::forward<T>(kmer)},
        is_reference {is_reference}
        {}
        
        std::size_t index;
        Kmer kmer;
        bool is_reference;
    };
    
    using KmerGraph = boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
    GraphNode, GraphEdge>;
    
    using Vertex = boost::graph_traits<KmerGraph>::vertex_descriptor;
    using Edge   = boost::graph_traits<KmerGraph>::edge_descriptor;
    
    using VertexIterator = boost::graph_traits<KmerGraph>::vertex_iterator;
    using EdgeIterator   = boost::graph_traits<KmerGraph>::edge_iterator;
    
    unsigned k_;
    
    std::deque<Kmer> reference_kmers_;
    
    KmerGraph graph_;
    
    std::unordered_map<Kmer, Vertex> vertex_cache_;
    
    bool contains_kmer(const Kmer& kmer) const noexcept;
    std::size_t count_kmer(const Kmer& kmer) const noexcept;
    
    void set_vertex_reference(const Kmer& kmer);
    void regenerate_vertex_indices();
    
    boost::optional<Vertex> add_vertex(const Kmer& kmer, bool is_reference = false);
    void remove_vertex(Vertex v);
    void clear_and_remove_vertex(Vertex v);
    void add_edge(Vertex u, Vertex v, unsigned weight);
    void remove_edge(Vertex u, Vertex v);
    void remove_edge(Edge e);
    void increment_weight(Edge e);
    
    bool is_reference(Vertex v) const;
    bool is_source_reference(Edge e) const;
    bool is_target_reference(Edge e) const;
    bool is_reference(Edge e) const;
    Vertex reference_head() const;
    Vertex reference_tail() const;
    
    Vertex next_reference(Vertex u) const;
    Vertex prev_reference(Vertex v) const;
    SequenceType make_reference(Vertex from, Vertex to) const;
    
    void remove_low_weight_edges(unsigned min_weight);
    void remove_disconnected_vertices();
    void cleanup_reference();
    void prune_reference_head();
    void prune_reference_tail();
    std::unordered_set<Vertex> find_reachable_kmers(Vertex src) const;
    void prune_disconnected_subgraphs();
    std::deque<Vertex> find_dangling_tails() const;
    void prune_backwards_until_bifurcation(Vertex v);
    void prune_dangling_paths();
    
    std::pair<Vertex, unsigned> find_bifurcation(Vertex from, Vertex to) const;
    
    std::unordered_map<Vertex, Vertex> build_dominator_tree(Vertex from) const;
    std::deque<Vertex> find_nondominants(Vertex from) const;
    std::deque<Vertex> find_nondominant_reference(Vertex from) const;
    
    void set_out_edge_log_probabilities(Vertex v);
    void set_all_edge_log_probabilities_from(Vertex src);
    
    void set_nondominant_reference_paths_impossible(Vertex from);
    
    void extract_highest_probability_bubbles(std::deque<Variant>& result);
    
    void print(Edge e) const;
    
    friend struct boost::property_map<KmerGraph, boost::vertex_index_t>;
    template <typename G>
    friend decltype(auto) boost::get(boost::vertex_index_t, G&);
    template <typename G>
    friend decltype(auto) boost::get(boost::vertex_index_t, const G&);
};

// Hack to make some older boost algorithms work with bundled properties
namespace boost
{
    template <>
    struct property_map<Assembler::KmerGraph, vertex_index_t>
    {
        using type       = property_map<Assembler::KmerGraph, std::size_t Assembler::GraphNode::*>::type;
        using const_type = property_map<Assembler::KmerGraph, std::size_t Assembler::GraphNode::*>::const_type;
    };
    
    template <typename G>
    static decltype(auto) get(vertex_index_t, G& g)
    {
        return get(&Assembler::GraphNode::index, g);
    }
    template <typename G>
    static decltype(auto)
    get(vertex_index_t, G const & g)
    {
        return get(&Assembler::GraphNode::index, g);
    }
} // namespace boost

#endif /* assembler_hpp */
