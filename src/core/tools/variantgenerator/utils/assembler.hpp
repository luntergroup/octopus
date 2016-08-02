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
#include <tuple>
#include <stdexcept>
#include <iosfwd>

#include <boost/graph/adjacency_list.hpp>
#include <boost/optional.hpp>

#include <interfaces/comparable.hpp>

namespace octopus { namespace coretools { class Assembler; }}

namespace boost {
    template <typename G>
    static decltype(auto) get(vertex_index_t, G& g);
    template <typename G>
    static decltype(auto) get(vertex_index_t, const G& g);
}

namespace octopus { namespace coretools {

namespace debug {
    void print(const Assembler& assembler);
}

class Assembler
{
public:
    using NucleotideSequence = std::string;
    
    struct Variant
    {
        Variant() = default;
        template <typename S1, typename S2>
        Variant(std::size_t pos, S1&& ref, S2&& alt);
        std::size_t begin_pos;
        NucleotideSequence ref, alt;
    };
    
    class BadReferenceSequence : public std::invalid_argument
    {
    public:
        BadReferenceSequence(NucleotideSequence reference_sequence);
        ~BadReferenceSequence() noexcept = default;
        const char* what() const noexcept override;
    private:
        NucleotideSequence reference_sequence_;
    };
    
    Assembler() = delete;
    
    Assembler(unsigned kmer_size);
    Assembler(unsigned kmer_size, const NucleotideSequence& reference);
    
    Assembler(const Assembler&)            = delete;
    Assembler& operator=(const Assembler&) = delete;
    Assembler(Assembler&&)                 = default;
    Assembler& operator=(Assembler&&)      = default;
    
    ~Assembler() = default;
    
    unsigned kmer_size() const noexcept;
    
    void insert_reference(const NucleotideSequence& sequence);
    void insert_read(const NucleotideSequence& sequence);
    
    std::size_t num_kmers() const noexcept;
    bool is_empty() const noexcept;
    
    bool is_acyclic() const;
    
    bool is_all_reference() const;
    
    bool prune(unsigned min_weight);
    void clear();
    
    std::deque<Variant> extract_variants(unsigned max = 100);
    
    friend void debug::print(const Assembler& assembler);
    
private:
    class Kmer : public Comparable<Kmer>
    {
    public:
        using NucleotideSequence     = Assembler::NucleotideSequence;
        using SequenceIterator = NucleotideSequence::const_iterator;
        
        Kmer() = delete;
        
        Kmer(SequenceIterator first, SequenceIterator last) noexcept;
        
        Kmer(const Kmer&)            = default;
        Kmer& operator=(const Kmer&) = default;
        Kmer(Kmer&&)                 = default;
        Kmer& operator=(Kmer&&)      = default;
        
        ~Kmer() = default;
        
        char front() const noexcept;
        char back() const noexcept;
        
        SequenceIterator begin() const noexcept;
        SequenceIterator end() const noexcept;
        
        explicit operator NucleotideSequence() const;
        
        std::size_t hash() const noexcept;
        
        friend bool operator==(const Kmer& lhs, const Kmer& rhs) noexcept;
        friend bool operator<(const Kmer& lhs, const Kmer& rhs) noexcept;
    private:
        SequenceIterator first_, last_;
        std::size_t hash_;
    };
    
    friend bool operator==(const Kmer& lhs, const Kmer& rhs) noexcept;
    friend bool operator<(const Kmer& lhs, const Kmer& rhs) noexcept;
    
    struct KmerHash
    {
        std::size_t operator()(const Kmer& k) const noexcept { return k.hash(); }
    };
    
    struct GraphEdge
    {
        using WeightType = unsigned;
        using ScoreType  = double;
        
        GraphEdge() = default;
        
        explicit GraphEdge(WeightType weight, bool is_reference = false);
        
        WeightType weight;
        ScoreType transition_score;
        bool is_reference;
    };
    
    struct GraphNode
    {
        GraphNode() = default;
        template <typename T>
        explicit GraphNode(std::size_t index, T&& kmer, bool is_reference = false)
        : index {index}, kmer {std::forward<T>(kmer)}, is_reference {is_reference} {}
        
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
    
    using DominatorMap = std::unordered_map<Vertex, Vertex>;
    
    using Path = std::deque<Vertex>;
    using PredecessorMap = std::unordered_map<Vertex, Vertex>;
    
    static constexpr GraphEdge::ScoreType BlockedScore = 1000;
    
    unsigned k_;
    
    std::deque<Kmer> reference_kmers_;
    std::size_t reference_head_position_;
    
    KmerGraph graph_;
    
    std::unordered_map<Kmer, Vertex, KmerHash> vertex_cache_;
    
    // methods
    
    void insert_reference_into_empty_graph(const NucleotideSequence& reference);
    void insert_reference_into_populated_graph(const NucleotideSequence& reference);
    
    bool contains_kmer(const Kmer& kmer) const noexcept;
    std::size_t count_kmer(const Kmer& kmer) const noexcept;
    std::size_t reference_size() const noexcept;
    
    void regenerate_vertex_indices();
    bool is_reference_unique_path() const;
    
    Vertex null_vertex() const;
    boost::optional<Vertex> add_vertex(const Kmer& kmer, bool is_reference = false);
    void remove_vertex(Vertex v);
    void clear_and_remove_vertex(Vertex v);
    void clear_and_remove_all(const std::unordered_set<Vertex>& vertices);
    void add_edge(Vertex u, Vertex v, GraphEdge::WeightType weight, bool is_reference = false);
    void add_reference_edge(Vertex u, Vertex v);
    void remove_edge(Vertex u, Vertex v);
    void remove_edge(Edge e);
    void increment_weight(Edge e);
    void set_vertex_reference(Vertex v);
    void set_vertex_reference(const Kmer& kmer);
    void set_edge_reference(Edge e);
    const Kmer& kmer_of(Vertex v) const;
    char front_base_of(Vertex v) const;
    char back_base_of(Vertex v) const;
    const Kmer& source_kmer_of(Edge e) const;
    const Kmer& target_kmer_of(Edge e) const;
    bool is_reference(Vertex v) const;
    bool is_source_reference(Edge e) const;
    bool is_target_reference(Edge e) const;
    bool is_reference(Edge e) const;
    bool is_reference_empty() const noexcept;
    Vertex reference_head() const;
    Vertex reference_tail() const;
    Vertex next_reference(Vertex u) const;
    Vertex prev_reference(Vertex v) const;
    std::size_t num_reference_kmers() const;
    NucleotideSequence make_sequence(const Path& path) const;
    NucleotideSequence make_reference(Vertex from, Vertex to) const;
    void remove_path(const Path& path);
    bool is_bridge(Vertex v) const;
    Path::const_iterator is_bridge_until(Path::const_iterator first, Path::const_iterator last) const;
    Path::const_iterator is_bridge_until(const Path& path) const;
    bool is_bridge(Path::const_iterator first, Path::const_iterator last) const;
    bool is_bridge(const Path& path) const;
    bool joins_reference_only(Vertex v) const;
    bool is_trivial_cycle(Edge e) const;
    bool graph_has_trivial_cycle() const;
    bool is_simple_deletion(Edge e) const;
    bool is_on_path(Edge e, const Path& path) const;
    
    void remove_trivial_nonreference_cycles();
    void remove_low_weight_edges(unsigned min_weight);
    void remove_disconnected_vertices();
    std::unordered_set<Vertex> find_reachable_kmers(Vertex from) const;
    void remove_vertices_that_cant_be_reached_from(Vertex v);
    void remove_vertices_that_cant_reach(Vertex v);
    void remove_vertices_past(Vertex v);
    bool can_prune_reference_flanks() const;
    void prune_reference_flanks();
    
    std::pair<Vertex, unsigned> find_bifurcation(Vertex from, Vertex to) const;
    
    DominatorMap build_dominator_tree(Vertex from) const;
    std::unordered_set<Vertex> extract_nondominants(Vertex from) const;
    std::deque<Vertex> extract_nondominant_reference(const DominatorMap&) const;
    
    void set_out_edge_transition_scores(Vertex v);
    void set_all_edge_transition_scores_from(Vertex src);
    void set_all_in_edge_transition_scores(Vertex v, GraphEdge::ScoreType score);
    bool is_blocked(Edge e) const;
    void block_edge(Edge e);
    void block_all_in_edges(Vertex v);
    bool all_in_edges_are_blocked(Vertex v) const;
    void block_all_vertices(const std::deque<Vertex>& vertices);
    bool all_vertices_are_blocked(const std::deque<Vertex>& vertices) const;
    
    PredecessorMap find_shortest_scoring_paths(Vertex from) const;
    
    bool is_on_path(Vertex v, const PredecessorMap& predecessors, Vertex from) const;
    bool is_on_path(Edge e, const PredecessorMap& predecessors, Vertex from) const;
    Path extract_full_path(const PredecessorMap& predecessors, Vertex from) const;
    std::tuple<Assembler::Vertex, Assembler::Vertex, unsigned>
    backtrack_until_nonreference(const PredecessorMap& predecessors, Vertex from) const;
    Path extract_nonreference_path(const PredecessorMap& predecessors, Vertex from) const;
    
    std::deque<Variant> extract_k_highest_scoring_bubble_paths(unsigned k);
    
    // for debug
    friend std::ostream& operator<<(std::ostream& os, const Kmer& kmer);
    
    void print_reference_head() const;
    void print_reference_tail() const;
    void print(Edge e) const;
    void print(const Path& path) const;
    void print_dominator_tree() const;
    
    friend struct boost::property_map<KmerGraph, boost::vertex_index_t>;
    template <typename G>
    friend decltype(auto) boost::get(boost::vertex_index_t, G&);
    template <typename G>
    friend decltype(auto) boost::get(boost::vertex_index_t, const G&);
};

template <typename S1, typename S2>
Assembler::Variant::Variant(std::size_t pos, S1&& ref, S2&& alt)
:
begin_pos {pos},
ref {std::forward<S1>(ref)},
alt {std::forward<S2>(alt)}
{}

bool operator==(const Assembler::Variant& lhs, const Assembler::Variant& rhs);
bool operator!=(const Assembler::Variant& lhs, const Assembler::Variant& rhs);

} // namespace coretools
} // namespace octopus

// Hack to make some older boost algorithms work with bundled properties
namespace boost {
    using Assembler = octopus::coretools::Assembler;
    
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
