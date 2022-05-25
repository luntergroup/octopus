// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef assembler_hpp
#define assembler_hpp

#include <vector>
#include <deque>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <cstddef>
#include <utility>
#include <tuple>
#include <stdexcept>
#include <iosfwd>

#include <boost/graph/adjacency_list.hpp>
#include <boost/optional.hpp>

#include "concepts/equitable.hpp"
#include "concepts/comparable.hpp"

namespace octopus { namespace coretools { class Assembler; }}

namespace boost {

template <typename G> static decltype(auto) get(vertex_index_t, G& g);
template <typename G> static decltype(auto) get(vertex_index_t, const G& g);

} // namespace boost

namespace octopus { namespace coretools {

class Assembler
{
public:
    using NucleotideSequence = std::string;
    using BaseQualityVector = std::vector<std::uint8_t>;
    enum class Direction { forward, reverse };
    using SampleID = std::size_t;

    struct Variant;
    class NonCanonicalReferenceSequence;
    class NonUniqueReferenceSequence {};
    
    struct Parameters
    {
        unsigned kmer_size;
        bool use_strand_bias = true;
        bool use_base_qualities = true;
    };

    using BubbleScoreSetter = std::function<double(std::size_t, std::size_t)>;
    
    Assembler() = delete;
    
    Assembler(Parameters params);
    Assembler(Parameters params, const NucleotideSequence& reference);
    
    Assembler(const Assembler&);
    Assembler& operator=(const Assembler&) = delete;
    Assembler(Assembler&&)                 = default;
    Assembler& operator=(Assembler&&)      = default;
    
    ~Assembler() = default;
    
    unsigned kmer_size() const noexcept;
    Parameters params() const;
    
    // Threads the given reference sequence into the graph.
    // Throws an exception if there is already reference sequence present.
    void insert_reference(const NucleotideSequence& sequence);
    
    // Threads the given read sequence into the graph
    void insert_read(const NucleotideSequence& sequence,
                     const BaseQualityVector& base_qualities,
                     Direction strand,
                     SampleID sample);
    
    // Returns the current number of unique kmers in the graph
    std::size_t num_kmers() const noexcept;
    
    bool is_empty() const noexcept;
    
    bool is_acyclic() const;
    std::vector<SampleID> find_cyclic_samples(unsigned min_weight = 0) const;
    void remove_nonreference_cycles(bool break_chains = true);
    
    // Returns true if all the kmers in the graph are in the reference sequence
    bool is_all_reference() const;
    bool is_unique_reference() const;
    
    void try_recover_dangling_branches();
    
    // Removes edges between kmers with weight less than the given value
    void prune(unsigned min_weight);
    
    void cleanup();
    
    void clear(const std::vector<SampleID>& samples);
    void clear();
    
    std::deque<Variant> extract_variants(unsigned max_bubbles, BubbleScoreSetter min_bubble_scorer);
    
    void write_dot(std::ostream& out) const;
    
private:
    class Kmer : public Comparable<Kmer>
    {
    public:
        using NucleotideSequence = Assembler::NucleotideSequence;
        using SequenceIterator   = NucleotideSequence::const_iterator;
        
        Kmer() = default;
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
        struct Data
        {
            WeightType weight, forward_strand_weight;
            int base_quality_sum = 0;
        };
        std::map<SampleID, Data> samples;
        WeightType weight, forward_strand_weight;
        int base_quality_sum = 0;
        bool is_reference = false;
        ScoreType transition_score = 0;
    };
    struct GraphNode
    {
        std::size_t index;
        Kmer kmer;
        bool is_reference = false;
    };
    
    using KmerGraph = boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, GraphNode, GraphEdge>;
    
    using Vertex = boost::graph_traits<KmerGraph>::vertex_descriptor;
    using Edge   = boost::graph_traits<KmerGraph>::edge_descriptor;
    
    using VertexIterator = boost::graph_traits<KmerGraph>::vertex_iterator;
    using EdgeIterator   = boost::graph_traits<KmerGraph>::edge_iterator;
    
    using DominatorMap = std::unordered_map<Vertex, Vertex>;
    
    using Path = std::deque<Vertex>;
    using EdgePath = std::vector<Edge>;
    using PredecessorMap = std::unordered_map<Vertex, Vertex>;
    
    struct SubGraph
    {
        Vertex head, tail;
        std::size_t reference_offset;
    };
    
    using PathWeightDistribution = std::vector<float>;
    
    struct PathWeightStats
    {
        PathWeightDistribution distribution;
        unsigned total, total_forward, total_reverse, min, max;
        double mean, median, stdev;
    };

    struct EdgeHash
    {
        std::size_t operator()(const Edge&) const;
        EdgeHash(const Assembler& assembler) : assembler_ {assembler} {}
    private:
        const Assembler& assembler_;
    };
    using EdgeSet = std::unordered_set<Edge, EdgeHash>;
    
    Parameters params_;
    
    std::deque<Kmer> reference_kmers_;
    std::size_t reference_head_position_;
    
    KmerGraph graph_;
    
    std::unordered_map<Kmer, Vertex, KmerHash> vertex_cache_;
    Path reference_vertices_;
    std::deque<Edge> reference_edges_;
    
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
    Edge add_edge(Vertex u, Vertex v,
                  GraphEdge::WeightType weight, GraphEdge::WeightType forward_weight,
                  int base_quality_sum,
                  bool is_reference = false, 
                  boost::optional<SampleID> sample = boost::none);
    Edge add_reference_edge(Vertex u, Vertex v);
    void remove_edge(Vertex u, Vertex v);
    void remove_edge(Edge e);
    void increment_weight(Edge e, bool is_forward, int base_quality, SampleID sample);
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
    bool is_artificial(Edge e) const;
    bool is_reference_empty() const noexcept;
    Vertex reference_head() const;
    Vertex reference_tail() const;
    Vertex next_reference(Vertex u) const;
    Vertex prev_reference(Vertex v) const;
    bool is_dangling_branch(Vertex v) const;
    boost::optional<Vertex> find_joining_kmer(Vertex v) const;
    std::size_t num_reference_kmers() const;
    NucleotideSequence make_sequence(const Path& path) const;
    NucleotideSequence make_reference(Vertex from, Vertex to) const;
    void remove_path(const Path& path);
    bool is_bridge(Vertex v) const;
    bool is_reference_bridge(Vertex v) const;
    Path::const_iterator is_bridge_until(Path::const_iterator first, Path::const_iterator last) const;
    Path::const_iterator is_bridge_until(const Path& path) const;
    bool is_bridge(Path::const_iterator first, Path::const_iterator last) const;
    bool is_bridge(const Path& path) const;
    std::pair<bool, Vertex> is_bridge_to_reference(Vertex from) const;
    bool joins_reference_only(Vertex v) const;
    bool joins_reference_only(Path::const_iterator first, Path::const_iterator last) const;
    bool is_trivial_cycle(Edge e) const;
    bool graph_has_trivial_cycle() const;
    bool graph_has_nontrivial_cycle() const;
    void remove_trivial_nonreference_cycles();
    void remove_nontrivial_nonreference_cycles();
    void remove_all_nonreference_cycles(bool break_chains);
    bool is_simple_deletion(Edge e) const;
    bool is_on_path(Edge e, const Path& path) const;
    bool connects_to_path(Edge e, const Path& path) const;
    bool is_dependent_on_path(Edge e, const Path& path) const;
    GraphEdge::WeightType weight(const Path& path) const;
    GraphEdge::WeightType max_weight(const Path& path, const EdgeSet& scored_edges) const;
    unsigned count_low_weights(const Path& path, unsigned low_weight) const;
    bool has_low_weight_flanks(const Path& path, unsigned low_weight) const;
    unsigned count_low_weight_flanks(const Path& path, unsigned low_weight) const;
    template <typename UnaryFunction>
    void for_each_edge(const Path& path, UnaryFunction f) const;
    PathWeightStats compute_weight_stats(const Path& path, const EdgeSet& scored_edges) const;
    GraphEdge::WeightType sum_source_in_edge_weight(Edge e) const;
    GraphEdge::WeightType sum_target_out_edge_weight(Edge e) const;
    bool all_in_edges_low_weight(Vertex v, unsigned min_weight) const;
    bool all_out_edges_low_weight(Vertex v, unsigned min_weight) const;
    std::size_t low_weight_out_degree(Vertex v, unsigned min_weight) const;
    std::size_t low_weight_in_degree(Vertex v, unsigned min_weight) const;
    bool is_low_weight_source(Vertex v, unsigned min_weight) const;
    bool is_low_weight_sink(Vertex v, unsigned min_weight) const;
    bool is_low_weight(Vertex v, unsigned min_weight) const;
    void remove_low_weight_edges(unsigned min_weight);
    void remove_disconnected_vertices();
    std::unordered_set<Vertex> find_reachable_kmers(Vertex from) const;
    std::deque<Vertex> remove_vertices_that_cant_be_reached_from(Vertex v);
    void remove_vertices_that_cant_reach(Vertex v);
    void remove_vertices_past(Vertex v);
    bool can_prune_reference_flanks() const;
    void pop_reference_head();
    void pop_reference_tail();
    void prune_reference_flanks();
    std::pair<Vertex, unsigned> find_bifurcation(Vertex from, Vertex to) const;
    DominatorMap build_dominator_tree(Vertex from) const;
    std::unordered_set<Vertex> extract_nondominants(Vertex from) const;
    std::deque<Vertex> extract_nondominant_reference(const DominatorMap&) const;
    void set_out_edge_transition_scores(Vertex v);
    void set_all_edge_transition_scores_from(Vertex src);
    void set_all_in_edge_transition_scores(Vertex v, GraphEdge::ScoreType score);
    void block_all_in_edges(Vertex v);
    PredecessorMap find_shortest_scoring_paths(Vertex from, bool use_weights = false) const;
    bool is_on_path(Vertex v, const PredecessorMap& predecessors, Vertex from) const;
    bool is_on_path(Edge e, const PredecessorMap& predecessors, Vertex from) const;
    Path extract_full_path(const PredecessorMap& predecessors, Vertex from) const;
    std::tuple<Assembler::Vertex, Assembler::Vertex, unsigned>
    backtrack_until_nonreference(const PredecessorMap& predecessors, Vertex from) const;
    Path extract_nonreference_path(const PredecessorMap& predecessors, Vertex from) const;
    std::vector<EdgePath> extract_k_shortest_paths(Vertex src, Vertex dst, unsigned k) const;
    Edge head_edge(const Path& path) const;
    int head_mean_base_quality(const Path& path) const;
    int tail_mean_base_quality(const Path& path) const;
    double get_min_bubble_score(Vertex ref_head, Vertex ref_tail, BubbleScoreSetter min_bubble_scorer) const;
    double bubble_score(const Path& path, const EdgeSet& exclude) const;
    std::deque<Variant> extract_bubble_paths(unsigned max_bubbles, BubbleScoreSetter min_bubble_scorer);
    std::deque<SubGraph> find_independent_subgraphs() const;
    std::deque<Variant> extract_bubble_paths_with_ksp(unsigned k, BubbleScoreSetter min_bubble_scorer, EdgeSet& scored_edges);
    
    // for debug
    
    friend std::ostream& operator<<(std::ostream& os, const Kmer& kmer);
    void print_reference_head() const;
    void print_reference_tail() const;
    void print_reference_path() const;
    void print(Edge e) const;
    void print(const Path& path) const;
    void print_verbose(const Path& path) const;
    void print_dominator_tree() const;
    
    friend struct boost::property_map<KmerGraph, boost::vertex_index_t>;
    template <typename G>
    friend decltype(auto) boost::get(boost::vertex_index_t, G&);
    template <typename G>
    friend decltype(auto) boost::get(boost::vertex_index_t, const G&);
};

struct Assembler::Variant : public Equitable<Variant>
{
    Variant() = default;
    template <typename S1, typename S2>
    Variant(std::size_t pos, S1&& ref, S2&& alt);
    NucleotideSequence ref, alt;
    std::size_t begin_pos;
};

class Assembler::NonCanonicalReferenceSequence : public std::invalid_argument
{
public:
    NonCanonicalReferenceSequence(NucleotideSequence reference_sequence);
    ~NonCanonicalReferenceSequence() noexcept = default;
    const char* what() const noexcept override;
private:
    NucleotideSequence reference_sequence_;
};

template <typename S1, typename S2>
Assembler::Variant::Variant(std::size_t pos, S1&& ref, S2&& alt)
: ref {std::forward<S1>(ref)}
, alt {std::forward<S2>(alt)}
, begin_pos {pos}
{}

bool operator==(const Assembler::Variant& lhs, const Assembler::Variant& rhs) noexcept;

template <typename UnaryFunction>
void Assembler::for_each_edge(const Path& path, UnaryFunction f) const
{
    for (std::size_t u {0}, v {1}; v < path.size(); ++u, ++v) {
        Edge e; bool good;
        std::tie(e, good) = boost::edge(path[u], path[v], graph_);
        assert(good);
        f(e);
    }
}

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

#endif
