//
//  assembler.cpp
//  Octopus
//
//  Created by Daniel Cooke on 20/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "assembler.h"

#include <algorithm>
#include <set>

Assembler::Assembler(unsigned k)
:k_ {k},
 the_graph_ {},
 added_kmer_suffixes_and_prefixes_ {}
{}

void Assembler::add_reference_contig(GenomicRegion the_region, const std::string& a_contig)
{
    add_sequence(a_contig, Source::Reference);
}

void Assembler::add_read(const AlignedRead& a_read)
{
    add_sequence(a_read.get_sequence(), Source::Read);
}

void Assembler::add_sequence(const std::string& the_sequence, Source the_source)
{
    if (the_sequence.size() < k_) return;
    std::string kmer(k_, '\0');
    unsigned num_kmers = static_cast<unsigned>(the_sequence.size()) - (k_ - 1);
    for (unsigned i = 0; i < num_kmers; ++i) {
        add_kmer(the_sequence.substr(i, k_), the_source);
    }
}

void Assembler::add_kmer(std::string the_kmer, Source the_source)
{
    auto new_edge_vertices = get_vertices(the_kmer);
    if (!merge_with_existing_edge(new_edge_vertices.first, new_edge_vertices.second, the_source)) {
        add_edge(new_edge_vertices.first, new_edge_vertices.second, std::move(the_kmer), the_source);
    }
}

std::pair<Assembler::Vertex, Assembler::Vertex>
Assembler::get_vertices(const std::string& a_kmer)
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

std::pair<Assembler::Vertex, bool>
Assembler::get_vertex(const std::string& a_kmer_prefix_or_suffix) const
{
    EdgePair ep;
    for (ep = edges(the_graph_); ep.first != ep.second; ++ep.first) {
        if (is_prefix(a_kmer_prefix_or_suffix, the_graph_[*ep.first].the_kmer)) {
            return {boost::source(*ep.first, the_graph_), true};
        }
        if (is_suffix(a_kmer_prefix_or_suffix, the_graph_[*ep.first].the_kmer)) {
            return {boost::target(*ep.first, the_graph_), true};
        }
    }
    return {boost::source(*ep.first, the_graph_), false};
}

Assembler::Vertex Assembler::add_vertex(const std::string& a_k_minus_1_mer)
{
    added_kmer_suffixes_and_prefixes_.emplace(a_k_minus_1_mer);
    return boost::add_vertex(the_graph_);
}

bool Assembler::merge_with_existing_edge(Vertex& source, Vertex& target, Source the_source)
{
    auto parrallel_edge_range = boost::edge_range(source, target, the_graph_);
    for (auto begin = parrallel_edge_range.first; begin != parrallel_edge_range.second; ++begin) {
        if (is_from_another_source(the_source, the_graph_[*begin].the_source)) {
            the_graph_[*begin].the_source = Source::ReferenceAndRead;
            ++the_graph_[*begin].weight;
            return true;
        }
    }
    return false;
}

void Assembler::add_edge(Vertex& source, Vertex& target, std::string&& the_kmer, Source the_source)
{
    auto a_new_edge = boost::add_edge(source, target, the_graph_).first;
    the_graph_[a_new_edge].the_kmer   = std::move(the_kmer);
    the_graph_[a_new_edge].the_source = the_source;
    the_graph_[a_new_edge].weight     = 1; // TODO: use functor of inputs
}

bool Assembler::is_in_graph(EdgeIterator an_edge) const
{
    return an_edge != boost::edges(the_graph_).second;
}

bool Assembler::is_in_graph(const std::string& a_k_minus_1_mer) const
{
    return added_kmer_suffixes_and_prefixes_.count(a_k_minus_1_mer) > 0;
}

unsigned Assembler::num_parallel_edges(const OutEdgePair& an_edge_range) const noexcept
{
    return static_cast<unsigned>(std::distance(an_edge_range.first, an_edge_range.second));
}

bool Assembler::is_single_edge(const OutEdgePair& an_edge_range) const noexcept
{
    return num_parallel_edges(an_edge_range) == 1;
}

bool Assembler::is_from_another_source(Source lhs, Source rhs) const noexcept
{
    return (lhs == Source::Read && rhs == Source::Reference) ||
            (lhs == Source::Reference && rhs == Source::Read);
}

unsigned Assembler::get_next_vertex_index() const
{
    return static_cast<unsigned>(boost::num_vertices(the_graph_)) - 1;
}

std::vector<std::string> Assembler::get_contigs()
{
    return get_all_euler_paths(boost::vertex(0, the_graph_));
}

std::vector<std::string> Assembler::get_all_euler_paths(Vertex the_source)
{
    std::vector<std::string> euler_paths {};
    std::set<Edge> visited_edges {};
    
    std::string path {};
    OutEdgeIterator edge_range_begin, edge_range_end;
    std::tie(edge_range_begin, edge_range_end) = boost::out_edges(the_source, the_graph_);
    while (true) {
        for (; edge_range_begin != edge_range_end; ++edge_range_begin) {
            if (visited_edges.count(*edge_range_begin) == 0) break;
        }
        if (edge_range_begin == edge_range_end) break;
        visited_edges.insert(*edge_range_begin);
        if (path.empty()) {
            path += the_graph_[*edge_range_begin].the_kmer;
        } else {
            path.push_back(the_graph_[*edge_range_begin].the_kmer.back());
        }
        the_source = boost::target(*edge_range_begin, the_graph_);
        std::tie(edge_range_begin, edge_range_end) = boost::out_edges(the_source, the_graph_);
    }
    euler_paths.push_back(path);
    
    return euler_paths;
}

std::vector<Variant> Assembler::get_variants()
{
    std::vector<Variant> result {};
    Vertex the_source = boost::vertex(0, the_graph_);
    
    std::set<Edge> visited_edges {};
    OutEdgeIterator edge_range_begin, edge_range_end;
    std::tie(edge_range_begin, edge_range_end) = boost::out_edges(the_source, the_graph_);
    while (true) {
        for (; edge_range_begin != edge_range_end; ++edge_range_begin) {
            if (visited_edges.count(*edge_range_begin) == 0 &&
                the_graph_[*edge_range_begin].the_source != Source::ReferenceAndRead) break;
        }
        if (edge_range_begin == edge_range_end) break;
        visited_edges.insert(*edge_range_begin);
        
        
        the_source = boost::target(*edge_range_begin, the_graph_);
        std::tie(edge_range_begin, edge_range_end) = boost::out_edges(the_source, the_graph_);
    }
    
    return result;
}

void Assembler::clear()
{
    the_graph_.clear();
    added_kmer_suffixes_and_prefixes_.clear();
}

std::string Assembler::get_prefix(const std::string& a_kmer) const
{
    return a_kmer.substr(0, k_ - 1);
}

std::string Assembler::get_suffix(const std::string& a_kmer) const
{
    return a_kmer.substr(1, k_ - 1);
}

bool Assembler::is_prefix(const std::string& lhs, const std::string& rhs) const noexcept
{
    return std::equal(cbegin(lhs), cend(lhs), cbegin(rhs));
}

bool Assembler::is_suffix(const std::string& lhs, const std::string& rhs) const noexcept
{
    return std::equal(cbegin(lhs), cend(lhs), std::next(cbegin(rhs)));
}
