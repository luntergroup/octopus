//
//  haplotype_tree.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_tree.hpp"

#include <deque>
#include <stack>
#include <stdexcept>
#include <cassert>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/copy.hpp>

#include "reference_genome.hpp"
#include "mappable_algorithms.hpp"

namespace octopus { namespace coretools {

HaplotypeTree::HaplotypeTree(const GenomicRegion::ContigName& contig, const ReferenceGenome& reference)
:
reference_ {reference},
tree_ {},
root_ {boost::add_vertex(tree_)},
haplotype_leafs_ {root_},
contig_ {contig},
haplotype_leaf_cache_ {}
{
    if (!reference.has_contig(contig)) {
        throw std::invalid_argument {"HaplotypeTree: constructed with contig "
            + contig + " which is not in the reference " + reference.name()};
    }
}

namespace debug {
template <typename G, typename V, typename Container>
bool is_tree(const G& graph, const V& root, const Container& leafs)
{
    std::unordered_set<V> visited_vertices {};
    visited_vertices.reserve(boost::num_vertices(graph));
    
    auto vis = boost::make_bfs_visitor(boost::write_property(boost::typed_identity_property_map<V>(),
                                                             std::inserter(visited_vertices,
                                                                           std::begin(visited_vertices)),
                                                             boost::on_discover_vertex()));
    
    std::unordered_map<V, std::size_t> index_map {};
    index_map.reserve(boost::num_vertices(graph));
    
    const auto p = boost::vertices(graph);
    
    std::size_t i {0};
    std::for_each(p.first, p.second, [&i, &index_map] (const auto& v) { index_map.emplace(v, i++); } );
    
    boost::breadth_first_search(graph, root,
                                boost::visitor(vis)
                                .vertex_index_map(boost::make_assoc_property_map(index_map)));
    
    if (visited_vertices.size() != boost::num_vertices(graph)) {
        return false;
    }
    
    assert(visited_vertices.count(root) == 1);
    assert(std::all_of(std::cbegin(leafs), std::cend(leafs),
                       [&visited_vertices] (const auto& v) {
                           return visited_vertices.count(v) == 1;
                       }));
    
    return std::all_of(std::cbegin(leafs), std::cend(leafs),
                       [&graph, &root] (const auto& v) {
                           return v == root
                           || (boost::in_degree(v, graph) == 1 && boost::out_degree(v, graph) == 0);
                       });
}
} // namespace debug

HaplotypeTree::HaplotypeTree(const HaplotypeTree& other)
:
reference_ {other.reference_},
tree_ {},
root_ {},
haplotype_leafs_ {},
contig_ {other.contig_},
haplotype_leaf_cache_ {}
{
    std::unordered_map<Vertex, std::size_t> index_map {};
    index_map.reserve(boost::num_vertices(other.tree_));
    
    const auto p = boost::vertices(other.tree_);
    
    std::size_t i {0};
    std::for_each(p.first, p.second, [&i, &index_map] (const auto& v) { index_map.emplace(v, i++); } );
    
    std::unordered_map<Vertex, Vertex> vertex_copy_map {};
    vertex_copy_map.reserve(boost::num_vertices(other.tree_));
    
    boost::copy_graph(other.tree_, tree_,
                      boost::vertex_index_map(boost::make_assoc_property_map(index_map))
                      .orig_to_copy(boost::make_assoc_property_map(vertex_copy_map)));
    
    assert(vertex_copy_map.size() == boost::num_vertices(other.tree_));
    
    root_ = vertex_copy_map.at(other.root_);
    
    std::transform(std::cbegin(other.haplotype_leafs_), std::cend(other.haplotype_leafs_),
                   std::back_inserter(haplotype_leafs_),
                   [&vertex_copy_map] (const Vertex& v) {
                       return vertex_copy_map.at(v);
                   });
}

HaplotypeTree& HaplotypeTree::operator=(HaplotypeTree other)
{
    using std::swap;
    swap(reference_, other.reference_);
    swap(tree_, other.tree_);
    swap(root_, other.root_);
    swap(haplotype_leafs_, other.haplotype_leafs_);
    swap(contig_, other.contig_);
    swap(haplotype_leaf_cache_, other.haplotype_leaf_cache_);
    return *this;
}

bool HaplotypeTree::is_empty() const noexcept
{
    return haplotype_leafs_.front() == root_;
}

unsigned HaplotypeTree::num_haplotypes() const noexcept
{
    return (is_empty()) ? 0 : static_cast<unsigned>(haplotype_leafs_.size());
}

bool HaplotypeTree::contains(const Haplotype& haplotype) const
{
    if (haplotype_leaf_cache_.count(haplotype) > 0) return true;
    
    return std::any_of(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                       [this, &haplotype] (const Vertex leaf) {
                           return is_branch_equal_haplotype(leaf, haplotype);
                       });
}
    
bool HaplotypeTree::includes(const Haplotype& haplotype) const
{
    if (haplotype_leaf_cache_.count(haplotype) > 0) return true;
    
    return std::any_of(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                       [this, &haplotype] (const Vertex leaf) {
                           return is_branch_exact_haplotype(leaf, haplotype);
                       });
}

bool HaplotypeTree::is_unique(const Haplotype& haplotype) const
{
    if (haplotype_leaf_cache_.count(haplotype) > 0) {
        return haplotype_leaf_cache_.count(haplotype) == 1;
    }
    
    bool haplotype_seen {false};
    
    for (const Vertex& leaf : haplotype_leafs_) {
        if (is_branch_equal_haplotype(leaf, haplotype)) {
            if (haplotype_seen) {
                return false;
            } else {
                haplotype_seen = true;
            }
        }
    }
    
    return haplotype_seen;
}

HaplotypeTree& HaplotypeTree::extend(const ContigAllele& allele)
{
    for (auto it = std::cbegin(haplotype_leafs_), end = std::cend(haplotype_leafs_); it != end; ++it) {
        it = extend_haplotype(it, allele);
    }
    haplotype_leaf_cache_.clear();
    return *this;
}

HaplotypeTree& HaplotypeTree::extend(const Allele& allele)
{
    if (contig_name(allele) != contig_) {
        throw std::domain_error {"HaplotypeTree: trying to extend with Allele on different contig"};
    }
    return extend(demote(allele));
}

template <typename Container, typename V>
struct Splicer : public boost::default_dfs_visitor
{
    Splicer(const ContigAllele& allele, std::stack<V>& candidate_splice_sites, Container& splice_sites,
            V root)
    : allele_ {allele}, candidate_splice_sites_ {candidate_splice_sites},
    splice_sites_ {splice_sites}, root_ {root} {}
    
    template <typename G>
    void finish_vertex(const V v, const G& tree)
    {
        if (!candidate_splice_sites_.empty() && v == candidate_splice_sites_.top()) {
            candidate_splice_sites_.pop();
            if (v == root_ || is_after(allele_.get(), tree[v])) {
                splice_sites_.push_back(v);
            } else {
                const auto u = *boost::inv_adjacent_vertices(v, tree).first;
                if (candidate_splice_sites_.empty() || candidate_splice_sites_.top() != u) {
                    candidate_splice_sites_.push(u);
                }
            }
        }
    }
private:
    std::reference_wrapper<const ContigAllele> allele_;
    std::stack<V>& candidate_splice_sites_;
    Container& splice_sites_;
    V root_;
};

template <typename Container, typename V>
auto make_splicer(const ContigAllele& allele, std::stack<V>& candidate_splice_sites,
                  Container& splice_sites, V root)
{
    return Splicer<Container, V> {allele, candidate_splice_sites, splice_sites, root};
}

void HaplotypeTree::splice(const ContigAllele& allele)
{
    if (is_empty()) {
        extend(allele);
        return;
    }
    
    std::unordered_map<Vertex, boost::default_color_type> colours {};
    colours.reserve(boost::num_vertices(tree_));
    
    std::deque<Vertex> splice_sites {};
    std::stack<Vertex> candidate_splice_sites {};
    
    boost::depth_first_visit(tree_, root_,
                             make_splicer(allele, candidate_splice_sites, splice_sites, root_),
                             boost::make_assoc_property_map(colours),
                             [&] (const Vertex v, const Tree& tree) -> bool {
                                 if (v != root_ && begins_before(allele, tree[v])) {
                                     const auto u = *boost::inv_adjacent_vertices(v, tree).first;
                                     if (candidate_splice_sites.empty() || candidate_splice_sites.top() != u) {
                                         candidate_splice_sites.push(u);
                                     }
                                     return true;
                                 }
                                 return false;
                             });
    
    assert(candidate_splice_sites.empty());
    
    for (const auto v : splice_sites) {
        const auto spliced = boost::add_vertex(allele, tree_);
        boost::add_edge(v, spliced, tree_);
        haplotype_leafs_.push_back(spliced);
    }
}

void HaplotypeTree::splice(const Allele& allele)
{
    if (contig_name(allele) != contig_) {
        throw std::domain_error {"HaplotypeTree: trying to splicing with Allele on different contig"};
    }
    return splice(demote(allele));
}

GenomicRegion HaplotypeTree::encompassing_region() const
{
    if (is_empty()) {
        throw std::runtime_error {"HaplotypeTree::encompassing_region called on empty tree"};
    }
    
    const auto vertex_range = boost::adjacent_vertices(root_, tree_);
    
    const auto leftmost = *std::min_element(vertex_range.first, vertex_range.second,
                                            [this] (const auto& lhs, const auto& rhs) {
                                                return begins_before(tree_[lhs], tree_[rhs]);
                                            });
    
    const auto rightmost = *std::max_element(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                                             [this] (const auto& lhs, const auto& rhs) {
                                                 return ends_before(tree_[lhs], tree_[rhs]);
                                             });
    
    return GenomicRegion {contig_, ::encompassing_region(tree_[leftmost], tree_[rightmost])};
}

std::vector<Haplotype> HaplotypeTree::extract_haplotypes() const
{
    if (is_empty()) return std::vector<Haplotype> {};
    return extract_haplotypes(encompassing_region());
}

std::vector<Haplotype> HaplotypeTree::extract_haplotypes(const GenomicRegion& region) const
{
    haplotype_leaf_cache_.clear();
    haplotype_leaf_cache_.reserve(num_haplotypes());
    
    std::vector<Haplotype> result {};
    
    if (is_empty() || !overlaps(region, encompassing_region())) return result;
    
    result.reserve(num_haplotypes());
    
    for (const auto leaf : haplotype_leafs_) {
        auto haplotype = extract_haplotype(leaf, region);
        
        // recently retreived haplotypes are added to the cache as it is likely these
        // are the haplotypes that will be pruned next
        haplotype_leaf_cache_.emplace(haplotype, leaf);
        
        result.push_back(std::move(haplotype));
    }
    
    return result;
}

void HaplotypeTree::prune_all(const Haplotype& haplotype)
{
    using std::cbegin; using std::cend; using std::for_each; using std::find;
    
    if (is_empty() || contig_name(haplotype) != contig_) return;
    
    // If any of the haplotypes in cache match the query haplotype then the cache must contain
    // all possible leaves corrosponding to that haplotype. So we don't need to look through
    // the list of all leaves. Win.
    if (haplotype_leaf_cache_.count(haplotype) > 0) {
        const auto possible_leafs = haplotype_leaf_cache_.equal_range(haplotype);
        
        for_each(possible_leafs.first, possible_leafs.second,
                 [this, &haplotype] (const auto& leaf_pair) {
                     const auto p = clear(leaf_pair.second, contig_region(haplotype));
                     
                     auto leaf_itr = find(cbegin(haplotype_leafs_), cend(haplotype_leafs_),
                                          leaf_pair.second);
                     
                     leaf_itr = haplotype_leafs_.erase(leaf_itr);
                     
                     if (p.second) {
                         haplotype_leafs_.insert(leaf_itr, p.first);
                     }
                 });
        
        haplotype_leaf_cache_.erase(haplotype);
    } else {
        auto leaf_itr = cbegin(haplotype_leafs_);
        
        while (true) {
            leaf_itr = find_equal_haplotype_leaf(leaf_itr, cend(haplotype_leafs_), haplotype);
            
            if (leaf_itr == cend(haplotype_leafs_)) return;
            
            const auto p = clear(*leaf_itr, contig_region(haplotype));
            
            leaf_itr = haplotype_leafs_.erase(leaf_itr);
            
            if (p.second) {
                leaf_itr = haplotype_leafs_.insert(leaf_itr, p.first);
            }
        }
    }
}

void HaplotypeTree::prune_unique(const Haplotype& haplotype)
{
    using std::cbegin; using std::cend; using std::for_each;
    
    if (is_empty()) return;
    
    if (haplotype_leaf_cache_.count(haplotype) > 0) {
        const auto possible_leafs = haplotype_leaf_cache_.equal_range(haplotype);
        
        const auto it = std::find_if(possible_leafs.first, possible_leafs.second,
                                     [this, &haplotype] (const auto& leaf_pair) {
                                         return is_branch_exact_haplotype(leaf_pair.second, haplotype);
                                     });
        
        if (it == possible_leafs.second) {
            throw std::runtime_error {"HaplotypeTree::prune_unique called with matching Haplotype not in tree"};
        }
        
        const auto leaf_to_keep_itr = it->second;
        
        std::for_each(possible_leafs.first, possible_leafs.second,
                      [this, &haplotype, leaf_to_keep_itr] (auto& leaf_pair) {
                          if (leaf_pair.second != leaf_to_keep_itr) {
                              const auto p = clear(leaf_pair.second, contig_region(haplotype));
                              
                              auto leaf_itr = std::find(cbegin(haplotype_leafs_), cend(haplotype_leafs_),
                                                        leaf_pair.second);
                              
                              leaf_itr = haplotype_leafs_.erase(leaf_itr);
                              
                              if (p.second) {
                                  haplotype_leafs_.insert(leaf_itr, p.first);
                              }
                          }
                      });
        
        haplotype_leaf_cache_.erase(haplotype);
        haplotype_leaf_cache_.emplace(haplotype, leaf_to_keep_itr);
    } else {
        auto leaf_itr = cbegin(haplotype_leafs_);
        
        const auto leaf_to_keep_itr = find_exact_haplotype_leaf(leaf_itr, cend(haplotype_leafs_),
                                                                haplotype);
        
        while (true) {
            leaf_itr = find_equal_haplotype_leaf(leaf_itr, cend(haplotype_leafs_), haplotype);
            
            if (leaf_itr == cend(haplotype_leafs_)) return;
            
            if (leaf_itr == leaf_to_keep_itr) {
                std::advance(leaf_itr, 1);
                continue;
            }
            
            const auto p = clear(*leaf_itr, contig_region(haplotype));
            
            leaf_itr = haplotype_leafs_.erase(leaf_itr);
            
            if (p.second) {
                leaf_itr = haplotype_leafs_.insert(leaf_itr, p.first);
            }
        }
    }
}

void HaplotypeTree::clear(const GenomicRegion& region)
{
    if (is_empty()) return;
    
    const auto tree_region = encompassing_region();
    
    if (::contains(region, tree_region)) {
        clear();
    } else if (overlaps(region, tree_region)) {
        haplotype_leaf_cache_.clear();
        
        std::list<Vertex> new_leafs {};
        
        std::for_each(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                      [this, &region, &new_leafs] (const Vertex leaf) {
                          const auto p = clear(leaf, contig_region(region));
                          if (p.second) new_leafs.push_back(p.first);
                      });
        
        haplotype_leafs_ = new_leafs;
    }
}

void HaplotypeTree::clear() noexcept
{
    haplotype_leaf_cache_.clear();
    haplotype_leafs_.clear();
    tree_.clear();
    root_ = boost::add_vertex(tree_);
    haplotype_leafs_.push_back(root_);
}

// Private methods

HaplotypeTree::Vertex HaplotypeTree::get_previous_allele(const Vertex allele) const
{
    const auto p = boost::inv_adjacent_vertices(allele, tree_);
    assert(std::distance(p.first, p.second) == 1);
    return *p.first;
}

bool HaplotypeTree::is_bifurcating(const Vertex v) const
{
    return boost::out_degree(v, tree_) > 1;
}

HaplotypeTree::Vertex HaplotypeTree::remove_forward(Vertex u)
{
    const auto p = boost::adjacent_vertices(u, tree_);
    assert(p.first != p.second);
    const auto v = *p.first;
    boost::remove_edge(u, v, tree_);
    boost::remove_vertex(u, tree_);
    return v;
}

HaplotypeTree::Vertex HaplotypeTree::remove_backward(const Vertex v)
{
    const auto u = get_previous_allele(v);
    boost::remove_edge(u, v, tree_);
    boost::remove_vertex(v, tree_);
    return u;
}

bool HaplotypeTree::allele_exists(Vertex leaf, const ContigAllele& allele) const
{
    const auto vertex_range = boost::adjacent_vertices(leaf, tree_);
    return std::any_of(vertex_range.first, vertex_range.second,
                       [this, &allele] (const Vertex vertex) {
                           return tree_[vertex] == allele;
                       });
}

HaplotypeTree::Vertex HaplotypeTree::find_allele_before(Vertex v, const ContigAllele& allele) const
{
    while (v != root_ && overlaps(allele, tree_[v])) {
        if (is_same_region(allele, tree_[v])) { // for insertions
            v = get_previous_allele(v);
            break;
        }
        v = get_previous_allele(v);
    }
    return v;
}

bool can_add_to_branch(const ContigAllele& new_allele, const ContigAllele& leaf)
{
    return !are_adjacent(leaf, new_allele)
            || !((is_insertion(leaf) && is_deletion(new_allele))
             || (is_deletion(leaf) && is_insertion(new_allele)));
}

HaplotypeTree::LeafIterator
HaplotypeTree::extend_haplotype(LeafIterator leaf_itr, const ContigAllele& new_allele)
{
    if (*leaf_itr == root_) {
        const auto new_leaf = boost::add_vertex(new_allele, tree_);
        boost::add_edge(*leaf_itr, new_leaf, tree_);
        leaf_itr = haplotype_leafs_.erase(leaf_itr);
        return haplotype_leafs_.insert(leaf_itr, new_leaf);
    }
    
    const auto& leaf_allele = tree_[*leaf_itr];
    
    if (!can_add_to_branch(new_allele, leaf_allele)) return leaf_itr;
    
    if (is_after(new_allele, leaf_allele)) {
        const auto new_leaf = boost::add_vertex(new_allele, tree_);
        boost::add_edge(*leaf_itr, new_leaf, tree_);
        leaf_itr = haplotype_leafs_.erase(leaf_itr);
        leaf_itr = haplotype_leafs_.insert(leaf_itr, new_leaf);
    } else if (overlaps(new_allele, tree_[*leaf_itr])) {
        const auto branch_point = find_allele_before(*leaf_itr, new_allele);
        if ((branch_point == root_ || can_add_to_branch(new_allele, tree_[branch_point]))
            && !allele_exists(branch_point, new_allele)) {
            const auto new_leaf = boost::add_vertex(new_allele, tree_);
            boost::add_edge(branch_point, new_leaf, tree_);
            haplotype_leafs_.insert(leaf_itr, new_leaf);
        }
    }
    
    return leaf_itr;
}

Haplotype HaplotypeTree::extract_haplotype(Vertex leaf, const GenomicRegion& region) const
{
    Haplotype::Builder result {region, reference_};
    
    const auto& contig_region = region.contig_region();
    
    while (leaf != root_ && !overlaps(tree_[leaf], contig_region)) {
        leaf = get_previous_allele(leaf);
    }
    
    while (leaf != root_ && overlaps(tree_[leaf], contig_region)) {
        result.push_front(tree_[leaf]);
        leaf = get_previous_allele(leaf);
    }
    
    return result.build();
}

bool HaplotypeTree::define_same_haplotype(Vertex leaf1, Vertex leaf2) const
{
    if (leaf1 == leaf2) return true;
    
    while (leaf1 != root_) {
        if (leaf2 == root_ || tree_[leaf1] != tree_[leaf2]) return false;
        leaf1 = get_previous_allele(leaf1);
        leaf2 = get_previous_allele(leaf2);
    }
    
    return leaf2 == root_;
}

bool HaplotypeTree::is_branch_exact_haplotype(Vertex leaf, const Haplotype& haplotype) const
{
    if (leaf == root_ || !overlaps(tree_[leaf], contig_region(haplotype))) {
        return false;
    }
    
    while (leaf != root_) {
        if (!haplotype.includes(tree_[leaf]))
            return false;
        leaf = get_previous_allele(leaf);
    }
    
    return true;
}

bool HaplotypeTree::is_branch_equal_haplotype(const Vertex leaf, const Haplotype& haplotype) const
{
    // TODO: check if this is quicker than calling Haplotype::contains for each ContigAllele
    return leaf != root_ && overlaps(contig_region(haplotype), tree_[leaf])
            && extract_haplotype(leaf, haplotype.mapped_region()) == haplotype;
}

HaplotypeTree::LeafIterator
HaplotypeTree::find_exact_haplotype_leaf(const LeafIterator first, const LeafIterator last,
                                         const Haplotype& haplotype) const
{
    return std::find_if(first, last,
                        [this, &haplotype] (Vertex leaf) {
                            return is_branch_exact_haplotype(leaf, haplotype);
                        });
}

HaplotypeTree::LeafIterator
HaplotypeTree::find_equal_haplotype_leaf(const LeafIterator first, const LeafIterator last,
                                         const Haplotype& haplotype) const
{
    return std::find_if(first, last,
                        [this, &haplotype] (Vertex leaf) {
                            return is_branch_equal_haplotype(leaf, haplotype);
                        });
}

std::pair<HaplotypeTree::Vertex, bool>
HaplotypeTree::clear(Vertex leaf, const ContigRegion& region)
{
    if (ends_before(region, tree_[leaf])) {
        return clear_internal(leaf, region);
    } else {
        return clear_external(leaf, region);
    }
}

std::pair<HaplotypeTree::Vertex, bool>
HaplotypeTree::clear_external(Vertex leaf, const ContigRegion& region)
{
    assert(boost::out_degree(leaf, tree_) == 0);
    
    while (leaf != root_) {
        if (boost::out_degree(leaf, tree_) > 0) {
            return std::make_pair(leaf, false);
        } else if (begins_before(tree_[leaf], region)) {
            return std::make_pair(leaf, true);
        } else {
            leaf = remove_backward(leaf);
        }
    }
    
    // the root should only be indicated as a leaf node if there are no other nodes in the tree
    return std::make_pair(leaf, boost::num_vertices(tree_) == 1);
}

std::pair<HaplotypeTree::Vertex, bool>
HaplotypeTree::clear_internal(const Vertex leaf, const ContigRegion& region)
{
    // TODO: we can optimise this for cases where region overlaps the leftmost alleles in the tree
    
    if (leaf == root_ || is_after(region, tree_[leaf])) {
        return std::make_pair(leaf, true);
    }
    
    auto current_allele = leaf;
    auto allele_to_move = leaf;
    std::deque<const Vertex> alleles_to_copy {};
    bool is_bifurcating_branch {false};
    
    while (true) {
        current_allele = get_previous_allele(current_allele);
        
        if (current_allele == root_ || overlaps(tree_[current_allele], region)) break;
        
        is_bifurcating_branch = is_bifurcating_branch || is_bifurcating(current_allele);
        
        if (!is_bifurcating_branch) {
            allele_to_move = current_allele;
        } else {
            alleles_to_copy.push_front(current_allele);
        }
    }
    
    if (alleles_to_copy.empty()) {
        boost::remove_edge(current_allele, allele_to_move, tree_);
    } else {
        boost::remove_edge(alleles_to_copy.back(), allele_to_move, tree_);
    }
    
    while (current_allele != root_ && overlaps(region, tree_[current_allele])) {
        const auto previous_allele = get_previous_allele(current_allele);
        
        is_bifurcating_branch = is_bifurcating_branch || boost::out_degree(current_allele, tree_) > 0;
        
        if (!is_bifurcating_branch) {
            boost::remove_edge(previous_allele, current_allele, tree_);
            boost::remove_vertex(current_allele, tree_);
        }
        
        current_allele = previous_allele;
    }
    
    std::for_each(std::crbegin(alleles_to_copy), std::crend(alleles_to_copy),
                  [this, &allele_to_move] (const Vertex allele) {
                      const auto v = boost::add_vertex(tree_[allele], tree_);
                      boost::add_edge(v, allele_to_move, tree_);
                      allele_to_move = v;
                  });
    
    alleles_to_copy.clear();
    alleles_to_copy.shrink_to_fit();
    
    auto allele_to_move_to = current_allele;
    
    while (true) {
        const auto vertex_range = boost::adjacent_vertices(allele_to_move_to, tree_);
        
        const auto it = std::find_if(vertex_range.first, vertex_range.second,
                                     [this, allele_to_move] (const Vertex allele) {
                                         return tree_[allele_to_move] == tree_[allele];
                                     });
        
        if (it == vertex_range.second) break;
        
        allele_to_move_to = *it;
        
        if (boost::out_degree(allele_to_move, tree_) == 0) break;
        
        allele_to_move = remove_forward(allele_to_move);
    }
    
    if (allele_to_move_to == root_ || tree_[allele_to_move_to] != tree_[allele_to_move]) {
        boost::add_edge(allele_to_move_to, allele_to_move, tree_);
        return std::make_pair(leaf, true);
    } else {
        while (boost::out_degree(allele_to_move, tree_) > 0) {
            allele_to_move = remove_forward(allele_to_move);
        }
        
        boost::remove_vertex(allele_to_move, tree_);
        
        return std::make_pair(allele_to_move_to, false);
    }
}

} // namespace coretools
} // namespace octopus
