//
//  haplotype_tree.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_tree.hpp"

#include <deque>
#include <stdexcept>
#include <iterator>  // std::next, std::tie, std::cbegin etc
#include <algorithm> // std::find, std::find_if, std::any_of

#include "genomic_region.hpp"
#include "reference_genome.hpp"
#include "mappable_algorithms.hpp"

namespace Octopus
{

HaplotypeTree::HaplotypeTree(ReferenceGenome& reference)
:
reference_ {reference},
tree_ {},
root_ {boost::add_vertex(tree_)},
haplotype_leafs_ {root_},
haplotype_leaf_cache_ {},
recently_removed_haplotypes_ {}
{}

bool HaplotypeTree::empty() const noexcept
{
    return haplotype_leafs_.front() == root_;
}

unsigned HaplotypeTree::num_haplotypes() const noexcept
{
    return static_cast<unsigned>(haplotype_leafs_.size());
}

bool HaplotypeTree::contains(const Haplotype& haplotype) const
{
    if (haplotype_leaf_cache_.count(haplotype) > 0) return true;
    
    return std::any_of(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                       [this, &haplotype] (const Vertex leaf) {
                           return is_branch_equal_haplotype(leaf, haplotype);
                       });
}

bool HaplotypeTree::is_unique(const Haplotype& haplotype) const
{
    if (haplotype_leaf_cache_.count(haplotype) > 0) {
        return haplotype_leaf_cache_.count(haplotype) == 1;
    } else {
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
}

void HaplotypeTree::extend(const Allele& allele)
{
    for (auto leaf_it = std::cbegin(haplotype_leafs_), end = std::cend(haplotype_leafs_); leaf_it != end; ++leaf_it) {
        leaf_it = extend_haplotype(leaf_it, allele);
    }
    
    haplotype_leaf_cache_.clear();
    recently_removed_haplotypes_.clear();
}

GenomicRegion HaplotypeTree::get_region() const
{
    if (empty()) {
        throw std::runtime_error {"empty HaplotypeTree does not have a defined region"};
    }
    
    auto vertex_range = boost::adjacent_vertices(root_, tree_);
    
    auto leftmost = *std::min_element(vertex_range.first, vertex_range.second,
                                      [this] (const auto& lhs, const auto& rhs) {
                                          return tree_[lhs] < tree_[rhs];
                                      });
    
    auto rightmost = *std::max_element(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                                       [this] (const auto& lhs, const auto& rhs) {
                                           return tree_[lhs] < tree_[rhs];
                                       });
    
    return get_encompassing(tree_[leftmost], tree_[rightmost]);
}
    
std::vector<Haplotype> HaplotypeTree::get_haplotypes() const
{
    return get_haplotypes(get_region());
}

std::vector<Haplotype> HaplotypeTree::get_haplotypes(const GenomicRegion& region) const
{
    haplotype_leaf_cache_.clear();
    haplotype_leaf_cache_.reserve(haplotype_leafs_.size());
    
    std::vector<Haplotype> result {};
    result.reserve(haplotype_leafs_.size());
    
    for (auto leaf : haplotype_leafs_) {
        auto haplotype = get_haplotype(leaf, region);
        
        // recently retreived haplotypes are added to the cache as it is likely these
        // are the haplotypes that will be pruned next
        haplotype_leaf_cache_.emplace(haplotype, leaf);
        
        result.emplace_back(std::move(haplotype));
    }
    
    haplotype_leaf_cache_.rehash(haplotype_leaf_cache_.size());
    
    result.shrink_to_fit();
    
    return result;
}

void HaplotypeTree::prune_all(const Haplotype& haplotype)
{
    if (empty()) return;
    if (recently_removed_haplotypes_.count(haplotype) > 1) return;
    
    Vertex new_haplotype_end;
    LeafIterator leaf_it;
    bool new_end_is_leaf {};
    
    // if any of the haplotypes in cache match the query haplotype then the cache must contain
    // all possible leaves corrosponding to that haplotype. So we don't need to look through
    // the list of all leaves. Win.
    if (haplotype_leaf_cache_.count(haplotype) > 0) {
        auto leaf_range = haplotype_leaf_cache_.equal_range(haplotype);
        
        std::for_each(leaf_range.first, leaf_range.second,
              [this, &haplotype, &leaf_it, &new_haplotype_end, &new_end_is_leaf] (auto& leaf_pair) {
                  std::tie(new_haplotype_end, new_end_is_leaf) = prune_branch(leaf_pair.second,
                                                                              haplotype.get_region());
                  
                  leaf_it = std::find(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                                      leaf_pair.second);
                  
                  leaf_it = haplotype_leafs_.erase(leaf_it);
                  
                  if (new_end_is_leaf) {
                      haplotype_leafs_.insert(leaf_it, new_haplotype_end);
                  }
              });
        
        haplotype_leaf_cache_.erase(haplotype);
    } else {
        leaf_it = std::cbegin(haplotype_leafs_);
        
        while (true) {
            leaf_it = find_equal_haplotype_leaf(leaf_it, std::cend(haplotype_leafs_), haplotype);
            
            if (leaf_it == std::cend(haplotype_leafs_)) return;
            
            std::tie(new_haplotype_end, new_end_is_leaf) = prune_branch(*leaf_it, haplotype.get_region());
            
            leaf_it = haplotype_leafs_.erase(leaf_it);
            
            if (new_end_is_leaf) {
                leaf_it = haplotype_leafs_.insert(leaf_it, new_haplotype_end);
            }
        }
    }
    
    recently_removed_haplotypes_.emplace(haplotype);
}

void HaplotypeTree::prune_unique(const Haplotype& haplotype)
{
    if (empty()) return;
    
    Vertex new_haplotype_end;
    Vertex leaf_to_keep;
    LeafIterator leaf_it;
    bool new_end_is_leaf {};
    
    if (haplotype_leaf_cache_.count(haplotype) > 0) {
        auto leaf_range = haplotype_leaf_cache_.equal_range(haplotype);
        
        leaf_to_keep = std::find_if(leaf_range.first, leaf_range.second,
                                    [this, &haplotype] (auto& leaf_pair) {
                                        return is_branch_exact_haplotype(leaf_pair.second, haplotype);
                                    })->second;
        
        std::for_each(leaf_range.first, leaf_range.second,
          [this, &haplotype, &leaf_it, &new_haplotype_end, &new_end_is_leaf, &leaf_to_keep] (auto& leaf_pair) {
              if (leaf_pair.second != leaf_to_keep) {
                  std::tie(new_haplotype_end, new_end_is_leaf) = prune_branch(leaf_pair.second,
                                                                              haplotype.get_region());
                  
                  leaf_it = std::find(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                                      leaf_pair.second);
                  
                  leaf_it = haplotype_leafs_.erase(leaf_it);
                  
                  if (new_end_is_leaf) {
                      haplotype_leafs_.insert(leaf_it, new_haplotype_end);
                  }
              }
          });
        
        haplotype_leaf_cache_.erase(haplotype);
        haplotype_leaf_cache_.emplace(haplotype, leaf_to_keep);
    } else {
        leaf_to_keep = *find_exact_haplotype_leaf(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                                                  haplotype);
        
        leaf_it = std::cbegin(haplotype_leafs_);
        
        while (true) {
            leaf_it = find_equal_haplotype_leaf(leaf_it, std::cend(haplotype_leafs_), haplotype);
            
            if (*leaf_it == leaf_to_keep) continue;
            
            if (leaf_it == std::cend(haplotype_leafs_)) return;
            
            std::tie(new_haplotype_end, new_end_is_leaf) = prune_branch(*leaf_it, haplotype.get_region());
            
            leaf_it = haplotype_leafs_.erase(leaf_it);
            
            if (new_end_is_leaf) {
                leaf_it = haplotype_leafs_.insert(leaf_it, new_haplotype_end);
            }
        }

    }
}

void HaplotypeTree::clear(const GenomicRegion& region)
{
    if (empty()) return;
    
    const auto tree_region = get_region();
    
    if (::contains(region, tree_region)) {
        clear();
    } else {
        haplotype_leaf_cache_.clear();
        recently_removed_haplotypes_.clear();
        
        std::deque<Vertex> new_leafs {};
        
        std::for_each(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                      [this, &region, &new_leafs] (Vertex leaf) {
                          auto p = prune_branch(leaf, region);
                          if (p.second) {
                              new_leafs.push_back(p.first);
                          }
                      });
        
        haplotype_leafs_.clear();
        
        std::deque<Vertex> leaves_to_prune {};
        
        std::for_each(std::begin(new_leafs), std::end(new_leafs),
                      [this, &leaves_to_prune] (Vertex new_leaf) {
                          auto seen_before = std::any_of(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                                                         [this, new_leaf] (Vertex leaf) {
                                                             return define_same_haplotype(new_leaf, leaf);
                                                         });
                          if (seen_before) {
                              leaves_to_prune.push_back(new_leaf);
                          } else {
                              haplotype_leafs_.push_back(new_leaf);
                          }
        });
        
        std::for_each(std::begin(leaves_to_prune), std::end(leaves_to_prune),
                      [this, &tree_region] (Vertex leaf) {
                          prune_branch(leaf, tree_region);
                      });
    }
}

void HaplotypeTree::clear()
{
    haplotype_leaf_cache_.clear();
    recently_removed_haplotypes_.clear();
    haplotype_leafs_.clear();
    tree_.clear();
    root_ = boost::add_vertex(tree_);
    haplotype_leafs_.push_back(root_);
}

// Private methods

HaplotypeTree::Vertex HaplotypeTree::get_previous_allele(Vertex allele) const
{
    return *boost::inv_adjacent_vertices(allele, tree_).first;
}

bool HaplotypeTree::allele_exists(Vertex leaf, const Allele& allele) const
{
    const auto vertex_range = boost::adjacent_vertices(leaf, tree_);
    return std::any_of(vertex_range.first, vertex_range.second,
                       [this, &allele] (Vertex vertex) {
                           return tree_[vertex] == allele;
                       });
}

HaplotypeTree::LeafIterator HaplotypeTree::extend_haplotype(LeafIterator haplotype_leaf, const Allele& new_allele)
{
    const Allele& the_leaf_allele = tree_[*haplotype_leaf];
    
    if (are_adjacent(the_leaf_allele, new_allele) &&
        ((is_insertion(the_leaf_allele) && is_deletion(new_allele)) ||
        (is_deletion(the_leaf_allele) && is_insertion(new_allele)))) {
        return haplotype_leaf;
    }
    
    if (*haplotype_leaf == root_ || is_after(new_allele, the_leaf_allele)) {
        Vertex new_leaf = boost::add_vertex(tree_);
        tree_[new_leaf] = new_allele;
        boost::add_edge(*haplotype_leaf, new_leaf, tree_);
        
        haplotype_leaf = haplotype_leafs_.erase(haplotype_leaf);
        haplotype_leaf = haplotype_leafs_.insert(haplotype_leaf, new_leaf);
    } else if (overlaps(new_allele, the_leaf_allele)) {
        Vertex previous_allele = get_previous_allele(*haplotype_leaf);
        
        const Allele& the_previous_allele = tree_[previous_allele];
        
        if (are_adjacent(the_previous_allele, new_allele) &&
            ((is_insertion(the_previous_allele) && is_deletion(new_allele)) ||
             (is_deletion(the_previous_allele) && is_insertion(new_allele)))) {
            return haplotype_leaf;
        }
        
        if (!allele_exists(previous_allele, new_allele)) {
            Vertex new_branch = boost::add_vertex(tree_);
            tree_[new_branch] = new_allele;
            boost::add_edge(previous_allele, new_branch, tree_);
            
            haplotype_leafs_.insert(haplotype_leaf, new_branch);
        }
    }
    
    return haplotype_leaf;
}

Haplotype HaplotypeTree::get_haplotype(Vertex haplotype_leaf, const GenomicRegion& region) const
{
    Haplotype result {reference_, region};
    
    while (haplotype_leaf != root_ && is_after(tree_[haplotype_leaf], region)) {
        haplotype_leaf = get_previous_allele(haplotype_leaf);
    }
    
    while (haplotype_leaf != root_ && overlaps(tree_[haplotype_leaf], region)) {
        result.push_front(tree_[haplotype_leaf]);
        haplotype_leaf = get_previous_allele(haplotype_leaf);
    }
    
    return result;
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

bool HaplotypeTree::is_branch_exact_haplotype(Vertex branch_vertex, const Haplotype& haplotype) const
{
    if (branch_vertex == root_) return false;
    
    const Allele* this_allele {&tree_[branch_vertex]};
    
    if (!overlaps(haplotype, *this_allele)) return false;
    
    Vertex previous_vertex;
    const Allele* previous_allele;
    
    while (!is_before(*this_allele, haplotype)) {
        if (!haplotype.contains(*this_allele)) return false;
        
        previous_vertex = get_previous_allele(branch_vertex);
        
        if (previous_vertex == root_) break;
        
        previous_allele = &tree_[previous_vertex];
        
        if (is_before(*previous_allele, haplotype)) break;
        
        if (!are_adjacent(*previous_allele, *this_allele) &&
            !haplotype.contains(get_reference_allele(get_intervening(*previous_allele, *this_allele),
                                                       reference_))) {
            return false;
        }
        
        branch_vertex = previous_vertex;
        this_allele = previous_allele;
    }
    
    return true;
}

bool HaplotypeTree::is_branch_equal_haplotype(Vertex branch_vertex, const Haplotype& haplotype) const
{
    return branch_vertex != root_ && overlaps(haplotype, tree_[branch_vertex]) &&
            get_haplotype(branch_vertex, haplotype.get_region()) == haplotype;
}

HaplotypeTree::LeafIterator HaplotypeTree::find_exact_haplotype_leaf(LeafIterator first, LeafIterator last,
                                                                     const Haplotype& haplotype) const
{
    return std::find_if(first, last, [this, &haplotype] (Vertex leaf) { return is_branch_exact_haplotype(leaf, haplotype); });
}

HaplotypeTree::LeafIterator HaplotypeTree::find_equal_haplotype_leaf(LeafIterator first, LeafIterator last,
                                                                     const Haplotype& haplotype) const
{
    return std::find_if(first, last, [this, &haplotype] (Vertex leaf) { return is_branch_equal_haplotype(leaf, haplotype); });
}

std::pair<HaplotypeTree::Vertex, bool> HaplotypeTree::prune_branch(Vertex leaf, const GenomicRegion& region)
{
    if (ends_before(region, tree_[leaf])) {
        return splice_region(leaf, region);
    }
    
    while (leaf != root_) {
        if (boost::out_degree(leaf, tree_) > 0) {
            return std::make_pair(leaf, false);
        } else if (begins_before(tree_[leaf], region)) {
            return std::make_pair(leaf, true);
        }
        
        auto new_leaf = get_previous_allele(leaf);
        
        boost::remove_edge(new_leaf, leaf, tree_);
        boost::remove_vertex(leaf, tree_);
        
        leaf = new_leaf;
    }
    
    // the root should only be indicated as a leaf node if there are no other nodes in the tree
    return std::make_pair(leaf, num_haplotypes() == 1);
}

std::pair<HaplotypeTree::Vertex, bool> HaplotypeTree::splice_region(const Vertex leaf, const GenomicRegion& region)
{
    if (leaf == root_ || is_after(region, tree_[leaf])) {
        return std::make_pair(leaf, true);
    }
    
    //std::cout << "splicing" << std::endl;
    
    auto curr = leaf;
    auto prev = curr;
    
    bool is_single_branch {true};
    
    while (curr != root_) {
        prev = get_previous_allele(curr);
        is_single_branch = is_single_branch && boost::out_degree(curr, tree_) == 1;
        if (prev == root_ || overlaps(tree_[prev], region)) break;
        curr = prev;
    }
    
    auto splice_end = curr;
    
    //std::cout << "splice end: " << tree_[splice_end] << std::endl;
    //std::cout << "prev: " << tree_[prev] << std::endl;
    
    boost::remove_edge(prev, splice_end, tree_);
    
    curr = prev;
    
    while (curr != root_ && overlaps(region, tree_[curr])) {
        prev = get_previous_allele(curr);
        
        //std::cout << "curr: " << tree_[curr] << std::endl;
        //std::cout << "prev: " << tree_[prev] << std::endl;
        
        is_single_branch = is_single_branch && boost::out_degree(curr, tree_) == 0;
        
        if (is_single_branch) {
            //std::cout << "removing" << std::endl;
            boost::remove_edge(prev, curr, tree_);
            boost::remove_vertex(curr, tree_);
        }
        
        curr = prev;
    }
    
    //std::cout << "adding edge between " << tree_[curr] << " and " << tree_[splice_end] << std::endl;
    
    boost::add_edge(curr, splice_end, tree_);
    
    //exit(0);
    
    return std::make_pair(leaf, true);
}

} // namespace Octopus
