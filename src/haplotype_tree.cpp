//
//  haplotype_tree.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_tree.h"

#include <stdexcept>
#include <iterator>  // std::next, std::tie, std::cbegin etc
#include <algorithm> // std::find, std::find_if

#include "reference_genome.h"
#include "region_utils.h"

#include <iostream> // TEST

HaplotypeTree::HaplotypeTree(ReferenceGenome& the_reference)
:
the_reference_ {the_reference},
the_tree_ {},
the_root_ {boost::add_vertex(the_tree_)},
haplotype_leafs_ {the_root_},
haplotype_leaf_cache_ {},
recently_removed_haplotypes_ {}
{}

std::size_t HaplotypeTree::num_haplotypes() const
{
    return haplotype_leafs_.size();
}

void HaplotypeTree::extend(const Allele& an_allele)
{
    auto leaf_it = std::cbegin(haplotype_leafs_);
    
    while (leaf_it != std::cend(haplotype_leafs_)) {
        leaf_it = extend_haplotype(leaf_it, an_allele);
        ++leaf_it;
    }
    
    haplotype_leaf_cache_.clear();
    recently_removed_haplotypes_.clear();
}

HaplotypeTree::Haplotypes HaplotypeTree::get_haplotypes(const GenomicRegion& a_region)
{
    haplotype_leaf_cache_.clear();
    
    Haplotypes result {};
    
    for (auto leaf : haplotype_leafs_) {
        auto haplotype = get_haplotype(leaf, a_region);
        
        // recently retreived haplotypes are added to the cache as it is likely these
        // are the haplotypes that will be pruned next
        add_to_cache(haplotype, leaf);
        
        result.emplace_back(std::move(haplotype));
    }
    
    return result;
}

void HaplotypeTree::prune_all(const Haplotype& haplotype)
{
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
        //TODO: this is not quite right.. needs to delete all equal haplotypes, not just the ones
        // with the same alleles
        
        leaf_it = std::cbegin(haplotype_leafs_);
        
        while (true) {
            leaf_it = find_haplotype_leaf(leaf_it, std::cend(haplotype_leafs_), haplotype);
            
            if (leaf_it == std::cend(haplotype_leafs_)) {
                return;
            }
            
            std::tie(new_haplotype_end, new_end_is_leaf) = prune_branch(*leaf_it, haplotype.get_region());
            
            leaf_it = haplotype_leafs_.erase(leaf_it);
            
            if (new_end_is_leaf) {
                leaf_it = haplotype_leafs_.insert(leaf_it, new_haplotype_end);
            }
        }
    }
    
    recently_removed_haplotypes_.emplace(haplotype);
}

void HaplotypeTree::prune_unique(const Haplotype &haplotype)
{
    Vertex new_haplotype_end;
    Vertex leaf_to_save;
    LeafIterator leaf_it;
    bool new_end_is_leaf {};
    
    if (haplotype_leaf_cache_.count(haplotype) > 0) {
        auto leaf_range = haplotype_leaf_cache_.equal_range(haplotype);
        
        leaf_to_save = std::find_if(leaf_range.first, leaf_range.second,
                                    [this, &haplotype] (auto& leaf_pair) {
                                        return is_branch_the_haplotype(leaf_pair.second, haplotype);
                                    })->second;
        
        std::for_each(leaf_range.first, leaf_range.second,
          [this, &haplotype, &leaf_it, &new_haplotype_end, &new_end_is_leaf, &leaf_to_save] (auto& leaf_pair) {
              if (leaf_pair.second != leaf_to_save) {
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
        haplotype_leaf_cache_.emplace(haplotype, leaf_to_save);
    } else {
        // TODO
    }
}

HaplotypeTree::Vertex HaplotypeTree::get_previous_allele(Vertex allele) const
{
    return *boost::inv_adjacent_vertices(allele, the_tree_).first;
}

bool HaplotypeTree::allele_exists(Vertex allele, const Allele& an_allele) const
{
    auto vertex_range = boost::adjacent_vertices(allele, the_tree_);
    return std::any_of(vertex_range.first, vertex_range.second, [this, &an_allele] (Vertex vertex) {
        return the_tree_[vertex].the_allele == an_allele;
    });
}

HaplotypeTree::LeafIterator HaplotypeTree::extend_haplotype(LeafIterator haplotype_leaf,
                                                            const Allele& the_new_allele)
{
    const Allele& the_leaf_allele = the_tree_[*haplotype_leaf].the_allele;
    
    if (*haplotype_leaf == the_root_ || is_after(the_new_allele, the_leaf_allele)) {
        Vertex new_leaf = boost::add_vertex(the_tree_);
        the_tree_[new_leaf].the_allele = the_new_allele;
        boost::add_edge(*haplotype_leaf, new_leaf, the_tree_);
        
        haplotype_leaf = haplotype_leafs_.erase(haplotype_leaf);
        haplotype_leaf = haplotype_leafs_.insert(haplotype_leaf, new_leaf);
    } else if (overlaps(the_new_allele, the_leaf_allele)) {
        Vertex previous_allele = get_previous_allele(*haplotype_leaf);
        
        if (!allele_exists(previous_allele, the_new_allele)) {
            Vertex new_branch = boost::add_vertex(the_tree_);
            the_tree_[new_branch].the_allele = the_new_allele;
            boost::add_edge(previous_allele, new_branch, the_tree_);
            
            haplotype_leafs_.insert(haplotype_leaf, new_branch);
        }
    }
    
    return haplotype_leaf;
}

Haplotype HaplotypeTree::get_haplotype(Vertex haplotype_leaf, const GenomicRegion& a_region)
{
    Haplotype result {the_reference_, a_region};
    
    while (haplotype_leaf != the_root_ && is_after(the_tree_[haplotype_leaf].the_allele, a_region)) {
        haplotype_leaf = get_previous_allele(haplotype_leaf);
    }
    
    while (haplotype_leaf != the_root_ && overlaps(the_tree_[haplotype_leaf].the_allele, a_region)) {
        result.push_front(the_tree_[haplotype_leaf].the_allele);
        haplotype_leaf = get_previous_allele(haplotype_leaf);
    }
    
    return result;
}

bool HaplotypeTree::is_branch_the_haplotype(Vertex this_vertex, const Haplotype& haplotype) const
{
    if (this_vertex == the_root_) {
        return false;
    }
    
    const Allele* this_allele {&the_tree_[this_vertex].the_allele}; // pointers are just for efficiency
    
    if (!overlaps(haplotype, *this_allele)) {
        return false;
    }
    
    Vertex previous_vertex;
    const Allele* previous_allele;
    
    while (!is_before(*this_allele, haplotype)) {
        if (!haplotype.contains(*this_allele)) {
            return false;
        }
        
        previous_vertex = get_previous_allele(this_vertex);
        
        if (previous_vertex == the_root_) {
            break;
        }
        
        previous_allele = &the_tree_[previous_vertex].the_allele;
        
        if (is_before(*previous_allele, haplotype)) {
            break;
        }
        
        if (!are_adjacent(*previous_allele, *this_allele) &&
            !haplotype.contains(get_reference_allele(get_intervening_region(*previous_allele, *this_allele),
                                                     the_reference_))) {
            return false;
        }
        
        this_vertex = previous_vertex;
        this_allele = previous_allele;
    }
    
    return true;
}

HaplotypeTree::LeafIterator HaplotypeTree::find_haplotype_leaf(LeafIterator first, LeafIterator last,
                                                               const Haplotype& haplotype) const
{
    return std::find_if(first, last, [this, &haplotype] (Vertex leaf) {
        return is_branch_the_haplotype(leaf, haplotype);
    });
}

std::pair<HaplotypeTree::Vertex, bool> HaplotypeTree::prune_branch(Vertex leaf, const GenomicRegion& a_region)
{
    Vertex new_leaf;
    
    while (leaf != the_root_) {
        if (boost::out_degree(leaf, the_tree_) > 0) {
            return {leaf, false};
        } else if (begins_before(the_tree_[leaf].the_allele, a_region)) {
            return {leaf, true};
        }
        
        new_leaf = get_previous_allele(leaf);
        boost::remove_edge(new_leaf, leaf, the_tree_);
        boost::remove_vertex(leaf, the_tree_);
        leaf = new_leaf;
    }
    
    // the root should only be indicated as a leaf node if there are no other nodes in the tree
    return {leaf, num_haplotypes() == 1};
}

void HaplotypeTree::add_to_cache(const Haplotype& haplotype, Vertex leaf)
{
    haplotype_leaf_cache_.emplace(haplotype, leaf);
}

std::pair<HaplotypeTree::LeafIterator, bool> HaplotypeTree::get_leaf_from_cache(const Haplotype& haplotype)
{
    if (haplotype_leaf_cache_.count(haplotype) == 0) {
        return {std::cend(haplotype_leafs_), false};
    }
    
    auto possible_leaf_range = haplotype_leaf_cache_.equal_range(haplotype);
    
    while (possible_leaf_range.first != possible_leaf_range.second) {
        if (is_branch_the_haplotype(possible_leaf_range.first->second, haplotype)) {
            haplotype_leaf_cache_.erase(possible_leaf_range.first);
            
            return {std::find(std::cbegin(haplotype_leafs_), std::cend(haplotype_leafs_),
                              possible_leaf_range.first->second), true};
        }
        
        ++possible_leaf_range.first;
    }
    
    return {std::cend(haplotype_leafs_), false};
}
