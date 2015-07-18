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
#include <algorithm> // std::find, std::find_if, std::any_of

#include "genomic_region.h"
#include "reference_genome.h"
#include "region_algorithms.h"

namespace Octopus
{

HaplotypeTree::HaplotypeTree(ReferenceGenome& the_reference)
:
the_reference_ {the_reference},
the_tree_ {},
the_root_ {boost::add_vertex(the_tree_)},
the_haplotype_leafs_ {the_root_},
haplotype_leaf_cache_ {},
recently_removed_haplotypes_ {}
{}

bool HaplotypeTree::empty() const noexcept
{
    return the_haplotype_leafs_.front() == the_root_;
}

std::size_t HaplotypeTree::num_haplotypes() const noexcept
{
    return the_haplotype_leafs_.size();
}

bool HaplotypeTree::contains(const Haplotype& a_haplotype) const
{
    if (haplotype_leaf_cache_.count(a_haplotype) > 0) return true;
    
    return std::any_of(std::cbegin(the_haplotype_leafs_), std::cend(the_haplotype_leafs_),
                       [this, &a_haplotype] (const Vertex leaf) {
                           return is_branch_equal_haplotype(leaf, a_haplotype);
                       });
}

bool HaplotypeTree::is_unique(const Haplotype& a_haplotype) const
{
    if (haplotype_leaf_cache_.count(a_haplotype) > 0) {
        return haplotype_leaf_cache_.count(a_haplotype) == 1;
    } else {
        bool haplotype_seen {false};
        
        for (const Vertex& leaf : the_haplotype_leafs_) {
            if (is_branch_equal_haplotype(leaf, a_haplotype)) {
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

GenomicRegion HaplotypeTree::get_seperation_region(const Haplotype& the_first_haplotype,
                                                   const Haplotype& the_second_haplotype) const
{
    if (the_first_haplotype == the_second_haplotype) {
        return the_first_haplotype.get_region();
    }
    
    if (!(is_unique(the_first_haplotype) && is_unique(the_second_haplotype))) {
        throw std::runtime_error {"cannot find seperation region in non-unique haplotypes"};
    }
    
    Vertex first_allele, second_allele;
    
    if (haplotype_leaf_cache_.count(the_first_haplotype) > 0) {
        first_allele = haplotype_leaf_cache_.find(the_first_haplotype)->second;
    } else {
        first_allele = *find_equal_haplotype_leaf(std::cbegin(the_haplotype_leafs_),
                                                  std::cend(the_haplotype_leafs_), the_first_haplotype);
    }
    if (haplotype_leaf_cache_.count(the_second_haplotype) > 0) {
        second_allele = haplotype_leaf_cache_.find(the_second_haplotype)->second;
    } else {
        second_allele = *find_equal_haplotype_leaf(std::cbegin(the_haplotype_leafs_),
                                                   std::cend(the_haplotype_leafs_), the_second_haplotype);
    }
    
    while (first_allele != the_root_ && second_allele != the_root_) {
        if (first_allele == second_allele) return the_tree_[first_allele].get_region();
        
        first_allele  = get_previous_allele(first_allele);
        second_allele = get_previous_allele(second_allele);
    }
    
    return the_first_haplotype.get_region();
}

void HaplotypeTree::extend(const Allele& an_allele)
{
    auto leaf_it = std::cbegin(the_haplotype_leafs_);
    
    while (leaf_it != std::cend(the_haplotype_leafs_)) {
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
    result.reserve(the_haplotype_leafs_.size());
    
    for (auto leaf : the_haplotype_leafs_) {
        auto haplotype = get_haplotype(leaf, a_region);
        
        // recently retreived haplotypes are added to the cache as it is likely these
        // are the haplotypes that will be pruned next
        haplotype_leaf_cache_.emplace(haplotype, leaf);
        
        result.emplace_back(std::move(haplotype));
    }
    
    return result;
}

void HaplotypeTree::prune_all(const Haplotype& a_haplotype)
{
    if (recently_removed_haplotypes_.count(a_haplotype) > 1) return;
    
    Vertex new_haplotype_end;
    LeafIterator leaf_it;
    bool new_end_is_leaf {};
    
    // if any of the haplotypes in cache match the query haplotype then the cache must contain
    // all possible leaves corrosponding to that haplotype. So we don't need to look through
    // the list of all leaves. Win.
    if (haplotype_leaf_cache_.count(a_haplotype) > 0) {
        auto leaf_range = haplotype_leaf_cache_.equal_range(a_haplotype);
        
        std::for_each(leaf_range.first, leaf_range.second,
              [this, &a_haplotype, &leaf_it, &new_haplotype_end, &new_end_is_leaf] (auto& leaf_pair) {
                  std::tie(new_haplotype_end, new_end_is_leaf) = prune_branch(leaf_pair.second,
                                                                              a_haplotype.get_region());
                  
                  leaf_it = std::find(std::cbegin(the_haplotype_leafs_), std::cend(the_haplotype_leafs_),
                                      leaf_pair.second);
                  
                  leaf_it = the_haplotype_leafs_.erase(leaf_it);
                  
                  if (new_end_is_leaf) {
                      the_haplotype_leafs_.insert(leaf_it, new_haplotype_end);
                  }
              });
        
        haplotype_leaf_cache_.erase(a_haplotype);
    } else {
        leaf_it = std::cbegin(the_haplotype_leafs_);
        
        while (true) {
            leaf_it = find_equal_haplotype_leaf(leaf_it, std::cend(the_haplotype_leafs_), a_haplotype);
            
            if (leaf_it == std::cend(the_haplotype_leafs_)) {
                return;
            }
            
            std::tie(new_haplotype_end, new_end_is_leaf) = prune_branch(*leaf_it, a_haplotype.get_region());
            
            leaf_it = the_haplotype_leafs_.erase(leaf_it);
            
            if (new_end_is_leaf) {
                leaf_it = the_haplotype_leafs_.insert(leaf_it, new_haplotype_end);
            }
        }
    }
    
    recently_removed_haplotypes_.emplace(a_haplotype);
}

void HaplotypeTree::prune_unique(const Haplotype& a_haplotype)
{
    Vertex new_haplotype_end;
    Vertex leaf_to_keep;
    LeafIterator leaf_it;
    bool new_end_is_leaf {};
    
    if (haplotype_leaf_cache_.count(a_haplotype) > 0) {
        auto leaf_range = haplotype_leaf_cache_.equal_range(a_haplotype);
        
        leaf_to_keep = std::find_if(leaf_range.first, leaf_range.second,
                                    [this, &a_haplotype] (auto& leaf_pair) {
                                        return is_branch_exact_haplotype(leaf_pair.second, a_haplotype);
                                    })->second;
        
        std::for_each(leaf_range.first, leaf_range.second,
          [this, &a_haplotype, &leaf_it, &new_haplotype_end, &new_end_is_leaf, &leaf_to_keep] (auto& leaf_pair) {
              if (leaf_pair.second != leaf_to_keep) {
                  std::tie(new_haplotype_end, new_end_is_leaf) = prune_branch(leaf_pair.second,
                                                                              a_haplotype.get_region());
                  
                  leaf_it = std::find(std::cbegin(the_haplotype_leafs_), std::cend(the_haplotype_leafs_),
                                      leaf_pair.second);
                  
                  leaf_it = the_haplotype_leafs_.erase(leaf_it);
                  
                  if (new_end_is_leaf) {
                      the_haplotype_leafs_.insert(leaf_it, new_haplotype_end);
                  }
              }
          });
        
        haplotype_leaf_cache_.erase(a_haplotype);
        haplotype_leaf_cache_.emplace(a_haplotype, leaf_to_keep);
    } else {
        leaf_to_keep = *find_exact_haplotype_leaf(std::cbegin(the_haplotype_leafs_), std::cend(the_haplotype_leafs_),
                                                  a_haplotype);
        
        leaf_it = std::cbegin(the_haplotype_leafs_);
        
        while (true) {
            leaf_it = find_equal_haplotype_leaf(leaf_it, std::cend(the_haplotype_leafs_), a_haplotype);
            
            if (*leaf_it == leaf_to_keep) continue;
            
            if (leaf_it == std::cend(the_haplotype_leafs_)) {
                return;
            }
            
            std::tie(new_haplotype_end, new_end_is_leaf) = prune_branch(*leaf_it, a_haplotype.get_region());
            
            leaf_it = the_haplotype_leafs_.erase(leaf_it);
            
            if (new_end_is_leaf) {
                leaf_it = the_haplotype_leafs_.insert(leaf_it, new_haplotype_end);
            }
        }

    }
}

void HaplotypeTree::clear()
{
    haplotype_leaf_cache_.clear();
    recently_removed_haplotypes_.clear();
    the_haplotype_leafs_.clear();
    the_tree_.clear();
    the_root_ = boost::add_vertex(the_tree_);
    the_haplotype_leafs_.push_back(the_root_);
}

// Private methods

HaplotypeTree::Vertex HaplotypeTree::get_previous_allele(Vertex an_allele) const
{
    return *boost::inv_adjacent_vertices(an_allele, the_tree_).first;
}

bool HaplotypeTree::allele_exists(Vertex allele, const Allele& an_allele) const
{
    auto vertex_range = boost::adjacent_vertices(allele, the_tree_);
    return std::any_of(vertex_range.first, vertex_range.second, [this, &an_allele] (Vertex vertex) {
        return the_tree_[vertex] == an_allele;
    });
}

HaplotypeTree::LeafIterator HaplotypeTree::extend_haplotype(LeafIterator haplotype_leaf,
                                                            const Allele& the_new_allele)
{
    const Allele& the_leaf_allele = the_tree_[*haplotype_leaf];
    
    if (are_adjacent(the_leaf_allele, the_new_allele) &&
        ((is_insertion(the_leaf_allele) && is_deletion(the_new_allele)) ||
        (is_deletion(the_leaf_allele) && is_insertion(the_new_allele)))) {
        return haplotype_leaf;
    }
    
    if (*haplotype_leaf == the_root_ || is_after(the_new_allele, the_leaf_allele)) {
        Vertex new_leaf = boost::add_vertex(the_tree_);
        the_tree_[new_leaf] = the_new_allele;
        boost::add_edge(*haplotype_leaf, new_leaf, the_tree_);
        
        haplotype_leaf = the_haplotype_leafs_.erase(haplotype_leaf);
        haplotype_leaf = the_haplotype_leafs_.insert(haplotype_leaf, new_leaf);
    } else if (overlaps(the_new_allele, the_leaf_allele)) {
        Vertex previous_allele = get_previous_allele(*haplotype_leaf);
        
        const Allele& the_previous_allele = the_tree_[previous_allele];
        
        if (are_adjacent(the_previous_allele, the_new_allele) &&
            ((is_insertion(the_previous_allele) && is_deletion(the_new_allele)) ||
             (is_deletion(the_previous_allele) && is_insertion(the_new_allele)))) {
            return haplotype_leaf;
        }
        
        if (!allele_exists(previous_allele, the_new_allele)) {
            Vertex new_branch = boost::add_vertex(the_tree_);
            the_tree_[new_branch] = the_new_allele;
            boost::add_edge(previous_allele, new_branch, the_tree_);
            
            the_haplotype_leafs_.insert(haplotype_leaf, new_branch);
        }
    }
    
    return haplotype_leaf;
}

Haplotype HaplotypeTree::get_haplotype(Vertex haplotype_leaf, const GenomicRegion& a_region) const
{
    Haplotype result {the_reference_, a_region};
    
    while (haplotype_leaf != the_root_ && is_after(the_tree_[haplotype_leaf], a_region)) {
        haplotype_leaf = get_previous_allele(haplotype_leaf);
    }
    
    while (haplotype_leaf != the_root_ && overlaps(the_tree_[haplotype_leaf], a_region)) {
        result.push_front(the_tree_[haplotype_leaf]);
        haplotype_leaf = get_previous_allele(haplotype_leaf);
    }
    
    return result;
}

bool HaplotypeTree::is_branch_exact_haplotype(Vertex this_vertex, const Haplotype& a_haplotype) const
{
    if (this_vertex == the_root_) {
        return false;
    }
    
    const Allele* this_allele {&the_tree_[this_vertex]};
    
    if (!overlaps(a_haplotype, *this_allele)) {
        return false;
    }
    
    Vertex previous_vertex;
    const Allele* previous_allele;
    
    while (!is_before(*this_allele, a_haplotype)) {
        if (!a_haplotype.contains(*this_allele)) {
            return false;
        }
        
        previous_vertex = get_previous_allele(this_vertex);
        
        if (previous_vertex == the_root_) {
            break;
        }
        
        previous_allele = &the_tree_[previous_vertex];
        
        if (is_before(*previous_allele, a_haplotype)) {
            break;
        }
        
        if (!are_adjacent(*previous_allele, *this_allele) &&
            !a_haplotype.contains(get_reference_allele(get_intervening(*previous_allele, *this_allele),
                                                       the_reference_))) {
            return false;
        }
        
        this_vertex = previous_vertex;
        this_allele = previous_allele;
    }
    
    return true;
}

bool HaplotypeTree::is_branch_equal_haplotype(Vertex this_vertex, const Haplotype& a_haplotype) const
{
    if (this_vertex == the_root_) {
        return false;
    }
    
    if (!overlaps(a_haplotype, the_tree_[this_vertex])) {
        return false;
    }
    
    return get_haplotype(this_vertex, a_haplotype.get_region()) == a_haplotype;
}

HaplotypeTree::LeafIterator HaplotypeTree::find_exact_haplotype_leaf(LeafIterator first, LeafIterator last,
                                                                     const Haplotype& a_haplotype) const
{
    return std::find_if(first, last, [this, &a_haplotype] (Vertex leaf) {
        return is_branch_exact_haplotype(leaf, a_haplotype);
    });
}

HaplotypeTree::LeafIterator HaplotypeTree::find_equal_haplotype_leaf(LeafIterator first, LeafIterator last,
                                                                     const Haplotype& a_haplotype) const
{
    return std::find_if(first, last, [this, &a_haplotype] (Vertex leaf) {
        return is_branch_equal_haplotype(leaf, a_haplotype);
    });
}

std::pair<HaplotypeTree::Vertex, bool> HaplotypeTree::prune_branch(Vertex leaf, const GenomicRegion& a_region)
{
    Vertex new_leaf;
    
    while (leaf != the_root_) {
        if (boost::out_degree(leaf, the_tree_) > 0) {
            return {leaf, false};
        } else if (begins_before(the_tree_[leaf], a_region)) {
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

} // end namespace Octopus
